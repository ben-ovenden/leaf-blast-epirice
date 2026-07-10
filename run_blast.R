#!/usr/bin/env Rscript
################################################################################
# run_blast.R
#
# Weekly rice leaf blast risk using the EPIRICE model (Savary et al. 2012),
# driven by Open-Meteo weather. Run from git / GitHub Actions.
#
#   Rscript run_blast.R
#
# Pipeline:
#   1. For each site, fetch daily weather from emergence to the latest archive
#      day (Open-Meteo ERA5).
#   2. Run the EPIRICE leaf blast SEIR up to that day.
#   3. Report current modelled intensity, its recent trend, and a risk band.
#   4. Write a CSV, a text summary, and a PNG map; optionally email them.
################################################################################

SCRIPT_DIR <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile), winslash = "/"),
  error = function(e) normalizePath(getwd(), winslash = "/")
)

source(file.path(SCRIPT_DIR, "epirice_model.R"))
source(file.path(SCRIPT_DIR, "openmeteo_wth.R"))
source(file.path(SCRIPT_DIR, "blast_config.R"))

# A weather-fetch hook, so tests can inject synthetic weather. Defaults to the
# live Open-Meteo adapter.
if (!exists("WTH_FN")) WTH_FN <- get_openmeteo_wth

classify <- function(intensity) {
  if (is.na(intensity)) return("no data")
  if (intensity < INTENSITY_LOW_MAX) return("low")
  if (intensity < INTENSITY_MODERATE_MAX) return("moderate")
  "high"
}

run_site <- function(name, lat, lon, emergence, end_date) {
  na_row <- data.table(name = name, lat = lat, lon = lon,
                       intensity = NA_real_, peak = NA_real_,
                       trend7 = NA_real_, level = "no data",
                       last_date = NA_character_, days = 0L)

  if (end_date <= as.Date(emergence) + MIN_DAYS) {
    na_row[, level := "pre-season"]
    return(na_row)
  }

  wth <- tryCatch(
    WTH_FN(lat = lat, lon = lon, start_date = emergence, end_date = end_date),
    error = function(e) NULL
  )
  if (is.null(wth) || nrow(wth) < MIN_DAYS) return(na_row)

  duration <- as.integer(min(120L, nrow(wth)))

  lb <- tryCatch(
    predict_leaf_blast(wth, emergence = emergence, duration = duration),
    error = function(e) NULL
  )
  if (is.null(lb) || nrow(lb) == 0) return(na_row)

  cur  <- lb$intensity[nrow(lb)]
  peak <- max(lb$intensity, na.rm = TRUE)
  trend7 <- if (nrow(lb) > 7) cur - lb$intensity[nrow(lb) - 7] else NA_real_

  data.table(name = name, lat = lat, lon = lon,
             intensity = cur, peak = peak, trend7 = trend7,
             level = classify(cur),
             last_date = as.character(lb$dates[nrow(lb)]),
             days = duration)
}

# ---- Run all sites --------------------------------------------------------
cat("\nRice leaf blast risk (EPIRICE / Open-Meteo)\n")
cat(strrep("=", 60), "\n")

end_date <- Sys.Date() - ARCHIVE_LAG_DAYS
# Use the same rolling window as the heatmap: assume a crop of age CROP_AGE_DAYS
# everywhere, so the table reflects current potential risk rather than a finished
# past season. Per-site emergence overrides are still honoured if provided.
rolling_emergence <- as.character(end_date - CROP_AGE_DAYS)
sites <- as.data.table(SITES)
if (!"emergence" %in% names(sites)) sites[, emergence := rolling_emergence]

cat("Run date:   ", format(Sys.Date(), "%A %d %B %Y"), "\n")
cat("Weather to: ", format(end_date, "%Y-%m-%d"), " (archive lag ",
    ARCHIVE_LAG_DAYS, " days)\n", sep = "")
cat("Crop age:   ", CROP_AGE_DAYS, " days (emergence ", rolling_emergence, ")\n\n",
    sep = "")

results <- rbindlist(lapply(seq_len(nrow(sites)), function(k) {
  s <- sites[k]
  cat(sprintf("%-14s ... ", s$name))
  r <- run_site(s$name, s$lat, s$lon, s$emergence, end_date)
  msg <- if (is.na(r$intensity)) r$level else
    sprintf("%.1f%% (%s)", r$intensity * 100, r$level)
  cat(msg, "\n")
  r
}))

# ---- Summary --------------------------------------------------------------
counts <- results[, .N, by = level]
getn <- function(lv) { x <- counts[level == lv, N]; if (length(x)) x else 0L }

summary_lines <- c(
  "Rice leaf blast risk summary",
  strrep("=", 60),
  paste0("Generated:  ", format(Sys.Date(), "%A %d %B %Y")),
  paste0("Model:      EPIRICE leaf blast (Savary et al. 2012)"),
  paste0("Weather:    Open-Meteo ERA5 archive, to ", format(end_date, "%Y-%m-%d")),
  "",
  sprintf("High:      %d", getn("high")),
  sprintf("Moderate:  %d", getn("moderate")),
  sprintf("Low:       %d", getn("low")),
  sprintf("Pre-season:%d", getn("pre-season")),
  sprintf("No data:   %d", getn("no data")),
  "",
  "Site detail (current modelled leaf blast intensity):",
  strrep("-", 60)
)
for (k in seq_len(nrow(results))) {
  r <- results[k]
  pct <- if (is.na(r$intensity)) "  -  " else sprintf("%5.1f%%", r$intensity * 100)
  tr  <- if (is.na(r$trend7)) "" else sprintf("  7d %+0.1f pts", r$trend7 * 100)
  summary_lines <- c(summary_lines,
    sprintf("%-14s %s  %-10s %s%s", r$name, pct, r$level,
            ifelse(is.na(r$last_date), "", r$last_date), tr))
}
summary_lines <- c(summary_lines, strrep("-", 60),
  "Intensity is the EPIRICE proportion of diseased sites. Risk band cut points",
  "are provisional; calibrate against field observations before relying on them.",
  "",
  if (exists("CITATION")) CITATION else NULL)

summary_text <- paste(summary_lines, collapse = "\n")
cat("\n", summary_text, "\n", sep = "")

# ---- Write outputs --------------------------------------------------------
OUT <- file.path(SCRIPT_DIR, OUTPUT_DIR)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
run_tag <- format(Sys.Date(), "%Y-%m-%d")

fwrite(results, file.path(OUT, sprintf("blast_results_%s.csv", run_tag)))
writeLines(summary_text, file.path(OUT, sprintf("blast_summary_%s.txt", run_tag)))
# Stable "latest" copies so automation can reference a fixed filename
fwrite(results, file.path(OUT, "blast_results_latest.csv"))
writeLines(summary_text, file.path(OUT, "blast_summary_latest.txt"))

band_col <- function(level) {
  vapply(level, function(l) switch(l,
    low = COL_LOW, moderate = COL_MODERATE, high = COL_HIGH, COL_NODATA),
    character(1))
}

if (isTRUE(MAKE_MAP)) {
  map_file <- file.path(OUT, sprintf("blast_map_%s.png", run_tag))
  png(map_file, width = MAP_WIDTH, height = MAP_HEIGHT, res = 120)
  op <- par(mar = c(4, 4, 3, 1))
  xr <- range(results$lon) + c(-2, 2)
  yr <- range(results$lat) + c(-2, 2)
  plot(results$lon, results$lat, type = "n", xlim = xr, ylim = yr,
       xlab = "Longitude", ylab = "Latitude",
       main = sprintf("Rice leaf blast risk  %s", run_tag))
  grid(col = "#E4E7E9")
  points(results$lon, results$lat, pch = 21, cex = 3,
         bg = band_col(results$level), col = "#22272B", lwd = 1.5)
  text(results$lon, results$lat, labels = results$name, pos = 3,
       offset = 1, cex = 0.8, col = "#22272B")
  legend("topright", legend = c("Low", "Moderate", "High", "No data / pre-season"),
         pt.bg = c(COL_LOW, COL_MODERATE, COL_HIGH, COL_NODATA),
         pch = 21, pt.cex = 2, bty = "n")
  par(op); dev.off()
  file.copy(map_file, file.path(OUT, "blast_map_latest.png"), overwrite = TRUE)
  cat("\nMap:     ", map_file, "\n")
}

cat("Results: ", file.path(OUT, sprintf("blast_results_%s.csv", run_tag)), "\n")

# ---- Optional email -------------------------------------------------------
if (isTRUE(SEND_EMAIL) && nzchar(EMAIL_FROM) && nzchar(EMAIL_PASSWORD) &&
    nzchar(EMAIL_TO)) {
  if (!requireNamespace("mailR", quietly = TRUE)) {
    cat("\nEmail skipped: install.packages('mailR') to enable.\n")
  } else {
    att <- c(file.path(OUT, sprintf("blast_summary_%s.txt", run_tag)),
             file.path(OUT, sprintf("blast_results_%s.csv", run_tag)))
    if (isTRUE(MAKE_MAP)) att <- c(att, file.path(OUT, sprintf("blast_map_%s.png", run_tag)))
    tryCatch({
      mailR::send.mail(
        from = EMAIL_FROM, to = strsplit(EMAIL_TO, ",")[[1]],
        subject = sprintf("Leaf blast risk %s", run_tag),
        body = summary_text,
        smtp = list(host.name = EMAIL_SMTP_HOST, port = EMAIL_SMTP_PORT,
                    user.name = EMAIL_FROM, passwd = EMAIL_PASSWORD, ssl = TRUE),
        authenticate = TRUE, send = TRUE, attach.files = att)
      cat("\nEmail sent to", EMAIL_TO, "\n")
    }, error = function(e) cat("\nEmail failed:", conditionMessage(e), "\n"))
  }
}

cat("\nDone.\n")

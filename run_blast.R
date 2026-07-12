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
source(file.path(SCRIPT_DIR, "blastam_model.R"))
source(file.path(SCRIPT_DIR, "openmeteo_wth.R"))
source(file.path(SCRIPT_DIR, "blast_config.R"))

# Weather-fetch hook (tests inject synthetic data). Fetches hourly and returns
# daily rows carrying both EPIRICE inputs and the BLASTAM night judgement.
if (!exists("WTH_FN")) WTH_FN <- function(lat, lon, start_date, end_date) {
  hw <- get_openmeteo_hourly(lat, lon, start_date, end_date)
  if (is.null(hw) || nrow(hw) == 0) return(NULL)
  blastam_daily_from_hourly(hw)
}

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
                       last_date = NA_character_, days = 0L,
                       blast_events = NA_integer_, blast_recent = NA_integer_)

  if (end_date <= as.Date(emergence) + MIN_DAYS) {
    na_row[, level := "pre-season"]
    return(na_row)
  }

  dd <- tryCatch(
    WTH_FN(lat = lat, lon = lon, start_date = emergence, end_date = end_date),
    error = function(e) NULL
  )
  if (is.null(dd) || nrow(dd) < MIN_DAYS) return(na_row)

  # BLASTAM: infection-favoured days (total over window, and last 7 days)
  bs <- blastam_score(dd$infect, dd$semi, as.Date(dd$date), end_date,
                      window = if (exists("BLASTAM_WINDOW_DAYS")) BLASTAM_WINDOW_DAYS else 21L,
                      recent = 7L)

  # EPIRICE: build daily weather and run the SEIR
  wth <- data.table(YYYYMMDD = as.Date(dd$date),
                    DOY = as.integer(format(as.Date(dd$date), "%j")),
                    TEMP = dd$TEMP, RHUM = dd$RHUM, RAIN = dd$RAIN, LAT = lat, LON = lon)
  setorder(wth, YYYYMMDD)
  duration <- as.integer(min(120L, nrow(wth)))
  lb <- tryCatch(
    predict_leaf_blast(wth, emergence = emergence, duration = duration),
    error = function(e) NULL
  )
  if (is.null(lb) || nrow(lb) == 0) {
    na_row[, `:=`(blast_events = bs$events, blast_recent = bs$recent)]
    return(na_row)
  }

  cur  <- lb$intensity[nrow(lb)]
  peak <- max(lb$intensity, na.rm = TRUE)
  trend7 <- if (nrow(lb) > 7) cur - lb$intensity[nrow(lb) - 7] else NA_real_

  data.table(name = name, lat = lat, lon = lon,
             intensity = cur, peak = peak, trend7 = trend7,
             level = classify(cur),
             last_date = as.character(lb$dates[nrow(lb)]),
             days = duration,
             blast_events = bs$events, blast_recent = bs$recent)
}

# ---- Run all sites --------------------------------------------------------
cat("\nRice leaf blast risk (EPIRICE / Open-Meteo)\n")
cat(strrep("=", 60), "\n")

end_date <- Sys.Date() - ARCHIVE_LAG_DAYS
# Use the same rolling window as the heatmap: assume a crop of age CROP_AGE_DAYS
# everywhere, so the table reflects current potential risk rather than a finished
# past season. Per-site emergence overrides are still honoured if provided.
rolling_emergence <- as.character(end_date - CROP_AGE_DAYS)
sites <- as.data.table(MONITOR_TOWNS)
if (!"emergence" %in% names(sites)) sites[, emergence := rolling_emergence]

cat("Run date:   ", format(Sys.Date(), "%A %d %B %Y"), "\n")
cat("Weather to: ", format(end_date, "%Y-%m-%d"), " (archive lag ",
    ARCHIVE_LAG_DAYS, " days)\n", sep = "")
cat("Crop age:   ", CROP_AGE_DAYS, " days (emergence ", rolling_emergence, ")\n\n",
    sep = "")

options(timeout = max(120, getOption("timeout", 60)))  # slow archive requests

# Towns are independent, and each is a separate (slow, cold) hourly request, so
# fetch and model them concurrently. mclapply forks on Linux (the GitHub runner);
# on Windows it falls back to serial automatically. TOWN_FETCH_CORES caps the
# concurrency; 31 town requests * ~6.4 weighted calls stays well under the rate
# limit even fired together.
.cores <- if (exists("TOWN_FETCH_CORES")) TOWN_FETCH_CORES else 8L
one_site <- function(k) {
  s <- sites[k]
  run_site(s$name, s$lat, s$lon, s$emergence, end_date)
}
results_list <- tryCatch(
  parallel::mclapply(seq_len(nrow(sites)), one_site, mc.cores = .cores),
  error = function(e) lapply(seq_len(nrow(sites)), one_site))
# Serial retry (clean, non-forked connection) for towns whose concurrent fetch
# was dropped by the API ("no data") or whose fork failed outright.
needs_retry <- function(r) (!is.data.frame(r)) ||
  identical(as.character(r$level), "no data")
retry <- which(vapply(results_list, needs_retry, logical(1)))
if (length(retry) > 0) {
  cat(sprintf("Serial retry for %d town(s) that returned no data...\n", length(retry)))
  for (k in retry) results_list[[k]] <- one_site(k)
}
results <- rbindlist(results_list)
for (k in seq_len(nrow(results))) {
  r <- results[k]
  msg <- if (is.na(r$intensity)) r$level else sprintf("%.1f%% (%s)", r$intensity * 100, r$level)
  cat(sprintf("%-14s ... %s\n", r$name, msg))
}

# Attach state and order by state then town (alphabetical within state)
st <- as.data.table(MONITOR_TOWNS)[, .(name, state)]
results <- merge(results, st, by = "name", all.x = TRUE, sort = FALSE)
results[is.na(state), state := "--"]
setorder(results, state, name)

# ---- Summary --------------------------------------------------------------
counts <- results[, .N, by = level]
getn <- function(lv) { x <- counts[level == lv, N]; if (length(x)) x else 0L }

# Read the grid's map stats (written by run_blast_grid.R) into a one-line note so
# the reader can see the cache growing and the map sharpening week to week.
map_growth_line <- function() {
  f <- file.path(SCRIPT_DIR, OUTPUT_DIR, "map_stats.txt")
  if (!file.exists(f)) return(NULL)
  s <- tryCatch(strsplit(readLines(f, warn = FALSE)[1], "\\|")[[1]],
                error = function(e) NULL)
  if (length(s) < 4) return(NULL)
  now <- as.integer(s[1]); sp <- as.numeric(s[2])
  prev <- as.integer(s[3]); finest <- as.numeric(s[4])
  chg <- if (is.na(prev) || prev <= 0) ""
         else if (now > prev) sprintf(" (up from %d last week)", prev)
         else " (steady)"
  at_target <- !is.na(sp) && !is.na(finest) && sp <= finest * 1.05
  tail <- if (at_target) "; at target resolution"
          else sprintf("; sharpening toward ~%.2f deg", finest)
  sprintf("%d grid points at ~%.2f deg spacing%s%s", now, sp, chg, tail)
}
mg <- map_growth_line()

summary_lines <- c(
  "Blast risk summary",
  strrep("=", 60),
  paste0("Generated:  ", format(Sys.Date(), "%A %d %B %Y")),
  paste0("Models:     EPIRICE (Savary et al. 2012) + BLASTAM (Koshimizu 1988)"),
  paste0("Weather:    Open-Meteo ERA5 archive, to ", format(end_date, "%Y-%m-%d")),
  if (!is.null(mg)) paste0("Map:        ", mg) else NULL,
  "",
  sprintf("EPIRICE bands  High:%d  Moderate:%d  Low:%d  Pre-season:%d  No data:%d",
          getn("high"), getn("moderate"), getn("low"),
          getn("pre-season"), getn("no data")),
  "",
  sprintf("Town detail: EPIRICE intensity + BLASTAM favourable days (last %d), by state:", BLASTAM_WINDOW_DAYS),
  strrep("-", 60)
)
cur_state <- ""
for (k in seq_len(nrow(results))) {
  r <- results[k]
  if (r$state != cur_state) {
    cur_state <- r$state
    summary_lines <- c(summary_lines, "", paste0("[", cur_state, "]"))
  }
  pct <- if (is.na(r$intensity)) "  -  " else sprintf("%5.1f%%", r$intensity * 100)
  bl  <- if (is.na(r$blast_events)) "-" else
    sprintf("%d days (7d %d)", r$blast_events, r$blast_recent)
  summary_lines <- c(summary_lines,
    sprintf("  %-14s EPIRICE %s %-9s  BLASTAM %s", r$name, pct, r$level, bl))
}
summary_lines <- c(summary_lines, "",
  "About these estimates: two different models",
  strrep("-", 60),
  "EPIRICE (intensity %) mechanistically simulates the whole epidemic. It steps",
  sprintf("through a crop about %d days old, tracking healthy, latent, infectious", CROP_AGE_DAYS),
  "and removed leaf sites, and reports the proportion of leaf tissue diseased",
  "(0-100%). It answers: how much disease has the season built up? It is slow and",
  "cumulative, and reads near zero in cool, dry conditions.",
  "",
  "BLASTAM (infection days) is the Japanese infection-warning model (Koshimizu",
  "1988). For each day it judges whether conditions favoured a NEW infection, and",
  "counts the favourable days in the last 21 days. A day is favourable when all three",
  "hold: leaf wetness >=10 h; mean temperature during wetness 15-32 C; and the",
  "preceding 5-day mean temperature 20-30 C (upper bounds raised from the",
  "Japanese 25 C for warm conditions). (If wetness >=10 h but one",
  "temperature is out of range, the day is only semi-favourable.) It answers: how",
  "often have conditions recently favoured NEW infections? It is fast and",
  "event-based, and responds to humid, wet spells before disease builds up.",
  "",
  "Used together: BLASTAM flags when infection windows open (useful for timing a",
  "response), while EPIRICE estimates the disease that may result. Leaf wetness",
  "is estimated from hourly humidity (>=90%) and rain; both models read low in",
  "cool, dry weather.",
  "",
  "Bands and thresholds are provisional; calibrate against field observations.",
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

# ---- Pretty HTML email body ----------------------------------------------
band_bg <- function(l) switch(l, low = "#E8F1F8", moderate = "#FDB147",
                              high = "#E8492B", "#EFEFEF")
band_fg <- function(l) if (identical(l, "high")) "#FFFFFF" else "#22272B"
pill <- function(label, n, bg, fg) sprintf(
  paste0("<span style='display:inline-block;background:%s;color:%s;padding:4px 12px;",
         "border-radius:14px;margin-right:6px;font-size:13px;'>%s %d</span>"),
  bg, fg, label, n)

row_html <- character(0); cur <- ""
bl_bg <- function(n) if (is.na(n)) "#FFFFFF" else if (n == 0) "#F7F7FA" else
  if (n < 5) "#E5E1F0" else if (n < 10) "#CBC9E2" else if (n < 20) "#9E9AC8" else "#756BB1"
bl_fg <- function(n) if (!is.na(n) && n >= 20) "#FFFFFF" else "#22272B"
for (k in seq_len(nrow(results))) {
  r <- results[k]
  if (r$state != cur) {
    cur <- r$state
    row_html <- c(row_html, sprintf(
      paste0("<tr><td colspan='5' style='background:#F0F3F6;font-weight:bold;",
             "padding:6px 8px;border-top:1px solid #DDE2E6;'>%s</td></tr>"), cur))
  }
  pct <- if (is.na(r$intensity)) "-" else sprintf("%.2f%%", r$intensity * 100)
  tr  <- if (is.na(r$trend7)) "" else sprintf("%+0.2f", r$trend7 * 100)
  blv <- if (is.na(r$blast_events)) "-" else
    sprintf("%d<span style='color:#8a8f94;font-size:11px;'> (7d %d)</span>",
            r$blast_events, r$blast_recent)
  row_html <- c(row_html, sprintf(
    paste0("<tr><td style='padding:5px 8px;border-top:1px solid #EEE;'>%s</td>",
           "<td style='padding:5px 8px;border-top:1px solid #EEE;text-align:right;'>%s</td>",
           "<td style='padding:5px 8px;border-top:1px solid #EEE;'>",
           "<span style='background:%s;color:%s;padding:2px 9px;border-radius:10px;",
           "font-size:12px;'>%s</span></td>",
           "<td style='padding:5px 8px;border-top:1px solid #EEE;text-align:right;",
           "color:#6b7378;'>%s</td>",
           "<td style='padding:5px 8px;border-top:1px solid #EEE;text-align:right;",
           "background:%s;color:%s;'>%s</td></tr>"),
    r$name, pct, band_bg(r$level), band_fg(r$level), r$level, tr,
    bl_bg(r$blast_events), bl_fg(r$blast_events), blv))
}

html <- paste0(
"<div style='font-family:Arial,Helvetica,sans-serif;color:#22272B;max-width:720px;margin:auto;'>",
"<div style='background:#002664;color:#fff;padding:16px 20px;border-radius:6px 6px 0 0;'>",
"<div style='font-size:19px;font-weight:bold;'>Blast risk summary</div>",
sprintf("<div style='font-size:13px;opacity:.9;'>%s</div>", format(Sys.Date(), "%A %d %B %Y")),
"</div>",
"<div style='padding:16px 20px;border:1px solid #E0E0E0;border-top:none;border-radius:0 0 6px 6px;'>",
sprintf(paste0("<p style='margin:0 0 12px;'>Two models for %d monitoring towns, from ",
        "Open-Meteo ERA5 weather to %s. <b>EPIRICE intensity</b> is how much disease the ",
        "season has built up; <b>BLASTAM days</b> is how many of the last 21 days favoured a ",
        "new infection. Both are weather-driven potentials, not field measurements.</p>"),
        nrow(results), format(end_date, "%d %b %Y")),
if (!is.null(mg))
  sprintf(paste0("<p style='margin:-4px 0 12px;font-size:12px;color:#6b7378;'>",
                 "<b>Map:</b> %s.</p>"), mg) else "",
"<p style='margin:0 0 14px;'>",
"<span style='font-size:12px;color:#6b7378;margin-right:8px;'>EPIRICE bands:</span>",
pill("High", getn("high"), "#E8492B", "#fff"),
pill("Moderate", getn("moderate"), "#FDB147", "#22272B"),
pill("Low", getn("low"), "#E8F1F8", "#22272B"),
"</p>",
"<table style='border-collapse:collapse;width:100%;font-size:13px;'>",
"<thead><tr style='text-align:left;border-bottom:2px solid #002664;'>",
"<th style='padding:6px 8px;'>Town</th>",
"<th style='padding:6px 8px;text-align:right;'>EPIRICE %</th>",
"<th style='padding:6px 8px;'>Risk</th>",
"<th style='padding:6px 8px;text-align:right;'>7-day pts</th>",
"<th style='padding:6px 8px;text-align:right;'>BLASTAM days</th></tr></thead><tbody>",
paste(row_html, collapse = ""),
"</tbody></table>",
paste0("<p style='font-size:12px;color:#6b7378;margin:14px 0 0;'><b>EPIRICE</b> (intensity %) ",
       "mechanistically simulates the epidemic: the proportion of leaf tissue diseased. ",
       "<b>BLASTAM</b> (days) counts the days in the last 21 that were favourable for a NEW infection, when leaf wetness is ",
       "&ge;10&nbsp;h, the mean temperature during wetness is 15-32&deg;C, and the preceding ",
       "5-day mean temperature is 20-30&deg;C (upper bounds raised from Japan's 25&deg;C for warm conditions). EPIRICE says how much disease may build; BLASTAM ",
       "flags when infection windows open. Leaf wetness is estimated from humidity (&ge;90%) and ",
       "rain; both read low in cool, dry weather. Two maps and two trends CSVs are attached. ",
       "Values are provisional.</p>"),
paste0("<p style='font-size:11px;color:#9aa0a6;margin:10px 0 0;'>EPIRICE: Savary ",
       "<em>et al.</em> 2012 (Crop Prot. 34:6-17); epicrop (A.H. Sparks). BLASTAM: ",
       "Koshimizu 1988 (Bull. Tohoku Natl. Agric. Exp. Stn. 78:67-121); Hayashi &amp; ",
       "Koshimizu 1988. Weather: Open-Meteo ERA5 (CC BY 4.0).</p>"),
"</div></div>")
writeLines(html, file.path(OUT, "blast_summary_latest.html"))

# ---- Rolling trends CSVs (last HISTORY_RUNS runs, one column per run) ------
# Wide format: one row per town, one column per run date, so a row reads left to
# right as the trend. Persisted in the repo so they accumulate across runs.
write_trends <- function(values, file_name) {
  today <- data.table(town = results$name)
  today[[run_tag]] <- values
  f <- file.path(OUT, file_name)
  hist <- if (file.exists(f))
    fread(f, header = TRUE, colClasses = list(character = "town")) else
    data.table(town = character())
  if (run_tag %in% names(hist)) hist[, (run_tag) := NULL]   # re-run same day
  hist <- merge(hist, today, by = "town", all = TRUE)
  ord <- c(results$name, setdiff(hist$town, results$name))
  hist <- hist[match(ord, town)]
  keep_n <- if (exists("HISTORY_RUNS")) HISTORY_RUNS else 10L
  dcols <- setdiff(names(hist), "town")
  dcols <- dcols[order(as.Date(dcols))]
  if (length(dcols) > keep_n) dcols <- tail(dcols, keep_n)
  hist <- hist[, c("town", dcols), with = FALSE]
  fwrite(hist, f)
  cat(sprintf("Trends: %d towns x %d runs -> %s\n", nrow(hist), length(dcols), f))
}
write_trends(round(results$intensity * 100, 3), "town_trends.csv")     # EPIRICE %
write_trends(as.integer(results$blast_events),  "blastam_trends.csv")  # BLASTAM days

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

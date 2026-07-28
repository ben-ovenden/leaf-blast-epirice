#!/usr/bin/env Rscript
################################################################################
# run_blast.R
#
# Weekly rice leaf blast risk for the monitoring towns, using EPIRICE (Savary
# et al. 2012) and BLASTAM (Koshimizu 1988), driven by Open-Meteo weather.
#
#   Rscript run_blast.R
#
# Pipeline:
#   1. Fetch hourly weather for all towns in ONE batched request.
#   2. Run both models per town over the rolling CROP_AGE_DAYS window.
#   3. Report current intensity, its recent trend, a risk band, and BLASTAM days.
#   4. Write a CSV, a text summary and an HTML body for the email step.
#
# CHANGES
#   * Towns are fetched in a single batched request rather than 31 separate
#     forked requests, which removes 31 round trips and the fork pool.
#   * The mailR path is gone: the workflow sends mail with send_email.py.
#   * The map growth line now compares like with like (see map_stats.txt).
################################################################################

Sys.setenv(TZ = "Australia/Sydney")

SCRIPT_DIR <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile), winslash = "/"),
  error = function(e) normalizePath(getwd(), winslash = "/")
)

source(file.path(SCRIPT_DIR, "blast_config.R"))
source(file.path(SCRIPT_DIR, "epirice_model.R"))
source(file.path(SCRIPT_DIR, "blastam_model.R"))
source(file.path(SCRIPT_DIR, "openmeteo_wth.R"))
source(file.path(SCRIPT_DIR, "openmeteo_batch.R"))

suppressPackageStartupMessages(library(data.table))

classify <- function(intensity) {
  if (is.na(intensity)) return("no data")
  if (intensity < INTENSITY_LOW_MAX) return("low")
  if (intensity < INTENSITY_MODERATE_MAX) return("moderate")
  "high"
}

na_row_for <- function(name, lat, lon, level = "no data") {
  data.table(name = name, lat = lat, lon = lon,
             intensity = NA_real_, peak = NA_real_, trend7 = NA_real_,
             level = level, last_date = NA_character_, days = 0L,
             blast_events = NA_integer_, blast_recent = NA_integer_)
}

# Model one town from its daily rows.
model_town <- function(name, lat, lon, dd, emergence, end_date) {
  if (is.null(dd) || nrow(dd) < MIN_DAYS) return(na_row_for(name, lat, lon))

  bs <- blastam_score(dd$infect, dd$semi, as.Date(dd$date), end_date,
                      window = BLASTAM_WINDOW_DAYS, recent = 7L)

  wth <- data.table(YYYYMMDD = as.Date(dd$date),
                    DOY = as.integer(format(as.Date(dd$date), "%j")),
                    TEMP = dd$TEMP, RHUM = dd$RHUM, RAIN = dd$RAIN,
                    LAT = lat, LON = lon)
  setorder(wth, YYYYMMDD)
  duration <- as.integer(min(120L, nrow(wth)))
  lb <- tryCatch(predict_leaf_blast(wth, emergence = emergence, duration = duration),
                 error = function(e) NULL)
  if (is.null(lb) || nrow(lb) == 0) {
    r <- na_row_for(name, lat, lon)
    r[, `:=`(blast_events = bs$events, blast_recent = bs$recent)]
    return(r)
  }
  cur  <- lb$intensity[nrow(lb)]
  data.table(name = name, lat = lat, lon = lon,
             intensity = cur, peak = max(lb$intensity, na.rm = TRUE),
             trend7 = if (nrow(lb) > 7) cur - lb$intensity[nrow(lb) - 7] else NA_real_,
             level = classify(cur),
             last_date = as.character(lb$dates[nrow(lb)]),
             days = duration,
             blast_events = bs$events, blast_recent = bs$recent)
}

# ---- Run all towns ---------------------------------------------------------
cat("\nRice leaf blast risk (EPIRICE + BLASTAM / Open-Meteo)\n")
cat(strrep("=", 60), "\n")

end_date <- Sys.Date() - ARCHIVE_LAG_DAYS
# Same rolling window as the heatmap: assume a crop of age CROP_AGE_DAYS
# everywhere, so the table reflects current potential risk rather than a finished
# past season. Note this emergence date MOVES with every run, so the trends CSVs
# are a rolling window series, not a cumulative season total.
emergence <- end_date - CROP_AGE_DAYS
sites <- as.data.table(MONITOR_TOWNS)
sites[, pid := sprintf("%s", name)]

cat("Run date:   ", format(Sys.Date(), "%A %d %B %Y"), "\n")
cat("Weather to: ", format(end_date, "%Y-%m-%d"), " (archive lag ",
    ARCHIVE_LAG_DAYS, " days)\n", sep = "")
cat("Crop age:   ", CROP_AGE_DAYS, " days (rolling emergence ",
    format(emergence), ")\n\n", sep = "")

# One batched request covers all towns. Previously this was 31 separate forked
# requests, which cost 31 round trips for the same weighted quota.
town_daily <- new.env(parent = emptyenv())
on_town <- function(pid, lon, lat, hw) {
  w <- tryCatch(blastam_daily_from_hourly(hw), error = function(e) NULL)
  if (is.null(w) || nrow(w) == 0) return(NULL)
  assign(pid, as.data.table(w), envir = town_daily)
  data.table(pid = pid)          # non-empty return marks the point as ok
}
fr <- fetch_points_batched(sites[, .(pid, lon, lat)], emergence, end_date,
                           on_point = on_town, label = "towns")

# Serial fallback for towns the batch did not deliver, on a clean connection.
missing <- setdiff(sites$pid, ls(town_daily))
if (length(missing) > 0) {
  cat(sprintf("Serial fallback for %d town(s): %s\n",
              length(missing), paste(missing, collapse = ", ")))
  for (nm in missing) {
    s <- sites[pid == nm]
    hw <- tryCatch(get_openmeteo_hourly(s$lat, s$lon, emergence, end_date),
                   error = function(e) NULL)
    if (!is.null(hw) && nrow(hw) > 0) {
      w <- tryCatch(blastam_daily_from_hourly(hw), error = function(e) NULL)
      if (!is.null(w) && nrow(w) > 0) assign(nm, as.data.table(w), envir = town_daily)
    }
  }
}

results <- rbindlist(lapply(seq_len(nrow(sites)), function(k) {
  s <- sites[k]
  dd <- if (exists(s$pid, envir = town_daily, inherits = FALSE))
    get(s$pid, envir = town_daily) else NULL
  model_town(s$name, s$lat, s$lon, dd, emergence, end_date)
}))

for (k in seq_len(nrow(results))) {
  r <- results[k]
  msg <- if (is.na(r$intensity)) r$level else sprintf("%.1f%% (%s)", r$intensity * 100, r$level)
  cat(sprintf("%-14s ... %s\n", r$name, msg))
}

# Attach state and order by state then town
st <- as.data.table(MONITOR_TOWNS)[, .(name, state)]
results <- merge(results, st, by = "name", all.x = TRUE, sort = FALSE)
results[is.na(state), state := "--"]
setorder(results, state, name)

# ---- Summary ---------------------------------------------------------------
counts <- results[, .N, by = level]
getn <- function(lv) { x <- counts[level == lv, N]; if (length(x)) x else 0L }

# Read the grid's stats line (written by run_blast_grid.R).
# Fields: mapped | mean_land_spacing | mapped_last_run | finest | fmt | kb |
#         read_fmt | window_end | complete_res | weighted_spent
map_growth_line <- function() {
  f <- file.path(SCRIPT_DIR, OUTPUT_DIR, "map_stats.txt")
  if (!file.exists(f)) return(NULL)
  s <- tryCatch(strsplit(readLines(f, warn = FALSE)[1], "\\|")[[1]],
                error = function(e) NULL)
  if (length(s) < 4) return(NULL)
  now    <- as.integer(s[1])
  finest <- as.numeric(s[4])
  prev   <- as.integer(s[3])          # mapped LAST run, comparable with `now`
  fmt    <- if (length(s) >= 5) s[5] else NA
  kb     <- if (length(s) >= 6) as.numeric(s[6]) else NA
  rd     <- if (length(s) >= 7) s[7] else NA
  wend   <- if (length(s) >= 8) s[8] else NA
  cres   <- if (length(s) >= 9 && nzchar(s[9])) as.numeric(s[9]) else NA
  spent  <- if (length(s) >= 10) as.numeric(s[10]) else NA

  chg <- if (is.na(prev) || prev <= 0) ""
         else if (now > prev) sprintf(" (up from %d mapped last run)", prev)
         else if (now < prev) sprintf(" (down from %d mapped last run)", prev)
         else " (steady)"
  # State the COMPLETED lattice level, not a mean spacing. While the grid fills
  # coarse to fine, coverage is a complete coarse level plus a partial finer one,
  # and a single mean spacing implies a uniformity that is not there.
  cbit <- if (is.na(cres)) "" else
    sprintf(" Complete to %.2f deg (~%.0f km cells).", cres, cres * 111)
  wtxt <- if (is.na(wend)) "" else {
    d <- suppressWarnings(as.Date(wend))
    if (is.na(d)) "" else sprintf(" All cells modelled over one window to %s.", format(d, "%d %b"))
  }
  cache_bit <- if (is.na(fmt) || is.na(kb)) "" else {
    active <- if (!is.na(rd) && rd %in% c("gz", "csv")) rd
              else if (grepl("^gz", fmt)) "gz" else "csv"
    extra <- if (active == "gz" && grepl("csv", fmt)) ", csv backup kept"
             else if (active == "csv" && identical(fmt, "csv")) ", gz unavailable"
             else if (active == "csv" && grepl("gz", fmt)) ", gz restored this run"
             else ""
    sprintf(" Cache: %s active%s, %s.", active, extra,
            if (kb >= 1024) sprintf("%.1f MB", kb / 1024) else sprintf("%.0f KB", kb))
  }
  qbit <- if (is.na(spent)) "" else sprintf(" Used ~%.0f weighted API calls.", spent)
  sprintf("%d cells on a %.2f deg lattice%s.%s%s%s%s",
          now, finest, chg, cbit, wtxt, cache_bit, qbit)
}
mg <- map_growth_line()

midweek_line <- function() {
  f <- file.path(SCRIPT_DIR, OUTPUT_DIR, "midweek_status.txt")
  if (!file.exists(f)) return(NULL)
  s <- tryCatch(strsplit(readLines(f, warn = FALSE)[1], "\\|")[[1]],
                error = function(e) NULL)
  if (length(s) < 4) return(NULL)
  d <- suppressWarnings(as.Date(s[1])); pts <- as.integer(s[2]); added <- as.integer(s[4])
  if (is.na(d) || as.integer(Sys.Date() - d) > 2)
    sprintf("Daily top-up: no run in the last 2 days (last %s); check the top-up job.",
            if (is.na(d)) "never" else format(d, "%d %b"))
  else
    sprintf("Daily top-up ran %s: +%d points, cache now %d (ok).",
            format(d, "%d %b"), added, pts)
}
mw <- midweek_line()

summary_lines <- c(
  "Blast risk summary",
  strrep("=", 60),
  paste0("Generated:  ", format(Sys.Date(), "%A %d %B %Y")),
  paste0("Models:     EPIRICE (Savary et al. 2012) + BLASTAM (Koshimizu 1988)"),
  paste0("Weather:    Open-Meteo ", OPENMETEO_MODEL, " archive, to ", format(end_date, "%Y-%m-%d")),
  "",
  sprintf("EPIRICE bands  High:%d  Moderate:%d  Low:%d  No data:%d",
          getn("high"), getn("moderate"), getn("low"), getn("no data")),
  "",
  sprintf("Town detail: EPIRICE intensity + BLASTAM favourable days (last %d), by state:",
          BLASTAM_WINDOW_DAYS),
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
  "Reading the trends CSVs: emergence is taken as (latest weather date minus",
  sprintf("%d days) and therefore MOVES each run, so a row is a rolling %d day", CROP_AGE_DAYS, CROP_AGE_DAYS),
  "window through time, not a cumulative season total.",
  "",
  "Bands and thresholds are provisional; calibrate against field observations.",
  "",
  if (exists("CITATION")) CITATION else NULL,
  if (!is.null(mg) || !is.null(mw)) "" else NULL,
  if (!is.null(mg)) paste0("Map:     ", mg) else NULL,
  if (!is.null(mw)) paste0("Top-up:  ", mw) else NULL,
  "",
  "Sent by Ben Ovenden, ben.ovenden@dpird.nsw.gov.au")

summary_text <- paste(summary_lines, collapse = "\n")
cat("\n", summary_text, "\n", sep = "")

# ---- Write outputs ---------------------------------------------------------
OUT <- file.path(SCRIPT_DIR, OUTPUT_DIR)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
run_tag <- format(Sys.Date(), "%Y-%m-%d")

fwrite(results, file.path(OUT, sprintf("blast_results_%s.csv", run_tag)))
writeLines(summary_text, file.path(OUT, sprintf("blast_summary_%s.txt", run_tag)))
fwrite(results, file.path(OUT, "blast_results_latest.csv"))
writeLines(summary_text, file.path(OUT, "blast_summary_latest.txt"))

# ---- HTML email body -------------------------------------------------------
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
  tr  <- if (is.na(r$trend7)) "" else sprintf("%+0.2f", { v <- r$trend7 * 100; if (abs(v) < 0.005) 0 else v })
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
        "Open-Meteo %s weather to %s. <b>EPIRICE intensity</b> is how much disease the ",
        "rolling %d day window has built up; <b>BLASTAM days</b> is how many of the last %d days ",
        "favoured a new infection. Both are weather-driven potentials, not field measurements.</p>"),
        nrow(results), OPENMETEO_MODEL, format(end_date, "%d %b %Y"),
        CROP_AGE_DAYS, BLASTAM_WINDOW_DAYS),
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
"<th style='padding:6px 8px;text-align:right;'>7 day pts</th>",
"<th style='padding:6px 8px;text-align:right;'>BLASTAM days</th></tr></thead><tbody>",
paste(row_html, collapse = ""),
"</tbody></table>",
paste0("<p style='font-size:12px;color:#6b7378;margin:14px 0 0;'><b>EPIRICE</b> (intensity %) ",
       "mechanistically simulates the epidemic: the proportion of leaf tissue diseased. ",
       "<b>BLASTAM</b> (days) counts the days in the last 21 that were favourable for a NEW infection, when leaf wetness is ",
       "&ge;10&nbsp;h, the mean temperature during wetness is 15-32&deg;C, and the preceding ",
       "5-day mean temperature is 20-30&deg;C (upper bounds raised from Japan's 25&deg;C for warm conditions). EPIRICE says how much disease may build; BLASTAM ",
       "flags when infection windows open. Leaf wetness is estimated from humidity (&ge;90%) and ",
       "rain; both read low in cool, dry weather. Emergence moves with each run, so the trends ",
       "CSVs are a rolling window rather than a season total. Two maps and two trends CSVs are ",
       "attached. Values are provisional.</p>"),
paste0("<p style='font-size:11px;color:#9aa0a6;margin:10px 0 0;'>EPIRICE: Savary ",
       "<em>et al.</em> 2012 (Crop Prot. 34:6-17); epicrop (A.H. Sparks). BLASTAM: ",
       "Koshimizu 1988 (Bull. Tohoku Natl. Agric. Exp. Stn. 78:67-121); Hayashi &amp; ",
       "Koshimizu 1988. Weather: Open-Meteo ERA5 (CC BY 4.0), non-commercial research use.</p>"),
if (!is.null(mg))
  sprintf(paste0("<p style='font-size:11px;color:#9aa0a6;margin:8px 0 0;'>",
                 "<b>Map:</b> %s</p>"), mg) else "",
if (!is.null(mw))
  sprintf(paste0("<p style='font-size:11px;color:#9aa0a6;margin:2px 0 0;'>",
                 "<b>Top-up:</b> %s</p>"), mw) else "",
"<p style='font-size:11px;color:#b0b5ba;margin:8px 0 0;'>Sent by Ben Ovenden, ",
"<a href='mailto:ben.ovenden@dpird.nsw.gov.au' style='color:#b0b5ba;'>",
"ben.ovenden@dpird.nsw.gov.au</a></p>",
"</div></div>")
writeLines(html, file.path(OUT, "blast_summary_latest.html"))

# ---- Rolling trends CSVs ---------------------------------------------------
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
  dcols <- setdiff(names(hist), "town")
  # Guard against a stray non-date column corrupting the ordering.
  parsed <- suppressWarnings(as.Date(dcols))
  dcols <- dcols[!is.na(parsed)][order(parsed[!is.na(parsed)])]
  keep_n <- HISTORY_RUNS
  if (length(dcols) > keep_n) dcols <- tail(dcols, keep_n)
  hist <- hist[, c("town", dcols), with = FALSE]
  fwrite(hist, f)
  cat(sprintf("Trends: %d towns x %d runs -> %s\n", nrow(hist), length(dcols), f))
}
write_trends(round(results$intensity * 100, 3), "town_trends.csv")     # EPIRICE %
write_trends(as.integer(results$blast_events),  "blastam_trends.csv")  # BLASTAM days

# ---- Optional simple town point map ----------------------------------------
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
       main = sprintf("Rice leaf blast risk, weather to %s", format(end_date)))
  grid(col = "#E4E7E9")
  points(results$lon, results$lat, pch = 21, cex = 3,
         bg = band_col(results$level), col = "#22272B", lwd = 1.5)
  text(results$lon, results$lat, labels = results$name, pos = 3,
       offset = 1, cex = 0.8, col = "#22272B")
  legend("topright", legend = c("Low", "Moderate", "High", "No data"),
         pt.bg = c(COL_LOW, COL_MODERATE, COL_HIGH, COL_NODATA),
         pch = 21, pt.cex = 2, bty = "n")
  par(op); dev.off()
  file.copy(map_file, file.path(OUT, "blast_map_latest.png"), overwrite = TRUE)
  cat("\nMap:     ", map_file, "\n")
}

cat("Results: ", file.path(OUT, sprintf("blast_results_%s.csv", run_tag)), "\n")
cat("\nDone.\n")

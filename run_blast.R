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
#   3. Report current intensity, its recent trend, a risk band, BLASTAM days and
#      the number of nights that could not be judged.
#   4. Write a CSV, a text summary and an HTML body for the email step.
#
# CHANGES IN THIS REVISION
#   * The run date comes from blast_run_date() once, so the town table and the
#     maps can no longer be dated differently when a run straddles midnight.
#   * The BLASTAM window is bounded at both ends inside blastam_score(), so the
#     town figure is the same 21 day count the map reports, not 22.
#   * The daily series is checked for calendar gaps before SEIR sees it. SEIR
#     indexes by position; only the grid runner used to screen for this.
#   * The town fetch is budgeted against the shared weighted-spend ledger.
#   * Reporting fixes: the wetness threshold is described as the temperature
#     dependent curve actually used, not a fixed 10 h; intensity is printed to
#     three decimals so the 0.2% and 1% band edges are resolvable; the trend
#     column shows a dash rather than "+0.00" for no change; unjudged nights are
#     surfaced; bands use the NSW Government palette; trends columns are keyed on
#     the DATA end date and blanks are written as NA.
################################################################################

SCRIPT_DIR <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile), winslash = "/"),
  error = function(e) normalizePath(getwd(), winslash = "/")
)

source(file.path(SCRIPT_DIR, "blast_config.R"))
source(file.path(SCRIPT_DIR, "epirice_model.R"))
source(file.path(SCRIPT_DIR, "blastam_model.R"))
source(file.path(SCRIPT_DIR, "openmeteo_wth.R"))
source(file.path(SCRIPT_DIR, "openmeteo_batch.R"))

suppressPackageStartupMessages({library(data.table); library(methods)})

# ---- One run date, one window ----------------------------------------------
RUN_DATE  <- blast_run_date()
data_end  <- RUN_DATE - ARCHIVE_LAG_DAYS      # last day FETCHED
end_date  <- data_end - DAY_CUT_LAG_DAYS      # last day MODELLED
emergence <- end_date - CROP_AGE_DAYS
run_tag   <- format(RUN_DATE, "%Y-%m-%d")
data_tag  <- format(end_date, "%Y-%m-%d")

classify <- function(intensity) {
  if (is.na(intensity)) return("no data")
  if (intensity < INTENSITY_LOW_MAX) return("low")
  if (intensity < INTENSITY_MODERATE_MAX) return("moderate")
  "high"
}

# Intensity is a fraction of leaf sites. Three decimals as a percentage, because
# the band edges are 0.2% and 1% and one decimal cannot resolve either: the
# 2026-07-30 email printed "0.00%" beside a "low" band for every town.
fmt_pct <- function(x) if (is.na(x)) "-" else sprintf("%.3f%%", x * 100)
fmt_chg <- function(x) {
  if (is.na(x)) return("-")
  v <- x * 100
  if (abs(v) < 0.0005) "-" else sprintf("%+.3f", v)
}

na_row_for <- function(name, lat, lon, level = "no data") {
  data.table(name = name, lat = lat, lon = lon,
             intensity = NA_real_, peak = NA_real_, trend7 = NA_real_,
             level = level, last_date = NA_character_, days = 0L,
             blast_events = NA_integer_, blast_recent = NA_integer_,
             blast_semi = NA_integer_, blast_unjudged = NA_integer_,
             blast_window = NA_integer_, note = "no usable weather")
}

# SEIR indexes the weather vector by POSITION, so a missing calendar day shifts
# every later day against crop age. Return the longest continuous run ending at
# the last available date, and report what was dropped.
trim_to_continuous <- function(dd) {
  setorder(dd, date)
  n <- nrow(dd)
  if (n < 2L) return(dd)
  gaps <- which(as.integer(diff(dd$date)) != 1L)
  if (length(gaps) == 0L) return(dd)
  dd[(max(gaps) + 1L):n]
}

# Model one town from its daily rows.
model_town <- function(name, lat, lon, dd, emergence, end_date) {
  if (is.null(dd) || nrow(dd) == 0L) return(na_row_for(name, lat, lon))
  dd <- as.data.table(dd)[date <= end_date]
  if (nrow(dd) == 0L) return(na_row_for(name, lat, lon))

  # BLASTAM is scored on the full series, gaps and all: an unjudgeable night is
  # already NA, so a gap lowers `unjudged` rather than corrupting the count.
  lagd <- if (exists("BLASTAM_END_LAG_DAYS")) BLASTAM_END_LAG_DAYS else 0L
  bs <- blastam_score(dd$infect, dd$semi, dd$date, end_date - lagd,
                      window = BLASTAM_WINDOW_DAYS, recent = BLASTAM_RECENT_DAYS)
  semi_n <- sum(dd$semi[dd$date > (end_date - lagd - BLASTAM_WINDOW_DAYS) &
                        dd$date <= (end_date - lagd)], na.rm = TRUE)

  bl <- function(r) { r[, `:=`(blast_events = bs$events, blast_recent = bs$recent,
                               blast_semi = as.integer(semi_n),
                               blast_unjudged = bs$unjudged,
                               blast_window = as.integer(bs$n_days))]; r }

  cont <- trim_to_continuous(dd)
  dropped <- nrow(dd) - nrow(cont)
  if (nrow(cont) < MIN_DAYS)
    return(bl(na_row_for(name, lat, lon))[, note := "series too short for EPIRICE"])

  wth <- data.table(YYYYMMDD = cont$date,
                    DOY = as.integer(format(cont$date, "%j")),
                    TEMP = cont$TEMP, RHUM = cont$RHUM, RAIN = cont$RAIN,
                    LAT = lat, LON = lon)
  setorder(wth, YYYYMMDD)
  duration <- as.integer(min(120L, nrow(wth)))
  lb <- tryCatch(predict_leaf_blast(wth, emergence = wth$YYYYMMDD[1], duration = duration),
                 error = function(e) NULL)
  if (is.null(lb) || nrow(lb) == 0)
    return(bl(na_row_for(name, lat, lon))[, note := "EPIRICE did not converge"])

  cur <- lb$intensity[nrow(lb)]
  bl(data.table(name = name, lat = lat, lon = lon,
                intensity = cur, peak = max(lb$intensity, na.rm = TRUE),
                trend7 = if (nrow(lb) > 7) cur - lb$intensity[nrow(lb) - 7] else NA_real_,
                level = classify(cur),
                last_date = as.character(lb$dates[nrow(lb)]),
                days = duration,
                blast_events = NA_integer_, blast_recent = NA_integer_,
                blast_semi = NA_integer_, blast_unjudged = NA_integer_,
                blast_window = NA_integer_,
                note = if (dropped > 0)
                  sprintf("%d day(s) before a gap excluded", dropped) else ""))
}

# ---- Run all towns ---------------------------------------------------------
cat("\nRice leaf blast risk (EPIRICE + BLASTAM / Open-Meteo)\n")
cat(strrep("=", 60), "\n")
blastam_check_fetch_arithmetic()

sites <- as.data.table(MONITOR_TOWNS)
sites[, pid := sprintf("%s", name)]

cat("Run date:   ", format(RUN_DATE, "%A %d %B %Y"), "\n")
cat("Fetched to: ", format(data_end, "%Y-%m-%d"), " (archive lag ", ARCHIVE_LAG_DAYS,
    " days)\n", sep = "")
cat("Modelled to:", format(end_date, "%Y-%m-%d"), " (model day cut at ",
    BLASTAM_DAY_CUT_HOUR, ":00 local solar, so the last fetched day is partial)\n", sep = "")
cat("Crop age:   ", CROP_AGE_DAYS, " days (rolling emergence ",
    format(emergence), ")\n\n", sep = "")

# One batched request covers all towns, with BLASTAM_LEADIN_DAYS of lead-in that
# is computed and then discarded, and lon passed so the model day and night
# window are in local solar time.
fetch_from <- emergence - BLASTAM_LEADIN_DAYS
OUT <- file.path(SCRIPT_DIR, OUTPUT_DIR)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ledger_path <- file.path(OUT, SPEND_LEDGER_FILE)
already <- om_spend_read(ledger_path)
town_budget <- max(0, DAILY_WEIGHTED_HARD_CAP - already)
cat(sprintf("Weighted quota: %.0f already spent today, %.0f available to this run.\n",
            already, town_budget))

town_daily <- new.env(parent = emptyenv())
on_town <- function(pid, lon, lat, hw) {
  w <- tryCatch(blastam_daily_from_hourly(hw, lon = lon), error = function(e) NULL)
  if (is.null(w) || nrow(w) == 0) return(NULL)
  w <- w[date >= emergence & date <= end_date]
  if (nrow(w) == 0) return(NULL)
  assign(pid, as.data.table(w), envir = town_daily)
  data.table(pid = pid)          # non-empty return marks the point as ok
}
fr <- fetch_points_batched(sites[, .(pid, lon, lat)], fetch_from, data_end,
                           on_point = on_town, budget = town_budget, label = "towns")
om_spend_add(ledger_path, fr$spent, "towns")

# Serial fallback for towns the batch did not deliver, on a clean connection.
missing <- setdiff(sites$pid, ls(town_daily))
if (length(missing) > 0 && fr$spent < town_budget) {
  cat(sprintf("Serial fallback for %d town(s): %s\n",
              length(missing), paste(missing, collapse = ", ")))
  for (nm in missing) {
    s <- sites[pid == nm]
    hw <- tryCatch(get_openmeteo_hourly(s$lat, s$lon, fetch_from, data_end),
                   error = function(e) NULL)
    if (!is.null(hw) && nrow(hw) > 0) {
      w <- tryCatch(blastam_daily_from_hourly(hw, lon = s$lon), error = function(e) NULL)
      if (!is.null(w) && nrow(w) > 0) {
        w <- w[date >= emergence & date <= end_date]
        if (nrow(w) > 0) assign(nm, as.data.table(w), envir = town_daily)
      }
    }
  }
}

results <- rbindlist(lapply(seq_len(nrow(sites)), function(k) {
  s <- sites[k]
  dd <- if (exists(s$pid, envir = town_daily, inherits = FALSE))
    get(s$pid, envir = town_daily) else NULL
  model_town(s$name, s$lat, s$lon, dd, emergence, end_date)
}), fill = TRUE)

for (k in seq_len(nrow(results))) {
  r <- results[k]
  msg <- if (is.na(r$intensity)) r$level else sprintf("%s (%s)", fmt_pct(r$intensity), r$level)
  cat(sprintf("%-16s ... %s\n", r$name, msg))
}

# Attach state and order by state then town
st <- as.data.table(MONITOR_TOWNS)[, .(name, state)]
results <- merge(results, st, by = "name", all.x = TRUE, sort = FALSE)
results[is.na(state), state := "--"]
setorder(results, state, name)

# ---- Summary ---------------------------------------------------------------
counts <- results[, .N, by = level]
getn <- function(lv) { x <- counts[level == lv, N]; if (length(x)) x else 0L }
n_unjudged <- sum(results$blast_unjudged, na.rm = TRUE)

wet_rule <- if (isTRUE(BLASTAM_USE_BJ_THRESHOLD)) {
  paste0("reaches a temperature dependent minimum (about ",
         sprintf("%.0f", blastam_bj_min_hours(16)), " h at 16 C falling to about ",
         sprintf("%.0f", blastam_bj_min_hours(27)), " h at 27 C; Barksdale and Jones 1965)")
} else {
  sprintf("reaches at least %d h (Koshimizu's fixed threshold)", BLASTAM_WET_HOURS_FIXED)
}

# Read the grid's stats line (written by run_blast_grid.R).
# Fields: mapped | mean_land_spacing | mapped_last_run | finest | fmt | kb |
#         read_fmt | window_end | complete_res | weighted_spent | obs_max_epi |
#         obs_max_blastam
map_growth_line <- function() {
  f <- file.path(OUT, "map_stats.txt")
  if (!file.exists(f)) return(NULL)
  s <- tryCatch(strsplit(readLines(f, warn = FALSE)[1], "\\|")[[1]],
                error = function(e) NULL)
  if (length(s) < 4) return(NULL)
  fld <- function(i, f = identity) if (length(s) >= i && nzchar(s[i])) f(s[i]) else NA
  now    <- as.integer(s[1]);  finest <- as.numeric(s[4]); prev <- as.integer(s[3])
  fmt    <- fld(5); kb <- fld(6, as.numeric); rd <- fld(7); wend <- fld(8)
  cres   <- fld(9, as.numeric); spent <- fld(10, as.numeric)
  mx_epi <- fld(11, as.numeric); mx_bl <- fld(12, as.numeric)

  chg <- if (is.na(prev) || prev <= 0) ""
         else if (now > prev) sprintf(" (up from %d mapped last run)", prev)
         else if (now < prev) sprintf(" (down from %d mapped last run)", prev)
         else " (steady)"
  cbit <- if (is.na(cres)) "" else
    sprintf(" Complete to %.2f deg (~%.0f km cells).", cres, cres * 111)
  wtxt <- if (is.na(wend)) "" else {
    d <- suppressWarnings(as.Date(wend))
    if (is.na(d)) "" else sprintf(" All cells modelled over one window to %s.", format(d, "%d %b"))
  }
  # Observed maxima, so a flat map is distinguishable from a broken one.
  mbit <- if (is.na(mx_epi) && is.na(mx_bl)) "" else
    sprintf(" Observed maxima: EPIRICE %s, BLASTAM %s.",
            if (is.na(mx_epi)) "NA" else sprintf("%.3f%%", mx_epi),
            if (is.na(mx_bl))  "NA" else sprintf("%.0f days", mx_bl))
  cache_bit <- if (is.na(fmt) || is.na(kb)) "" else {
    active <- if (!is.na(rd) && rd %in% c("gz", "csv")) rd
              else if (grepl("^gz", fmt)) "gz" else "csv"
    sprintf(" Cache: %s active, %s.", active,
            if (kb >= 1024) sprintf("%.1f MB", kb / 1024) else sprintf("%.0f KB", kb))
  }
  qbit <- if (is.na(spent)) "" else sprintf(" Used ~%.0f weighted API calls.", spent)
  sprintf("%d cells on a %.2f deg lattice%s.%s%s%s%s%s",
          now, finest, chg, cbit, wtxt, mbit, cache_bit, qbit)
}
mg <- map_growth_line()

midweek_line <- function() {
  f <- file.path(OUT, "midweek_status.txt")
  if (!file.exists(f)) return(NULL)
  s <- tryCatch(strsplit(readLines(f, warn = FALSE)[1], "\\|")[[1]],
                error = function(e) NULL)
  if (length(s) < 4) return(NULL)
  d <- suppressWarnings(as.Date(s[1])); pts <- as.integer(s[2]); added <- as.integer(s[4])
  if (is.na(d) || as.integer(RUN_DATE - d) > 2)
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
  paste0("Generated:   ", format(RUN_DATE, "%A %d %B %Y")),
  paste0("Models:      EPIRICE (Savary et al. 2012) + BLASTAM (Koshimizu 1988)"),
  paste0("Weather:     Open-Meteo ERA5 archive, fetched to ", format(data_end, "%Y-%m-%d"),
         ", modelled to ", format(end_date, "%Y-%m-%d")),
  sprintf("EPIRICE:     %d day window, RcT infection optimum %d C",
          CROP_AGE_DAYS, EPIRICE_RCT_PEAK),
  sprintf("Cache schema: version %d", CACHE_SCHEMA_VERSION),
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
  bl <- if (is.na(r$blast_events)) "-" else
    sprintf("%d days (%dd %d)", r$blast_events, BLASTAM_RECENT_DAYS, r$blast_recent)
  uj <- if (!is.na(r$blast_unjudged) && r$blast_unjudged > 0)
    sprintf(" [%d night(s) not judged]", r$blast_unjudged) else ""
  summary_lines <- c(summary_lines,
    sprintf("  %-16s EPIRICE %8s %-9s  BLASTAM %s%s",
            r$name, fmt_pct(r$intensity), r$level, bl, uj))
}
summary_lines <- c(summary_lines, "",
  "About these estimates: two different models",
  strrep("-", 60),
  "EPIRICE (intensity %) mechanistically simulates the whole epidemic. It steps",
  sprintf("through a crop %d days old at the final day, tracking healthy, latent,", CROP_AGE_DAYS),
  "infectious and removed leaf sites, and reports the proportion of leaf tissue",
  "diseased. Intensity is (diseased - removed) / (total sites - removed) on the",
  "final day. It is slow and cumulative and reads near zero in cool, dry weather.",
  "",
  "The EPIRICE wetness gate opens on a daily MEAN RH of 90% or a daily rain SUM",
  "of 5 mm. A 24 hour mean of 90% requires an all day saturated air mass and is",
  "rarely reached, so in practice this configuration is driven by the rain branch.",
  "",
  "BLASTAM (infection days) is the Japanese infection-warning model (Koshimizu",
  "1988). For each night it judges whether conditions favoured a NEW infection and",
  sprintf("counts the favourable nights in the last %d days. A night is favourable when", BLASTAM_WINDOW_DAYS),
  paste0("all three hold: leaf wetness ", wet_rule, ";"),
  sprintf("mean temperature during wetness %d-%d C; and the PRECEDING 5 day mean",
          BLASTAM_TWET_MIN, BLASTAM_TWET_MAX),
  sprintf("temperature %d-%d C (upper bounds raised from Japan's 25 C for warm",
          BLASTAM_PREV5_MIN, BLASTAM_PREV5_MAX),
  "conditions). If the wetness bar is met but one temperature condition fails, the",
  "night is only semi-favourable and is reported separately in the CSV.",
  "",
  sprintf("The model day runs from %02d:00 local solar so a night's wet period and",
          BLASTAM_DAY_CUT_HOUR),
  "rainfall sit inside one day rather than being split at midnight.",
  "",
  "Nights that could not be judged, from missing hours or insufficient lead-in,",
  "are counted separately and are NOT scored as unfavourable.",
  if (n_unjudged > 0)
    sprintf("This run: %d unjudged night(s) across all towns.", n_unjudged)
  else "This run: no unjudged nights.",
  "",
  "Used together: BLASTAM flags when infection windows open (useful for timing a",
  "response), while EPIRICE estimates the disease that may result.",
  "",
  strwrap(CAVEAT_CANOPY, width = 78),
  "",
  "Reading the trends CSVs: emergence is taken as (last modelled date minus",
  sprintf("%d days) and therefore MOVES each run, so a row is a rolling %d day", CROP_AGE_DAYS, CROP_AGE_DAYS),
  "window through time, not a cumulative season total. Columns are keyed on the",
  "DATA end date, not the run date, and run_log.csv records the schema and model",
  "options behind each column so a change of method is not read as a change in",
  "the weather.",
  "",
  "Bands and thresholds are provisional; calibrate against field observations.",
  "",
  if (exists("CITATION")) CITATION else NULL,
  if (!is.null(mg) || !is.null(mw)) "" else NULL,
  if (!is.null(mg)) paste0("Map:     ", mg) else NULL,
  if (!is.null(mw)) paste0("Top-up:  ", mw) else NULL,
  "",
  paste0("Sent by ", CONTACT_NAME, ", ", CONTACT_EMAIL))

summary_text <- paste(summary_lines, collapse = "\n")
cat("\n", summary_text, "\n", sep = "")

# ---- Write outputs ---------------------------------------------------------
fwrite(results, file.path(OUT, sprintf("blast_results_%s.csv", run_tag)), na = "NA")
writeLines(summary_text, file.path(OUT, sprintf("blast_summary_%s.txt", run_tag)))
fwrite(results, file.path(OUT, "blast_results_latest.csv"), na = "NA")
writeLines(summary_text, file.path(OUT, "blast_summary_latest.txt"))

# ---- HTML email body -------------------------------------------------------
band <- function(l) { s <- BAND_STYLE[[l]]; if (is.null(s)) BAND_STYLE[["no data"]] else s }
chip <- function(l) {
  s <- band(l)
  sprintf(paste0("<span style='display:inline-block;background:%s;color:%s;padding:2px 9px;",
                 "border-radius:10px;font-size:12px;white-space:nowrap;'>",
                 "<span style='display:inline-block;width:7px;height:7px;border-radius:50%%;",
                 "background:%s;margin-right:5px;'></span>%s</span>"),
          s$bg, s$fg, s$dot, l)
}
pill <- function(l, n) {
  s <- band(l)
  sprintf(paste0("<span style='display:inline-block;background:%s;color:%s;padding:4px 12px;",
                 "border-radius:14px;margin-right:6px;font-size:13px;'>",
                 "<span style='display:inline-block;width:8px;height:8px;border-radius:50%%;",
                 "background:%s;margin-right:6px;'></span>%s %d</span>"),
          s$bg, s$fg, s$dot, tools::toTitleCase(l), n)
}
bl_bg <- function(n) {
  if (is.na(n)) return(BLASTAM_CELL_BG[1])
  i <- findInterval(n, c(0, 1, 3, 6, 11, 16)) 
  BLASTAM_CELL_BG[max(1L, min(length(BLASTAM_CELL_BG), i))]
}
bl_fg <- function(n) if (!is.na(n) && n >= 11) "#FFFFFF" else NSW_BRAND_DARK

row_html <- character(0); cur <- ""
for (k in seq_len(nrow(results))) {
  r <- results[k]
  if (r$state != cur) {
    cur <- r$state
    row_html <- c(row_html, sprintf(
      paste0("<tr><td colspan='5' style='background:%s;font-weight:bold;",
             "padding:6px 8px;border-top:1px solid %s;'>%s</td></tr>"),
      NSW_GREY_01, NSW_GREY_02, cur))
  }
  blv <- if (is.na(r$blast_events)) "-" else
    sprintf("%d<span style='color:#6b7378;font-size:11px;'> (%dd %d)</span>%s",
            r$blast_events, BLASTAM_RECENT_DAYS, r$blast_recent,
            if (!is.na(r$blast_unjudged) && r$blast_unjudged > 0)
              sprintf("<span title='%d night(s) could not be judged' style='color:%s;'>&nbsp;*</span>",
                      r$blast_unjudged, COL_MODERATE) else "")
  row_html <- c(row_html, sprintf(
    paste0("<tr><td style='padding:5px 8px;border-top:1px solid #EEE;'>%s</td>",
           "<td style='padding:5px 8px;border-top:1px solid #EEE;text-align:right;",
           "font-variant-numeric:tabular-nums;'>%s</td>",
           "<td style='padding:5px 8px;border-top:1px solid #EEE;'>%s</td>",
           "<td style='padding:5px 8px;border-top:1px solid #EEE;text-align:right;",
           "color:#6b7378;'>%s</td>",
           "<td style='padding:5px 8px;border-top:1px solid #EEE;text-align:right;",
           "background:%s;color:%s;'>%s</td></tr>"),
    r$name, fmt_pct(r$intensity), chip(r$level), fmt_chg(r$trend7),
    bl_bg(r$blast_events), bl_fg(r$blast_events), blv))
}

bl_legend <- paste0(
  "<span style='font-size:12px;color:#6b7378;margin-right:8px;'>BLASTAM days:</span>",
  paste(mapply(function(lab, bg, fg) sprintf(
    paste0("<span style='display:inline-block;background:%s;color:%s;padding:3px 9px;",
           "border-radius:3px;margin-right:4px;font-size:12px;'>%s</span>"), bg, fg, lab),
    c("0", "1-2", "3-5", "6-10", "11-15", "16-21"),
    BLASTAM_CELL_BG,
    c(rep(NSW_BRAND_DARK, 4), "#FFFFFF", "#FFFFFF")), collapse = ""))

html <- paste0(
"<div style='font-family:Arial,Helvetica,sans-serif;color:", NSW_BRAND_DARK,
";max-width:760px;margin:auto;'>",
"<div style='background:", NSW_BRAND_BLUE, ";color:#fff;padding:16px 20px;border-radius:6px 6px 0 0;'>",
"<div style='font-size:19px;font-weight:bold;'>Blast risk summary</div>",
sprintf("<div style='font-size:13px;opacity:.9;'>%s</div>", format(RUN_DATE, "%A %d %B %Y")),
"</div>",
"<div style='padding:16px 20px;border:1px solid #E0E0E0;border-top:none;border-radius:0 0 6px 6px;'>",
sprintf(paste0("<p style='margin:0 0 12px;'>Two models for %d monitoring towns, from ",
        "Open-Meteo ERA5 weather fetched to %s and modelled to <b>%s</b>. ",
        "<b>EPIRICE intensity</b> is how much disease the rolling %d day window has built up; ",
        "<b>BLASTAM days</b> is how many of the last %d days favoured a new infection. ",
        "Both are weather-driven potentials, not field measurements.</p>"),
        nrow(results), format(data_end, "%d %b %Y"), format(end_date, "%d %b %Y"),
        CROP_AGE_DAYS, BLASTAM_WINDOW_DAYS),
"<p style='margin:0 0 8px;'>",
"<span style='font-size:12px;color:#6b7378;margin-right:8px;'>EPIRICE bands:</span>",
pill("high", getn("high")), pill("moderate", getn("moderate")), pill("low", getn("low")),
if (getn("no data") > 0) pill("no data", getn("no data")) else "",
"</p>",
"<p style='margin:0 0 14px;'>", bl_legend, "</p>",
"<table style='border-collapse:collapse;width:100%;font-size:13px;'>",
sprintf("<thead><tr style='text-align:left;border-bottom:2px solid %s;'>", NSW_BRAND_BLUE),
"<th style='padding:6px 8px;'>Town</th>",
"<th style='padding:6px 8px;text-align:right;'>EPIRICE intensity</th>",
"<th style='padding:6px 8px;'>Risk</th>",
"<th style='padding:6px 8px;text-align:right;'>7 day change (pts)</th>",
sprintf("<th style='padding:6px 8px;text-align:right;'>BLASTAM days (last %d)</th></tr></thead><tbody>",
        BLASTAM_WINDOW_DAYS),
paste(row_html, collapse = ""),
"</tbody></table>",
if (n_unjudged > 0) sprintf(
  paste0("<p style='font-size:12px;color:%s;margin:10px 0 0;'>* %d night(s) could not be ",
         "judged this run, from missing hours or insufficient lead-in. These are NOT counted ",
         "as unfavourable.</p>"), COL_MODERATE, n_unjudged) else "",
sprintf(paste0("<p style='font-size:12px;color:#6b7378;margin:14px 0 0;'><b>EPIRICE</b> ",
       "mechanistically simulates the epidemic and reports the proportion of leaf tissue ",
       "diseased. Its wetness gate needs a daily mean RH of 90%% or a daily rain total of ",
       "5 mm; the humidity branch rarely opens, so in practice the signal is rain driven. ",
       "<b>BLASTAM</b> counts the nights in the last %d that favoured a NEW infection: leaf ",
       "wetness %s, mean temperature during wetness %d-%d&deg;C, and the preceding 5 day mean ",
       "temperature %d-%d&deg;C (upper bounds raised from Japan's 25&deg;C for warm ",
       "conditions). Wetness is estimated from hourly humidity (&ge;%d%%) and rain, over a ",
       "model day starting at %02d:00 local solar. Both read low in cool, dry weather. ",
       "Emergence moves with each run, so the trends CSVs are a rolling window rather than a ",
       "season total; their columns are keyed on the data end date and run_log.csv records the ",
       "method behind each one. Values are provisional.</p>"),
       BLASTAM_WINDOW_DAYS, wet_rule, BLASTAM_TWET_MIN, BLASTAM_TWET_MAX,
       BLASTAM_PREV5_MIN, BLASTAM_PREV5_MAX, BLASTAM_RH_WET, BLASTAM_DAY_CUT_HOUR),
sprintf("<p style='font-size:12px;color:#6b7378;margin:10px 0 0;'>%s</p>", CAVEAT_CANOPY),
paste0("<p style='font-size:11px;color:#9aa0a6;margin:10px 0 0;'>EPIRICE: Savary ",
       "<em>et al.</em> 2012 (Crop Prot. 34:6-17); epicrop (A.H. Sparks); RcT infection ",
       sprintf("optimum %d&deg;C. BLASTAM: Koshimizu 1988 (Bull. Tohoku Natl. Agric. Exp. ",
               EPIRICE_RCT_PEAK),
       "Stn. 78:67-121); Hayashi &amp; Koshimizu 1988; wetness threshold from Barksdale &amp; ",
       "Jones 1965. Prior Australian modelling: Lanoiselet <em>et al.</em> 2002 (Australas. ",
       "Plant Pathol. 31:1-7). Weather: Open-Meteo ERA5 (CC BY 4.0), non-commercial research ",
       "use.</p>"),
if (!is.null(mg))
  sprintf(paste0("<p style='font-size:11px;color:#9aa0a6;margin:8px 0 0;'>",
                 "<b>Map:</b> %s</p>"), mg) else "",
if (!is.null(mw))
  sprintf(paste0("<p style='font-size:11px;color:#9aa0a6;margin:2px 0 0;'>",
                 "<b>Top-up:</b> %s</p>"), mw) else "",
sprintf(paste0("<p style='font-size:11px;color:#b0b5ba;margin:8px 0 0;'>Sent by %s, ",
               "<a href='mailto:%s' style='color:#b0b5ba;'>%s</a></p>"),
        CONTACT_NAME, CONTACT_EMAIL, CONTACT_EMAIL),
"</div></div>")
writeLines(html, file.path(OUT, "blast_summary_latest.html"))

# ---- Rolling trends CSVs ---------------------------------------------------
# Wide format: one row per town, one column per DATA end date.
#
# Columns used to be keyed on the RUN date, so three test runs on 28, 29 and 30
# July each took a column while describing almost the same weather, and with a
# short history that evicts the genuinely older columns. Keying on the data end
# date means a re-run over the same window replaces its column instead.
#
# Blanks are written as NA. In the delivered 2026-07-11 column an absent town
# read as an empty cell, which is indistinguishable from a zero in a spreadsheet.
write_trends <- function(values, file_name) {
  today <- data.table(town = results$name)
  today[[data_tag]] <- values
  f <- file.path(OUT, file_name)
  hist <- if (file.exists(f))
    fread(f, header = TRUE, colClasses = list(character = "town"), na.strings = "NA") else
    data.table(town = character())
  if (data_tag %in% names(hist)) hist[, (data_tag) := NULL]
  hist <- merge(hist, today, by = "town", all = TRUE)
  ord <- c(results$name, setdiff(hist$town, results$name))
  hist <- hist[match(ord, town)]
  dcols <- setdiff(names(hist), "town")
  parsed <- suppressWarnings(as.Date(dcols))
  dcols <- dcols[!is.na(parsed)][order(parsed[!is.na(parsed)])]
  if (length(dcols) > HISTORY_RUNS) dcols <- tail(dcols, HISTORY_RUNS)
  hist <- hist[, c("town", dcols), with = FALSE]
  fwrite(hist, f, na = "NA")
  cat(sprintf("Trends: %d towns x %d data windows -> %s\n", nrow(hist), length(dcols), f))
}
write_trends(round(results$intensity * 100, 4), "town_trends.csv")     # EPIRICE %
write_trends(as.integer(results$blast_events),  "blastam_trends.csv")  # BLASTAM days

# ---- Run log ---------------------------------------------------------------
# One row per run recording the method behind each trends column, so a change in
# schema or parameters is visible in the series rather than looking like a real
# epidemiological collapse. This is what was missing when Malanda dropped from
# 0.374% to 0.006% between the 2026-07-28 and 2026-07-29 columns.
log_row <- data.table(
  run_date = format(RUN_DATE), data_end = data_tag, fetched_to = format(data_end),
  crop_age_days = CROP_AGE_DAYS, cache_schema = CACHE_SCHEMA_VERSION,
  rct_peak_c = EPIRICE_RCT_PEAK, day_cut_hour = BLASTAM_DAY_CUT_HOUR,
  blastam_window = BLASTAM_WINDOW_DAYS,
  bj_threshold = BLASTAM_USE_BJ_THRESHOLD,
  twet = sprintf("%d-%d", BLASTAM_TWET_MIN, BLASTAM_TWET_MAX),
  prev5 = sprintf("%d-%d", BLASTAM_PREV5_MIN, BLASTAM_PREV5_MAX),
  towns_modelled = sum(!is.na(results$intensity)),
  unjudged_nights = n_unjudged,
  weighted_spent_towns = round(fr$spent, 1))
log_f <- file.path(OUT, RUN_LOG_FILE)
# Append safely. fread() infers types, so a column written as character (a date
# formatted with format()) comes back as Date and rbind() then refuses with
# "Class attribute on column 3 ... does not match". Align the old rows to the new
# row's classes before binding, and never let a metadata write kill a run whose
# real outputs are already on disk.
append_run_log <- function(f, row) {
  old <- if (file.exists(f)) tryCatch(fread(f), error = function(e) NULL) else NULL
  if (!is.null(old) && nrow(old) > 0) {
    if ("data_end" %in% names(old)) old <- old[as.character(data_end) != data_tag]
    for (cl in intersect(names(row), names(old)))
      if (!identical(class(old[[cl]]), class(row[[cl]])))
        set(old, j = cl, value = methods::as(as.character(old[[cl]]), class(row[[cl]])[1]))
  }
  fwrite(if (is.null(old) || nrow(old) == 0) row else rbind(old, row, fill = TRUE),
         f, na = "NA")
}
tryCatch(append_run_log(log_f, log_row),
         error = function(e) cat("Run log not updated:", conditionMessage(e), "\n"))
cat("Run log: ", log_f, "\n")

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
         bg = band_col(results$level), col = NSW_BRAND_DARK, lwd = 1.5)
  text(results$lon, results$lat, labels = results$name, pos = 3,
       offset = 1, cex = 0.8, col = NSW_BRAND_DARK)
  legend("topright", legend = c("Low", "Moderate", "High", "No data"),
         pt.bg = c(COL_LOW, COL_MODERATE, COL_HIGH, COL_NODATA),
         pch = 21, pt.cex = 2, bty = "n")
  par(op); dev.off()
  file.copy(map_file, file.path(OUT, "blast_map_latest.png"), overwrite = TRUE)
  cat("\nMap:     ", map_file, "\n")
}

cat("Results: ", file.path(OUT, sprintf("blast_results_%s.csv", run_tag)), "\n")
cat("\nDone.\n")

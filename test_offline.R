#!/usr/bin/env Rscript
################################################################################
# test_offline.R
#
# Offline regression tests. No network, no API quota. Run from the repo root:
#
#   Rscript test_offline.R
#
# Each test guards a bug that was actually shipped at some point, so a failure
# here means a real regression rather than a style complaint. Several of the
# tests below guard a specific class of failure that reached a delivered email.
################################################################################
suppressPackageStartupMessages({library(data.table); library(methods)})
SCRIPT_DIR <- normalizePath(".", winslash = "/")
source("blast_config.R"); source("epirice_model.R"); source("blastam_model.R")

fails <- 0L; n <- 0L
ok <- function(label, cond, extra = "") {
  n <<- n + 1L
  if (isTRUE(cond)) cat(sprintf("  PASS  %s\n", label))
  else { fails <<- fails + 1L; cat(sprintf("  FAIL  %s %s\n", label, extra)) }
}
have <- function(f) exists(f, mode = "function")

# Synthetic hourly UTC series whose humidity follows a LOCAL diurnal cycle.
mkseries <- function(lon, days = 20, wet_from = 19, wet_to = 8, tbase = 24,
                     rain_mm = 0, rain_from = 22, rain_to = 3, rain_every = 0L) {
  dt <- as.POSIXct("2026-06-01 00:00", tz = "UTC") + (0:(days * 24 - 1)) * 3600
  loc <- (as.numeric(format(dt, "%H")) + lon / 15) %% 24
  d <- as.integer(as.Date(format(dt, "%Y-%m-%d")))
  wet <- loc >= wet_from | loc < wet_to
  rr <- if (rain_every > 0L)
    ifelse(d %% rain_every == 0L & (loc >= rain_from | loc < rain_to), rain_mm, 0) else 0
  data.table(dt = dt,
             temp = tbase + 4 * cos((loc - 15) / 24 * 2 * pi),
             rh   = ifelse(wet, 95, 55),
             rain = rr)
}

cat("\n1. BLASTAM night window is LOCAL SOLAR, not UTC\n")
# Regression: applying the 15:00-09:00 window to UTC timestamps put it at roughly
# 01:00-19:00 local over Australia, returning 8 wet hours instead of 13.
for (lon in c(115, 130, 145, 153)) {
  d <- blastam_daily_from_hourly(mkseries(lon), lon = lon)
  j <- d[!is.na(infect)]
  ok(sprintf("lon %d: wetness 12-13 h detected", lon),
     nrow(j) > 0 && all(j$wet_hours >= 12) && all(j$wet_hours <= 13),
     sprintf("(got %s)", paste(unique(j$wet_hours), collapse = "/")))
  ok(sprintf("lon %d: every judged night favourable", lon),
     nrow(j) > 0 && all(j$infect == 1L),
     sprintf("(got %d of %d)", sum(j$infect), nrow(j)))
}

cat("\n2. The preceding 5-day mean genuinely precedes\n")
# Regression: prev5 was frollmean(TEMP, 5, align = "right"), which covers days
# i-4 to i and so includes the night's own day. Five leading NAs, not four, is
# the signature of a correctly lagged window.
d <- blastam_daily_from_hourly(mkseries(145), lon = 145)
ok("first 5 days NA (no PRECEDING 5-day mean yet)",
   sum(cumprod(is.na(d$infect))) == 5L,
   sprintf("(got %d leading NAs)", sum(cumprod(is.na(d$infect)))))
ok("lead-in constant covers the shift plus the lagged mean",
   BLASTAM_LEADIN_DAYS >= 6L)
ok("short wetness gives 0, not NA", {
  s <- blastam_daily_from_hourly(mkseries(145, wet_from = 23, wet_to = 3), lon = 145)
  sj <- s[!is.na(infect)]; nrow(sj) > 0 && all(sj$infect == 0L)
})

cat("\n3. Unjudgeable nights are NA, and partial model days are dropped\n")
ok("no day carries a truncated aggregate",
   all(diff(as.integer(d$date)) == 1L),
   sprintf("(got %d rows spanning %d days)", nrow(d),
           as.integer(diff(range(d$date))) + 1L))
# Regression: complete was n_eve > 0 & n_morn > 0, so a night observed for two
# hours was judged rather than left NA.
ok("a sparsely observed night is NA, not 0", {
  h <- mkseries(145, days = 12)
  # keep only two hours of one night, drop the rest of that night
  target <- as.Date("2026-06-06")
  loc <- (as.numeric(format(h$dt, "%H")) + 145 / 15) %% 24
  md <- as.Date(h$dt + (145 / 15 - BLASTAM_DAY_CUT_HOUR) * 3600)
  drop <- md == target & (loc >= 15 | loc < 9)
  keepers <- which(drop)[1:2]
  drop[keepers] <- FALSE
  h2 <- h[!drop]
  r <- blastam_daily_from_hourly(h2, lon = 145)
  is.null(r) || !isTRUE(target %in% r$date) || is.na(r[date == target, infect])
})
# Regression: an hour with missing humidity counted as dry.
ok("a location with no humidity at all is rejected, not scored dry", {
  source("openmeteo_batch.R")
  el <- list(hourly = list(time = as.list(format(
    as.POSIXct("2026-06-01 00:00", tz = "UTC") + (0:47) * 3600, "%Y-%m-%dT%H:%M")),
    temperature_2m = as.list(rep(25, 48)),
    relative_humidity_2m = vector("list", 48),
    precipitation = as.list(rep(0, 48))))
  is.null(.om_hourly_dt(el))
})

cat("\n4. The model day keeps a night's rain in one day\n")
# Regression: schema 2 cut the model day at local midnight, which split a
# nocturnal rain event across two days and halved the peak daily total. EPIRICE's
# rainlim gate is a daily SUM, so days reaching 5 mm went from 24 to 0 and
# Malanda fell from 0.374% to 0.006% between the 2026-07-28 and 2026-07-29 runs.
h <- mkseries(145.6, days = 40, rain_mm = 1.4, rain_every = 3L)
dd <- blastam_daily_from_hourly(h, lon = 145.6)
ok("nocturnal rain lands in one model day, above rainlim",
   max(dd$RAIN) >= 5, sprintf("(max daily rain %.1f mm)", max(dd$RAIN)))
ok("the day cut is not midnight", BLASTAM_DAY_CUT_HOUR != 0L)
ok("splitting at midnight would have hidden it", {
  x <- copy(h); x[, sdt := dt + 145.6 / 15 * 3600]; x[, cd := as.Date(sdt)]
  keep <- x[, .N, by = cd][N >= 24L, cd]
  max(x[cd %in% keep, .(r = sum(rain)), by = cd]$r) < 5
})

cat("\n5. BLASTAM scoring window is bounded at BOTH ends\n")
# Regression: inwin was `dates > (end_date - window)` with no upper bound, so
# run_blast.R reported a 22 day count and an 8 day "7d" count while the map used
# the correct form. The two products in one email disagreed.
dts <- seq(as.Date("2026-06-01"), as.Date("2026-07-24"), by = "day")
inf <- rep(1L, length(dts)); sem <- rep(0L, length(dts))
bs <- blastam_score(inf, sem, dts, as.Date("2026-07-23"), window = 21L, recent = 7L)
ok("21 day window counts 21 days", bs$events == 21L, sprintf("(got %d)", bs$events))
ok("7 day window counts 7 days",  bs$recent == 7L,  sprintf("(got %d)", bs$recent))
ok("rows after end_date are excluded",
   blastam_score(inf, sem, dts, as.Date("2026-07-10"), window = 21L)$events == 21L)
ok("the window width is reported", bs$n_days == 21L)

cat("\n6. EPIRICE RcT curve matches the configured optimum\n")
# Regression: the README argued at length for the published 25 C peak while this
# file shipped epicrop's 20 C curve, a factor of about two at 28 C.
ok(sprintf("configured peak is %d C", EPIRICE_RCT_PEAK),
   EPIRICE_RCT_PEAK %in% c(20L, 25L))
ok("the configured curve peaks where it says it does",
   epirice_rct()[which.max(epirice_rct()[, 2]), 1] == EPIRICE_RCT_PEAK,
   sprintf("(peaks at %d)", epirice_rct()[which.max(epirice_rct()[, 2]), 1]))
ok("both curves are available and differ at 28 C",
   abs(.fn_Rc(epirice_rct(25), 28) - 0.76) < 1e-9 &&
   abs(.fn_Rc(epirice_rct(20), 28) - 0.36) < 1e-9,
   sprintf("(25C %.2f, 20C %.2f)", .fn_Rc(epirice_rct(25), 28), .fn_Rc(epirice_rct(20), 28)))

cat("\n7. EPIRICE date alignment and gap safety\n")
end_date <- as.Date("2026-07-22"); emergence <- end_date - CROP_AGE_DAYS
mkw <- function(from, to) {
  w <- data.table(YYYYMMDD = seq(from, to, by = "day"))
  w[, `:=`(DOY = as.integer(format(YYYYMMDD, "%j")), TEMP = 24, RHUM = 92,
           RAIN = 1, LAT = -25, LON = 148)][]
}
w <- mkw(emergence, end_date)
ok("inclusive window models without error",
   !is.null(tryCatch(predict_leaf_blast(w, emergence, nrow(w)), error = function(e) NULL)))
ok("exclusive window is the failure mode this guards",
   is.null(tryCatch(predict_leaf_blast(w[-1], emergence, nrow(w) - 1L),
                    error = function(e) NULL)))
# Regression: run_blast_grid.R passed the run's global emergence while truncating
# the weather to an earlier model_end, so SEIR threw for every point and the
# EPIRICE map rendered empty while BLASTAM rendered normally.
model_end <- end_date - 5L
wc <- mkw(model_end - CROP_AGE_DAYS, model_end)
ok("coverage-mode window fails with the run's global emergence",
   is.null(tryCatch(predict_leaf_blast(wc, emergence, nrow(wc)), error = function(e) NULL)))
ok("coverage-mode window works with emergence derived from model_end",
   !is.null(tryCatch(predict_leaf_blast(wc, model_end - CROP_AGE_DAYS, nrow(wc)),
                     error = function(e) NULL)))
# Regression: SEIR indexes the weather BY POSITION.
ok("SEIR refuses a series with a calendar gap",
   is.null(tryCatch(predict_leaf_blast(w[-30], emergence, nrow(w) - 1L),
                    error = function(e) NULL)))

cat("\n8. Weighted cost model and the 14 day arithmetic\n")
if (!have("om_weight_per_location")) {
  ok("om_weight_per_location() is defined", FALSE, "(openmeteo_batch.R stubbed?)")
} else {
  ok("14 day floor: 1 day costs the same as 14",
     om_weight_per_location(1, 3) == om_weight_per_location(14, 3))
  ok("11+ variables trigger the multiplier", om_weight_per_location(14, 11) > 1)
  # A refresh must cost EXACTLY 1.00, or a 7% surcharge on every point eats the
  # headroom for adding new cells once the grid is full.
  ref_days <- REFRESH_TAIL_DAYS + BLASTAM_LEADIN_DAYS + DAY_CUT_LAG_DAYS
  ok("tail + lead-in + day-cut lag comes to exactly 14 days", ref_days == 14L,
     sprintf("(got %d)", ref_days))
  ok("so a refresh costs exactly 1.00 weighted",
     abs(om_weight_per_location(ref_days, 3) - 1) < 1e-12)
  add_days <- CROP_AGE_DAYS + 1L + BLASTAM_LEADIN_DAYS + DAY_CUT_LAG_DAYS
  ok(sprintf("a new point costs %.2f weighted over %d days",
             om_weight_per_location(add_days, 3), add_days),
     abs(om_weight_per_location(add_days, 3) - add_days / 14) < 1e-9)
}

cat("\n9. Pacer and the shared spend ledger\n")
if (!have("om_pacer")) {
  ok("om_pacer() is defined", FALSE, "(openmeteo_batch.R stubbed?)")
} else {
  p <- om_pacer(6000); p(6000)
  t0 <- Sys.time(); for (i in 1:5) p(100)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok("500 weighted at 6000/min takes ~5 s", abs(el - 5) < 2.0, sprintf("(got %.2fs)", el))
}
# Regression: DAILY_WEIGHTED_CAP was applied per RUN, and run_blast.R fetched its
# towns with budget = Inf on top of whatever the grid run had already spent.
if (have("om_spend_add")) {
  lf <- tempfile(fileext = ".csv")
  om_spend_add(lf, 8667, "grid"); om_spend_add(lf, 151, "towns")
  ok("the spend ledger accumulates within a UTC day",
     abs(om_spend_read(lf) - 8818) < 1e-6, sprintf("(got %.0f)", om_spend_read(lf)))
  ok("and the combined total stays under the free daily ceiling",
     om_spend_read(lf) < FREE_DAILY_CALLS)
  unlink(lf)
} else {
  ok("om_spend_add() is defined", FALSE)
}

cat("\n10. The run log appends across runs\n")
# Regression: the run log wrote dates with format() (character) but fread() read
# them back as Date, so rbind() refused on the SECOND run with "Class attribute
# on column 3 does not match" and the log never gained a row.
{
  lf <- tempfile(fileext = ".csv")
  mkrow <- function(rd, de) data.table(run_date = rd, data_end = de,
    fetched_to = de, cache_schema = CACHE_SCHEMA_VERSION, rct_peak_c = EPIRICE_RCT_PEAK,
    bj_threshold = BLASTAM_USE_BJ_THRESHOLD, twet = "15-32")
  appender <- function(f, row) {
    old <- if (file.exists(f)) tryCatch(fread(f), error = function(e) NULL) else NULL
    if (!is.null(old) && nrow(old) > 0) {
      old <- old[as.character(data_end) != row$data_end[1]]
      for (cl in intersect(names(row), names(old)))
        if (!identical(class(old[[cl]]), class(row[[cl]])))
          set(old, j = cl, value = methods::as(as.character(old[[cl]]), class(row[[cl]])[1]))
    }
    fwrite(if (is.null(old) || nrow(old) == 0) row else rbind(old, row, fill = TRUE),
           f, na = "NA")
  }
  e1 <- tryCatch({ appender(lf, mkrow("2026-07-30", "2026-07-23")); NULL },
                 error = function(e) conditionMessage(e))
  e2 <- tryCatch({ appender(lf, mkrow("2026-08-06", "2026-07-30")); NULL },
                 error = function(e) conditionMessage(e))
  got <- if (file.exists(lf)) nrow(fread(lf)) else 0L
  ok("a second run appends rather than erroring", is.null(e1) && is.null(e2) && got == 2L,
     sprintf("(rows %d; %s)", got, paste(c(e1, e2), collapse = "; ")))
  ok("re-running the same data window replaces its row", {
    appender(lf, mkrow("2026-08-07", "2026-07-30")); nrow(fread(lf)) == 2L })
  unlink(lf)
}

cat("\n11. Cache gz is really gzipped\n")
# Regression: the atomic write used a ".tmp" temp name, so fwrite() stopped
# compressing and an 8x larger plain file was committed under a .gz name.
f <- tempfile(fileext = ".tmp.gz"); fwrite(data.table(a = 1:500, b = strrep("x", 20)), f)
magic <- as.integer(readBin(f, "raw", 2L))
ok("fwrite compresses when the extension survives", identical(magic, c(31L, 139L)),
   sprintf("(magic %s)", paste(magic, collapse = " ")))
f2 <- tempfile(fileext = ".tmp"); fwrite(data.table(a = 1:500, b = strrep("x", 20)), f2)
ok("and does NOT when it is stripped (the bug)",
   !identical(as.integer(readBin(f2, "raw", 2L)), c(31L, 139L)))
unlink(c(f, f2))

cat("\n12. Documented controls are actually wired up\n")
# Regression: HEAT_STRETCH was documented in the README as the fix for the flat
# map and no script read it, so both delivered maps rendered as one pale blue.
# Comments are stripped, so a note that MENTIONS a banned call is not a hit.
# Stripping from the first "#" can only remove text, never create a match.
code_of <- function(f) paste(sub("#.*$", "", readLines(f, warn = FALSE)), collapse = "\n")
gsrc <- code_of("run_blast_grid.R")
tsrc <- code_of("run_blast.R")
ok("HEAT_STRETCH is read by the renderer", grepl("HEAT_STRETCH", gsrc, fixed = TRUE))
ok("BLASTAM_STRETCH is read by the renderer", grepl("BLASTAM_STRETCH", gsrc, fixed = TRUE))
ok("COAST_MASK_KM is read by the renderer", grepl("COAST_MASK_KM", gsrc, fixed = TRUE))
ok("OVERLAY_MAX_SEGMENT_DEG is read by the renderer",
   grepl("OVERLAY_MAX_SEGMENT_DEG", gsrc, fixed = TRUE))
# Regression: the email footnote stated a fixed 10 h wetness threshold while the
# temperature dependent curve was in use.
ok("no hard-coded wetness threshold in the email prose",
   !grepl("wetness is (&ge;|>=)\\s*10", tsrc) && !grepl("wetness >=10", tsrc))
ok("both runners take the run date from blast_run_date()",
   grepl("blast_run_date()", gsrc, fixed = TRUE) &&
   grepl("blast_run_date()", tsrc, fixed = TRUE))
# Regression: run_tag came from a second Sys.Date() called after a multi-hour
# fetch, so a run straddling midnight dated the map and the table differently.
ok("neither runner calls Sys.Date() directly",
   !grepl("Sys.Date()", gsrc, fixed = TRUE) && !grepl("Sys.Date()", tsrc, fixed = TRUE))

cat("\n13. A starved run still models from cache\n")
# Regression: under GRID_WINDOW_MODE = "latest", model_end was pinned to
# end_date unconditionally. The 2026-07-31 run was correctly capped at 4 weighted
# calls by the shared ledger after an earlier run the same UTC day, so no cell
# reached end_date, nothing was modelled, no map rendered and the email step
# failed on a missing attachment, with 1,874 usable cached points in the repo.
source("grid_window.R")    # the real pick_window_end(), not a copy
{
  ed <- as.Date("2026-07-24")
  # every cached cell is one day behind, exactly the failing case
  stale <- data.table(pid = sprintf("p%03d", 1:200), mx = ed - 1L)
  reach <- sum(stale$mx >= ed) / nrow(stale)
  ok("the failing case is detected", reach < GRID_WINDOW_MIN_COVERAGE,
     sprintf("(%.0f%% reach end_date)", 100 * reach))
  fb <- pick_window_end(stale, ed, GRID_WINDOW_MIN_COVERAGE)
  ok("the window falls back one day", fb == ed - 1L, sprintf("(got %s)", format(fb)))
  ok("and every cached cell then reaches it", sum(stale$mx >= fb) == nrow(stale))
  # a healthy run must NOT fall back
  fresh <- data.table(pid = sprintf("p%03d", 1:200), mx = ed)
  ok("a healthy run keeps end_date",
     sum(fresh$mx >= ed) / nrow(fresh) >= GRID_WINDOW_MIN_COVERAGE)
  # a mixed cache: 95% current, 5% stale. Should keep end_date, not fall back.
  mixed <- data.table(pid = sprintf("p%03d", 1:200),
                      mx = c(rep(ed, 190), rep(ed - 3L, 10)))
  ok("a mostly current cache keeps end_date",
     sum(mixed$mx >= ed) / nrow(mixed) >= GRID_WINDOW_MIN_COVERAGE)
  # empty cache must not error
  ok("an empty cache does not error",
     pick_window_end(data.table(pid=character(), mx=as.Date(character())), ed, 0.9) == ed)
}
# The heatmaps must be OPTIONAL attachments, so a degraded run still emails.
esrc <- paste(readLines("send_email.py", warn = FALSE), collapse = "\n")
opt_block <- sub(".*optional = \\[", "", esrc)
opt_block <- sub("\\].*", "", opt_block)
ok("heatmaps are optional attachments, not required",
   grepl("heatmap", opt_block, fixed = TRUE))
ok("trends CSVs are still required",
   grepl("town_trends.csv", sub("\\].*", "", sub(".*required = \\[", "", esrc)), fixed = TRUE))

cat("\n14. Grid lattice covers the requested extent\n")
fin <- GRID_RES_FINEST
lat_top <- GRID_EXTENT[3] + ceiling((GRID_EXTENT[4] - GRID_EXTENT[3]) / fin) * fin
ok("the northernmost row is inside the lattice", lat_top >= GRID_EXTENT[4],
   sprintf("(lattice reaches %.2f, extent asks for %.2f)", lat_top, GRID_EXTENT[4]))
ok("the old seq() would have dropped it",
   max(seq(GRID_EXTENT[3], GRID_EXTENT[4], by = fin)) < GRID_EXTENT[4] ||
   isTRUE(all.equal((GRID_EXTENT[4] - GRID_EXTENT[3]) %% fin, 0)))

cat("\n15. Overlay artefacts and label declutter (terra)\n")
if (!requireNamespace("terra", quietly = TRUE)) {
  cat("  SKIP  terra not installed\n")
} else {
  suppressPackageStartupMessages(library(terra))
  # THE MASKING REGRESSION. run_blast_grid.R loads terra after data.table, and
  # terra::shift masks data.table::shift. The bare call failed, the caller's
  # tryCatch turned it into every point coming back "empty", and the run produced
  # a blank map with no error in the log.
  dt2 <- blastam_daily_from_hourly(mkseries(145), lon = 145)
  ok("the aggregator still works with terra attached",
     !is.null(dt2) && nrow(dt2) > 0 && any(!is.na(dt2$infect)))

  # australia_roads.geojson holds a feature with a 3.25 deg step from Victoria to
  # Tasmania, which drew as a line across Bass Strait on every map.
  if (file.exists("australia_roads.geojson")) {
    v <- terra::vect("australia_roads.geojson")
    g <- as.data.frame(terra::geom(v))
    longest <- max(vapply(split(g, paste(g$geom, g$part)), function(s)
      if (nrow(s) < 2) 0 else max(sqrt(diff(s$x)^2 + diff(s$y)^2)), numeric(1)))
    ok("the bundled roads layer really does contain a long jump",
       longest > OVERLAY_MAX_SEGMENT_DEG, sprintf("(longest %.2f deg)", longest))
  } else {
    ok("australia_roads.geojson present", FALSE)
  }

  declutter_labels <- function(lon, lat, minsep) {
    keep <- logical(length(lon)); px <- numeric(0); py <- numeric(0)
    for (i in order(lat)) {
      if (length(px) == 0L || all(sqrt((lon[i] - px)^2 + (lat[i] - py)^2) >= minsep)) {
        keep[i] <- TRUE; px <- c(px, lon[i]); py <- c(py, lat[i])
      }
    }
    keep
  }
  tw <- as.data.frame(MONITOR_TOWNS)
  kp <- declutter_labels(tw$lon, tw$lat, LABEL_MIN_SEP_DEG)
  ok("declutter drops overprinting town labels", sum(!kp) > 0 && sum(kp) > 20,
     sprintf("(kept %d of %d)", sum(kp), nrow(tw)))
  ok("declutter is deterministic",
     identical(kp, declutter_labels(tw$lon, tw$lat, LABEL_MIN_SEP_DEG)))
}

cat("\n16. The data window follows the UTC date, not only the Sydney date\n")
# Regression: the Monday job fires at ~22:30 UTC on Sunday, so "Sydney run date
# minus 6" was only five days behind the archive's clock, and that day was not in
# the archive yet. With the model day cut at 10:00 local solar, every cell west
# of 120 E needed more than the two permitted hours from the missing day to
# complete end_date and landed one day short: 841 cells, "only 89% of cached
# cells reached end_date" on every scheduled run, a Degraded run banner on every
# Monday email, and a week later a full-price refetch for each of them.
{
  old <- Sys.getenv(c("BLAST_RUN_DATE", "BLAST_UTC_DATE"), unset = NA)
  Sys.setenv(BLAST_RUN_DATE = "2026-09-21", BLAST_UTC_DATE = "2026-09-20")
  ok("the Monday run (fired while it is still Sunday UTC) counts the lag from the UTC date",
     blast_data_end(blast_run_date()) == as.Date("2026-09-20") - ARCHIVE_LAG_DAYS,
     sprintf("(got %s)", format(blast_data_end(blast_run_date()))))
  Sys.setenv(BLAST_UTC_DATE = "2026-09-21")
  ok("a run after 10:00 Sydney (same UTC date) is unchanged",
     blast_data_end(blast_run_date()) == as.Date("2026-09-21") - ARCHIVE_LAG_DAYS)
  Sys.setenv(BLAST_UTC_DATE = "2026-09-25")
  ok("a manual re-run of an earlier run date keeps that date's window",
     blast_data_end(blast_run_date()) == as.Date("2026-09-21") - ARCHIVE_LAG_DAYS)
  for (v in names(old))
    if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, as.list(setNames(old[[v]], v)))
  ok("both runners take the window from blast_data_end()",
     grepl("blast_data_end(", gsrc, fixed = TRUE) && grepl("blast_data_end(", tsrc, fixed = TRUE))
  ok("neither runner subtracts ARCHIVE_LAG_DAYS from the run date itself",
     !grepl("RUN_DATE\\s*-\\s*ARCHIVE_LAG_DAYS", gsrc) &&
     !grepl("RUN_DATE\\s*-\\s*ARCHIVE_LAG_DAYS", tsrc))
  # The geometry behind it. With the final UTC day absent, a cell at 150 E still
  # completes its last model day (that day ends at 23:59 UTC), while a cell at
  # 115 E needs three hours of the missing day and may lose only two.
  h <- mkseries(115, days = 20)
  h[dt >= as.POSIXct("2026-06-20 00:00", tz = "UTC"),
    `:=`(temp = NA_real_, rh = NA_real_, rain = NA_real_)]
  d150 <- blastam_daily_from_hourly(copy(h), lon = 150)
  d115 <- blastam_daily_from_hourly(copy(h), lon = 115)
  ok("150 E completes model day 19 June without the 20th",
     max(d150$date) == as.Date("2026-06-19"), sprintf("(got %s)", format(max(d150$date))))
  ok("115 E does not: it needs 3 hours of the missing day and may lose only 2",
     max(d115$date) == as.Date("2026-06-18"), sprintf("(got %s)", format(max(d115$date))))
}

cat("\n17. Stale cells are refetched over the days they miss, not the whole crop window\n")
# Regression: any cell more than REFRESH_TAIL_DAYS behind was refetched over the
# full window at ~4.86 weighted. After the 2026-09-07 run was cut short at 5,000
# weighted, 3,243 cells were in that class; the weekly budget recovered ~800 of
# them per run and the map sat on 29 August for three weeks (65%, 74%, 76%).
source("openmeteo_batch.R")
{
  ed <- as.Date("2026-09-14"); de <- ed + DAY_CUT_LAG_DAYS
  akf <- ed - CROP_AGE_DAYS
  lb <- data.table(pid = sprintf("p%04d", 1:1000),
                   last = c(rep(ed, 50),          # already current
                            rep(ed - 7L, 800),    # the weekly cadence
                            rep(ed - 8L, 100),    # one day short (the western cells)
                            rep(ed - 22L, 40),    # cut off on 09-07
                            rep(ed - 90L, 10)))   # older than the crop window
  pl <- om_plan_refresh(lb, ed, de, akf, lead = BLASTAM_LEADIN_DAYS, n_vars = 3)
  ok("cells already at end_date are not eligible", nrow(pl) == 950L, sprintf("(got %d)", nrow(pl)))
  cost_of <- function(m) pl[missing == m, unique(cost)]
  ok("7 days behind costs exactly 1.00 (the 14 day floor)", abs(cost_of(7) - 1) < 1e-12)
  ok("8 days behind costs 15/14, not 4.86", abs(cost_of(8) - 15 / 14) < 1e-9,
     sprintf("(got %.3f)", cost_of(8)))
  ok("22 days behind costs 29/14", abs(cost_of(22) - 29 / 14) < 1e-9,
     sprintf("(got %.3f)", cost_of(22)))
  add_cost <- om_weight_per_location(CROP_AGE_DAYS + 1L + BLASTAM_LEADIN_DAYS + DAY_CUT_LAG_DAYS, 3)
  ok("older than the crop window costs the same as a new cell, never more",
     abs(cost_of(90) - add_cost) < 1e-9 && all(pl$cost <= add_cost + 1e-9))
  ok("each cohort starts the day after its newest row (no gap, no overlap)",
     pl[missing < 90, all(keep_from == last + 1L)])
  ok("and fetches BLASTAM_LEADIN_DAYS of lead-in before that",
     all(pl$fetch_from == pl$keep_from - BLASTAM_LEADIN_DAYS))
  ok("cheapest cohort first", !is.unsorted(pl$cost))
  # Budget for the 800 weekly cells, the 100 one-day-short cells and half of the
  # 22-day cohort: the prefix must stop mid-cohort, not skip a cohort.
  bud <- 800 + 100 * 15 / 14 + 20 * 29 / 14 + 0.01
  pb <- om_plan_refresh(lb, ed, de, akf, lead = BLASTAM_LEADIN_DAYS, n_vars = 3, budget = bud)
  ok("the budget buys an affordable prefix", pb[take == TRUE, .N] == 920L,
     sprintf("(took %d)", pb[take == TRUE, .N]))
  ok("spending stays within it", pb[take == TRUE, sum(cost)] <= bud)
  ok("what is left is the dearest cohort", pb[take == FALSE, all(missing >= 22L)])
  pc <- om_plan_refresh(lb, ed, de, akf, lead = BLASTAM_LEADIN_DAYS, n_vars = 3, max_n = 100)
  ok("the fetch-count cap is honoured too", pc[take == TRUE, .N] == 100L)
  # Under the old planner every one of the 150 cells past the tail cost 4.86.
  old_cost <- 150 * add_cost; new_cost <- pl[missing > 7, sum(cost)]
  ok("a September-type backlog costs a fraction of what it did", new_cost < old_cost / 2,
     sprintf("(%.0f vs %.0f weighted)", new_cost, old_cost))
  ok("an empty cache plans nothing, without error",
     nrow(om_plan_refresh(data.table(pid = character(), last = as.Date(character())),
                          ed, de, akf)) == 0L)
}

cat("\n18. An hourly HTTP 429 is waited out; a daily one stops the run\n")
# Regression: every 429 was "quota spent, stop". The 2026-09-07 grid run hit the
# HOURLY ceiling at exactly 5,000 weighted about an hour in, abandoned the rest
# of the grid with 140 minutes of deadline left, and the town fetch in the same
# hour was refused too: 31 towns of "no data" in that email.
{
  ok("the ceiling is read from the body",
     om_quota_kind("Hourly API request limit exceeded. Please try again in the next hour.") == "hour" &&
     om_quota_kind("Minutely API request limit exceeded. Please try again in one minute.") == "minute" &&
     om_quota_kind("Daily API request limit exceeded. Please try again tomorrow.") == "day" &&
     om_quota_kind("") == "unknown")
  # Stub the HTTP layer: fetch_points_batched() resolves om_request at call time.
  real_om_request <- om_request
  old_wait <- OM_QUOTA_WAIT_S
  hourly <- list(status = "quota", code = 429L, retry_after = NA_real_, body = NULL,
                 msg = "Hourly API request limit exceeded. Please try again in the next hour.")
  daily  <- modifyList(hourly, list(msg = "Daily API request limit exceeded. Please try again tomorrow."))
  http   <- list(status = "http", code = 503L, retry_after = NA_real_, body = NULL, msg = "boom")
  okbody <- function(lats) list(status = "ok", code = 200L, retry_after = NA_real_, msg = "",
    body = lapply(seq_along(lats), function(i) list(hourly = list(
      time = as.list(format(as.POSIXct("2026-06-01 00:00", tz = "UTC") + (0:23) * 3600,
                            "%Y-%m-%dT%H:%M")),
      temperature_2m = as.list(rep(25, 24)), relative_humidity_2m = as.list(rep(80, 24)),
      precipitation = as.list(rep(0, 24))))))
  pts <- data.table(pid = sprintf("t%02d", 1:30), lon = 145 + (1:30) / 10, lat = -25)
  on_pt <- function(pid, lon, lat, hw) data.table(pid = pid)
  script <- NULL; calls <- 0L
  om_request <- function(lats, lons, start_date, end_date, timeout_s = 60) {
    calls <<- calls + 1L
    r <- script[[min(calls, length(script))]]
    if (identical(r, "ok")) okbody(lats) else r
  }
  OM_QUOTA_WAIT_S <- 0.05          # read through .cfg() at call time
  paced <- 0L; count_pacer <- function(weight) { paced <<- paced + 1L; invisible(NULL) }
  d0 <- as.Date("2026-06-01"); d1 <- as.Date("2026-06-14")   # 14 days: 1.00 each

  script <- list(hourly, hourly, "ok")   # first batch refused twice, then served
  r <- fetch_points_batched(pts, d0, d1, on_pt, budget = Inf, label = "test-hourly",
                            pacer = count_pacer, quota_wait_max_s = 10)
  ok("an hourly 429 is waited out and the same batch is sent again",
     r$n_ok == 30L && r$stopped == "",
     sprintf("(ok %d, stopped '%s', %d requests)", r$n_ok, r$stopped, calls))
  ok("the wait is recorded", r$waited_s > 0)
  ok("refused requests are never charged", abs(r$spent - 30) < 1e-9, sprintf("(spent %.2f)", r$spent))

  calls <- 0L; script <- list(daily)
  r <- fetch_points_batched(pts, d0, d1, on_pt, budget = Inf, label = "test-daily",
                            pacer = count_pacer, quota_wait_max_s = 10)
  ok("a daily 429 stops the call at once", r$stopped == "quota" && r$n_ok == 0L && calls == 1L,
     sprintf("(stopped '%s', %d requests)", r$stopped, calls))
  ok("and the refused points are on the ledger as quota",
     nrow(r$ledger) == 25L && r$ledger[, all(status == "quota")])

  calls <- 0L; script <- list(hourly)    # refused every time
  r <- fetch_points_batched(pts, d0, d1, on_pt, budget = Inf, label = "test-cap",
                            pacer = count_pacer, quota_wait_max_s = 0.12)
  ok("the wait allowance bounds it", r$stopped == "quota" && r$waited_s <= 0.12 + 1e-9 && calls <= 4L,
     sprintf("(waited %.2f s over %d requests)", r$waited_s, calls))

  calls <- 0L; paced <- 0L; script <- list(http, "ok", "ok")
  r <- fetch_points_batched(pts, d0, d1, on_pt, budget = Inf, label = "test-retry",
                            pacer = count_pacer, quota_wait_max_s = 10)
  ok("every retry goes through the pacer", paced == calls && calls == 3L,
     sprintf("(%d requests, %d paced)", calls, paced))
  ok("and every non-429 attempt is charged", abs(r$spent - 55) < 1e-9, sprintf("(spent %.2f)", r$spent))

  om_request <- real_om_request; OM_QUOTA_WAIT_S <- old_wait
}

cat("\n19. Stale cells are drawn grey, not dropped, and the banner says why\n")
# Regression: with 76% of the grid current and 24% stale, the 90th-percentile
# fallback drew EVERY cell at a 16 day old window three weeks running (7, 14 and
# 21 Sep 2026), and the email blamed a daily quota that had not been spent.
{
  ed <- as.Date("2026-09-21")
  mk <- function(...) { v <- c(...); data.table(pid = sprintf("p%04d", seq_along(v)), mx = v) }
  # The September cache, scaled by ten: 5890 current; 286, 555 and 990 behind.
  sept <- mk(rep(ed, 589), rep(ed - 15L, 29), rep(ed - 22L, 55), rep(ed - 23L, 99))
  w <- grid_choose_window(sept, ed, "latest", min_cov = 0.9, draw_cov = 0.5, max_fallback = 3L)
  ok("the September cache is drawn at end_date", w$model_end == ed && w$rule == "draw-stale",
     sprintf("(got %s, rule %s)", format(w$model_end), w$rule))
  ok("with the stale quarter counted rather than hidden",
     abs(w$reach_now - 589 / 772) < 1e-9 && !nzchar(w$fallback_note))
  want <- sprintf("29 last updated %s, 55 %s, 99 %s", format(ed - 15L, "%d %b"),
                  format(ed - 22L, "%d %b"), format(ed - 23L, "%d %b"))
  ok("and its cohorts summarised newest first", grid_stale_summary(sept, ed) == want,
     sprintf("(got '%s')", grid_stale_summary(sept, ed)))
  # The pre-fix Monday: 89% current, 11% one day short. A one day step back that
  # brings every cell in beats greying out Western Australia.
  w <- grid_choose_window(mk(rep(ed, 89), rep(ed - 1L, 11)), ed, "latest", 0.9, 0.5, 3L)
  ok("a short step back that brings everyone in is still taken",
     w$model_end == ed - 1L && w$rule == "short-fallback" && grepl("fell back 1 day", w$fallback_note))
  # 7 Sep: 58% current, the rest two days back.
  w <- grid_choose_window(mk(rep(ed, 58), rep(ed - 2L, 42)), ed, "latest", 0.9, 0.5, 3L)
  ok("two days back is still a short step", w$model_end == ed - 2L && w$rule == "short-fallback")
  # Mostly stale: 30% current, 70% sixteen days back. The old complete map wins.
  w <- grid_choose_window(mk(rep(ed, 30), rep(ed - 16L, 70)), ed, "latest", 0.9, 0.5, 3L)
  ok("a mostly stale grid still falls back to the old complete map",
     w$model_end == ed - 16L && w$rule == "mostly-stale" && grepl("16 days", w$fallback_note))
  w <- grid_choose_window(mk(rep(ed, 100)), ed, "latest")
  ok("a healthy cache is untouched", w$model_end == ed && w$rule == "current" && !nzchar(w$fallback_note))
  ok("an empty cache does not error",
     grid_choose_window(mk(as.Date(character())), ed, "latest")$model_end == ed)
  ok("coverage mode is unchanged",
     grid_choose_window(mk(rep(ed, 95), rep(ed - 4L, 5)), ed, "coverage", coverage = 0.98)$model_end == ed - 4L)
  # The stated reason comes from what the fetch recorded, never a fixed sentence.
  ok("a 429 stop is named", grepl("HTTP 429", grid_fetch_reason(stops = "quota", spent = 5000, waited_s = 600)))
  ok("an unaffordable backlog is named",
     grepl("did not fit", grid_fetch_reason(n_left = 990, left_cost = 2050, plan_cap = 8550)))
  ok("a ledger cap from an earlier run is named",
     grepl("earlier run today", grid_fetch_reason(already = 8700, wt_cap = 800)))
  ok("cells still behind with nothing recorded is said plainly",
     grepl("incomplete or empty", grid_fetch_reason(held_out = 5)))
  ok("nothing to explain gives an empty string", grid_fetch_reason() == "")
  ok("the email no longer asserts a spent quota it cannot know about",
     !grepl("usually because the daily weather-API quota was already spent", tsrc, fixed = TRUE))
  ok("the renderer is handed the stale cells", grepl("stale_pts = stale_pts", gsrc, fixed = TRUE))
  ok("and masks them out of the value raster before the GeoTIFF is written",
     regexpr("terra::mask(r, sm, inverse = TRUE)", gsrc, fixed = TRUE) <
     regexpr("WRITE_GEOTIFF))", gsrc, fixed = TRUE))
}

cat("\n20. A degraded run is called degraded: verdict, status file, subject, workflow\n")
# Regression: the 2026-09-07 run delivered 31 towns of "no data" and the 09-14
# and 09-21 runs a map three weeks stale; all three exited 0, showed a green tick
# and went out under a subject that said nothing was wrong.
source("run_health.R")
{
  v <- health_verdict(0L, 31L)
  ok("31 towns of no data is degraded", v$degraded && grepl("only 0 of 31 towns", v$reasons[1]))
  ok("29 of 31 towns (94%) is not", !health_verdict(29L, 31L)$degraded)
  ok("27 of 31 towns is", health_verdict(27L, 31L)$degraded)
  v <- health_verdict(31L, 31L, map_behind_days = 16L)
  ok("a map 16 days behind the table is degraded", v$degraded && grepl("16 days behind", v$reasons))
  ok("one day behind is tolerated", !health_verdict(31L, 31L, map_behind_days = 1L)$degraded)
  v <- health_verdict(31L, 31L, map_cells = 5890L, map_grey = 1831L)
  ok("a quarter of the grid grey is degraded", v$degraded && grepl("1831 map cells \\(24%\\)", v$reasons))
  v <- health_verdict(31L, 31L, map_cells = 7700L, map_grey = 21L)
  ok("a few grey cells are a note, not a failure", !v$degraded && length(v$warnings) == 1L)
  v <- health_verdict(31L, 31L, map_cells = 5000L, map_cells_prev = 7721L)
  ok("a 35% drop in mapped cells is degraded", v$degraded && grepl("fell from 7721 to 5000", v$reasons))
  ok("growth is not", !health_verdict(31L, 31L, map_cells = 7721L, map_cells_prev = 5000L)$degraded)
  v <- health_verdict(31L, 31L, fetch_reason = "the weather API refused further requests (HTTP 429)")
  ok("a fetch reason alone is a note, not a failure", !v$degraded && grepl("HTTP 429", v$warnings))
  ok("missing map stats neither error nor fail",
     !health_verdict(31L, 31L, map_cells = NA, map_cells_prev = NA, map_grey = NA)$degraded)
  f <- tempfile(fileext = ".txt")
  writeLines(run_status_lines(list(run_date = "2026-09-21", degraded = 1L,
                                   reasons = c("a | b", "c\nd"), warnings = character(0))), f)
  st <- read_run_status(f)
  ok("the status file round-trips and is pipe- and newline-safe",
     st$run_date == "2026-09-21" && st$degraded == "1" && st$reasons == "a   b; c d" &&
     identical(st$warnings, ""), sprintf("(reasons '%s')", st$reasons))
  unlink(f)
  ok("run_blast.R judges the run and writes the status file",
     grepl("health_verdict(", tsrc, fixed = TRUE) && grepl("RUN_STATUS_FILE", tsrc, fixed = TRUE))
  py <- paste(readLines("send_email.py", warn = FALSE), collapse = "\n")
  ok("send_email.py reads the status and marks the subject",
     grepl("run_status.txt", py, fixed = TRUE) && grepl("[DEGRADED]", py, fixed = TRUE))
  ok("and puts the town window first, adding the map window only when it differs",
     grepl("town_window_end", py, fixed = TRUE) && grepl("maps to", py, fixed = TRUE))
  yml <- paste(readLines(".github/workflows/weekly_blast.yml", warn = FALSE), collapse = "\n")
  ok("the workflow ends with a Verdict step that fails on degraded=1",
     grepl("name: Verdict", yml, fixed = TRUE) && grepl("degraded=1", yml, fixed = TRUE))
  ok("which runs after the email and the artifact upload",
     regexpr("name: Verdict", yml, fixed = TRUE) > regexpr("name: Upload artifact", yml, fixed = TRUE) &&
     regexpr("name: Verdict", yml, fixed = TRUE) > regexpr("name: Email summary", yml, fixed = TRUE))
  ok("and clears the previous run's status at the start",
     grepl("rm -f blast_outputs/run_status.txt", yml, fixed = TRUE))
  ok("and commits the status file with the other state",
     grepl("run_status.txt", sub("Email summary.*$", "", sub("^.*Commit trends and run state", "", yml)), fixed = TRUE))
}

cat("\n21. A rerun over the same window cannot make the record worse\n")
# Regression: on 2026-09-21 a manual rerun the same UTC day, correctly capped by
# the spend ledger at 42 weighted calls for the towns, fell through to a serial
# fallback that was unpaced, uncharged and invisible to the ledger, fetched 29
# towns off the books, and then overwrote the scheduled run's 31-town trends
# column and run log row with its own (Moree and Borroloola became NA).
{
  # (a) the second pass goes through fetch_points_batched(), one town per request
  real_om_request <- om_request; calls <- 0L
  om_request <- function(lats, lons, start_date, end_date, timeout_s = 60) {
    calls <<- calls + 1L
    list(status = "ok", code = 200L, retry_after = NA_real_, msg = "",
         body = lapply(seq_along(lats), function(i) list(hourly = list(
           time = as.list(format(as.POSIXct("2026-06-01 00:00", tz = "UTC") + (0:23) * 3600,
                                 "%Y-%m-%dT%H:%M")),
           temperature_2m = as.list(rep(25, 24)), relative_humidity_2m = as.list(rep(80, 24)),
           precipitation = as.list(rep(0, 24))))))
  }
  pts5 <- data.table(pid = sprintf("s%d", 1:5), lon = 145 + (1:5) / 10, lat = -25)
  one <- function(pid, lon, lat, hw) data.table(pid = pid)
  r <- fetch_points_batched(pts5, as.Date("2026-06-01"), as.Date("2026-06-14"), one,
                            budget = Inf, label = "test-single", batch_size = 1L,
                            pacer = function(w) invisible(NULL))
  ok("batch_size = 1 sends one request per town", calls == 5L && r$n_ok == 5L, sprintf("(%d requests)", calls))
  r <- fetch_points_batched(pts5, as.Date("2026-06-01"), as.Date("2026-06-14"), one,
                            budget = 3, label = "test-single-budget", batch_size = 1L,
                            pacer = function(w) invisible(NULL))
  ok("and it stops at the budget like any other fetch",
     r$n_ok == 3L && r$stopped == "budget" && abs(r$spent - 3) < 1e-9)
  om_request <- real_om_request
  ok("run_blast.R no longer fetches towns outside the batch path",
     !grepl("get_openmeteo_hourly", tsrc, fixed = TRUE) && grepl("batch_size = 1L", tsrc, fixed = TRUE))
  # (b) trends: this run's blanks are filled from the earlier run's column
  hist <- data.table(town = c("Moree", "Borroloola", "Dubbo"),
                     `2026-09-07` = c(1, 2, 3), `2026-09-14` = c(0.0009, 3, 0.0059))
  today <- data.table(town = c("Moree", "Borroloola", "Dubbo")); today[["2026-09-14"]] <- c(NA, NA, 0.0061)
  m <- trends_merge_rerun(today, hist, "2026-09-14")
  ok("towns the rerun could not model keep the earlier value",
     m$n_kept == 2L && m$today[town == "Moree"][["2026-09-14"]] == 0.0009 &&
     m$today[town == "Borroloola"][["2026-09-14"]] == 3)
  ok("towns it did model take the new value", m$today[town == "Dubbo"][["2026-09-14"]] == 0.0061)
  ok("the earlier column is dropped so the merged one replaces it",
     !"2026-09-14" %in% names(m$hist) && "2026-09-07" %in% names(m$hist))
  m2 <- trends_merge_rerun(today, hist[, .(town, `2026-09-07`)], "2026-09-14")
  ok("a first run over a window merges nothing",
     m2$n_kept == 0L && identical(m2$today[["2026-09-14"]], c(NA, NA, 0.0061)))
  # (c) run log: the row that modelled the most towns is the one kept
  old <- data.table(run_date = c("2026-09-14", "2026-09-21"),
                    data_end = c("2026-09-07", "2026-09-14"), towns_modelled = c(31L, 31L))
  k <- runlog_keep_better(old, data.table(run_date = "2026-09-21", data_end = "2026-09-14",
                                          towns_modelled = 29L), "2026-09-14")
  ok("a rerun that modelled fewer towns does not replace the row", !k$write_new && k$prev_towns == 31L)
  k <- runlog_keep_better(old, data.table(run_date = "2026-09-21", data_end = "2026-09-14",
                                          towns_modelled = 31L), "2026-09-14")
  ok("an equal or better rerun does, and the old row goes", k$write_new && nrow(k$old) == 1L)
  k <- runlog_keep_better(old, data.table(run_date = "2026-09-28", data_end = "2026-09-21",
                                          towns_modelled = 5L), "2026-09-21")
  ok("a new window is always written", k$write_new && nrow(k$old) == 2L)
  ok("an empty log is always written", runlog_keep_better(NULL, data.table(towns_modelled = 1L), "x")$write_new)
  # (d) the spend ledger is keyed on the run's pinned UTC date, not the clock: the
  # Sunday grid run ends after 00:00 UTC, and booking it to Monday starved the
  # 2026-09-21 rerun of budget the API had mostly counted against Sunday.
  oldu <- Sys.getenv("BLAST_UTC_DATE", unset = NA)
  Sys.setenv(BLAST_UTC_DATE = "2026-09-20")
  ok("the spend ledger day is the pinned UTC date", om_spend_utc_day() == as.Date("2026-09-20"),
     sprintf("(got %s)", format(om_spend_utc_day())))
  lf <- tempfile(fileext = ".csv"); om_spend_add(lf, 8550, "grid")
  Sys.setenv(BLAST_UTC_DATE = "2026-09-21")
  ok("so a run on the next pinned day starts with a clean budget", om_spend_read(lf) == 0)
  Sys.setenv(BLAST_UTC_DATE = "2026-09-20")
  ok("and a rerun on the same pinned day still shares it", om_spend_read(lf) == 8550)
  unlink(lf); if (is.na(oldu)) Sys.unsetenv("BLAST_UTC_DATE") else Sys.setenv(BLAST_UTC_DATE = oldu)
}

cat("\n22. The midweek top-up: a second fetch day, no maps, no email\n")
# The grid takes about twelve weekly runs to fill from cold and, after an
# interrupted run, several to recover; a second fetch day halves both. A
# BLAST_MIDWEEK branch existed in run_blast_grid.R but ran the whole pipeline,
# nothing invoked it, and the email line assumed a daily job that never existed.
{
  ok("MIDWEEK_MAX_AGE_DAYS is set so that a missed Sunday top-up is reported on Monday",
     exists("MIDWEEK_MAX_AGE_DAYS") && MIDWEEK_MAX_AGE_DAYS >= 1L && MIDWEEK_MAX_AGE_DAYS <= 3L)
  ok("in midweek mode the grid runner saves the cache and stops before modelling",
     regexpr("quit(save = \"no\", status = 0)", gsrc, fixed = TRUE) > 0 &&
     regexpr("quit(save = \"no\", status = 0)", gsrc, fixed = TRUE) <
     regexpr("writeLines(run_tag, file.path(OUT, \"run_date.txt\"))", gsrc, fixed = TRUE))
  ok("the same cache writer serves both paths",
     lengths(regmatches(gsrc, gregexpr("save_cache_now(cache)", gsrc, fixed = TRUE))) == 2L)
  ymlm <- if (file.exists(".github/workflows/midweek_topup.yml"))
    paste(readLines(".github/workflows/midweek_topup.yml", warn = FALSE), collapse = "\n") else ""
  ok("a midweek workflow exists", nzchar(ymlm))
  ok("it sets BLAST_MIDWEEK and runs the grid script",
     grepl("BLAST_MIDWEEK", ymlm, fixed = TRUE) && grepl("run_blast_grid.R", ymlm, fixed = TRUE))
  ok("it shares the weekly run's concurrency group", grepl("group: blast-grid", ymlm, fixed = TRUE))
  ok("it sends no email and runs no town table",
     !grepl("send_email.py", ymlm, fixed = TRUE) && !grepl("run_blast.R", ymlm, fixed = TRUE))
  ok("it saves the cache and commits the top-up status",
     grepl("weather_cache.csv.gz", ymlm, fixed = TRUE) && grepl("midweek_status.txt", ymlm, fixed = TRUE))
  # Thursday recovers from a bad Monday; Sunday is Monday's insurance (cells the
  # Monday fetch does not reach are then one day behind, not three). Neither may
  # share the weekly run's UTC day, which is Monday, or its daily quota.
  tl <- strsplit(ymlm, "\n")[[1]]
  tcr <- unlist(regmatches(tl, gregexpr("(?<=cron: ')[^']+", tl, perl = TRUE)))
  tday <- vapply(strsplit(tcr, " "), function(x) x[5], character(1))
  thr  <- vapply(strsplit(tcr, " "), function(x) as.integer(x[2]), integer(1))
  ok("there are two top-ups a week, Thursday and Sunday UTC", setequal(tday, c("4", "0")),
     sprintf("(crons: %s)", paste(tcr, collapse = " | ")))
  ok("neither on the weekly run's UTC day", !"1" %in% tday)
  ok("both after the archive's daily update, so they fetch the newest day", all(thr >= 2L))
  ok("the Monday email reports the top-up on a weekly cadence",
     grepl("MIDWEEK_MAX_AGE_DAYS", tsrc, fixed = TRUE) && !grepl("Daily top-up", tsrc, fixed = TRUE))
}

cat("\n23. The email sender runs end to end, minus the SMTP conversation\n")
# Regression: the 2026-09-28 email was SENT and the step still failed. A variable
# renamed in main() (window_end to map_end) was left behind in the final "Email
# sent" print, so send_email.py died with a NameError one line after
# server.send_message(). Nothing had executed the script before the live run: the
# suite only grepped it. BLAST_EMAIL_DRY_RUN=1 runs every line except the SMTP
# block, here against fixture files and in CI before any fetching starts.
{
  py_src <- paste(readLines("send_email.py", warn = FALSE), collapse = "\n")
  py_code <- paste(sub("#.*$", "", readLines("send_email.py", warn = FALSE)), collapse = "\n")
  ok("no code path refers to the removed window_end variable",
     !grepl("(?<![A-Za-z_])window_end", py_code, perl = TRUE))
  ok("the sender has a dry-run mode and an output-directory override",
     grepl("BLAST_EMAIL_DRY_RUN", py_src, fixed = TRUE) && grepl("BLAST_OUT_DIR", py_src, fixed = TRUE))
  ok("the weekly workflow installs python before the tests, so the dry run is not skipped there",
     regexpr("name: Ensure python3", yml, fixed = TRUE) > 0 &&
     regexpr("name: Ensure python3", yml, fixed = TRUE) < regexpr("name: Offline tests", yml, fixed = TRUE))
  # Take the first interpreter that actually RUNS. On Windows `python3` is often
  # the Microsoft Store stub, which exists on PATH and exits non-zero.
  cands <- Sys.which(c("python3", "python", "py")); cands <- unname(cands[nzchar(cands)])
  runs <- function(p) isTRUE(tryCatch(
    suppressWarnings(system2(p, "--version", stdout = FALSE, stderr = FALSE)) == 0L,
    error = function(e) FALSE))
  py <- NA_character_
  for (p in cands) if (runs(p)) { py <- p; break }
  if (is.na(py)) {
    cat("  SKIP  no working python here; the dry run is exercised in CI\n")
  } else {
    fx <- tempfile("blast_email_"); dir.create(fx)
    writeLines("plain body", file.path(fx, "blast_summary_latest.txt"))
    writeLines("<p>html body</p>", file.path(fx, "blast_summary_latest.html"))
    writeLines(c("town,2026-09-20", "Dubbo,0"), file.path(fx, "town_trends.csv"))
    writeLines(c("town,2026-09-20", "Dubbo,0"), file.path(fx, "blastam_trends.csv"))
    writeLines("7272|0.31|7721|0.30|gz|6496|gz|2026-08-29|0.60|8550|0.0164|8||epirice+blastam|449|x|y",
               file.path(fx, "map_stats.txt"))
    run_dry <- function(status_lines) {
      sf <- file.path(fx, "run_status.txt")
      if (is.null(status_lines)) unlink(sf) else writeLines(status_lines, sf)
      vars <- c("MAIL_USERNAME", "MAIL_PASSWORD", "MAIL_TO", "BLAST_RUN_DATE", "RUN_DATE",
                "BLAST_OUT_DIR", "BLAST_EMAIL_DRY_RUN")
      old <- Sys.getenv(vars, unset = NA)
      Sys.setenv(MAIL_USERNAME = "sender@example.org", MAIL_PASSWORD = "x",
                 MAIL_TO = "a@example.org, b@example.org", BLAST_RUN_DATE = "2026-09-21",
                 BLAST_OUT_DIR = fx, BLAST_EMAIL_DRY_RUN = "1")
      on.exit(for (v in names(old))
        if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, as.list(setNames(old[[v]], v))),
        add = TRUE)
      out <- suppressWarnings(system2(py, "send_email.py", stdout = TRUE, stderr = TRUE))
      st <- attr(out, "status")
      list(status = if (is.null(st)) 0L else as.integer(st), out = paste(out, collapse = "\n"))
    }
    r <- run_dry(c("town_window_end=2026-09-14", "degraded=1",
                   "reasons=1831 map cells (24%) not refreshed, drawn grey"))
    ok("a degraded dry run exits 0 and reaches the final line",
       r$status == 0L && grepl("DRY RUN", r$out, fixed = TRUE),
       sprintf("(status %d: %s)", r$status, substr(r$out, 1, 400)))
    ok("with the subject September should have had",
       grepl("[DEGRADED] Blast risk summary 2026-09-21 (weather to 2026-09-14; maps to 2026-08-29)",
             r$out, fixed = TRUE))
    ok("and the attachments counted", grepl("with 2 attachments", r$out, fixed = TRUE))
    r <- run_dry(c("town_window_end=2026-08-29", "degraded=0", "reasons="))
    ok("a healthy dry run has a plain subject carrying one window",
       r$status == 0L && grepl("subject 'Blast risk summary 2026-09-21 (weather to 2026-08-29)'",
                                r$out, fixed = TRUE), sprintf("(%s)", substr(r$out, 1, 300)))
    r <- run_dry(NULL)
    ok("without a status file it falls back to the map window",
       r$status == 0L && grepl("(weather to 2026-08-29)", r$out, fixed = TRUE) &&
       !grepl("[DEGRADED]", r$out, fixed = TRUE))
    unlink(fx, recursive = TRUE)
  }
}

cat("\n24. The weekly run fires once, on Monday UTC\n")
# History. Two seasonal crons at 06:30 local, with a gate job choosing between
# them by daylight-saving offset; GitHub fired them up to 2 h 40 late, so they
# moved to 01:47. Both were SUNDAY in UTC, which capped the map at eight days old.
# The archive publishes one more day at about 00:30 to 01:00 UTC, so a run on
# MONDAY UTC models a day later. That constraint is in UTC: one cron, no gate.
{
  wl <- readLines(".github/workflows/weekly_blast.yml", warn = FALSE)
  crons <- unlist(regmatches(wl, gregexpr("(?<=cron: ')[^']+", wl, perl = TRUE)))
  f <- function(s) strsplit(s, " ")[[1]]
  ok("there is exactly one weekly cron", length(crons) == 1L,
     sprintf("(found: %s)", paste(crons, collapse = " | ")))
  ok("on Monday in UTC, so the UTC date equals the Sydney run date",
     length(crons) == 1L && f(crons[1])[5] == "1")
  ok("early in the UTC day, so the email is still Monday afternoon in Sydney",
     length(crons) == 1L && as.integer(f(crons[1])[2]) <= 2L)
  ok("off the hour and half hour, where the scheduler is busiest",
     length(crons) == 1L && !as.integer(f(crons[1])[1]) %in% c(0L, 30L))
  ok("no seasonal gate is left to disagree with it",
     !any(grepl("$SCHED\" = \"", wl, fixed = TRUE)) && !any(grepl("needs: gate", wl, fixed = TRUE)))
  vars <- c("BLAST_RUN_DATE", "BLAST_UTC_DATE", "BLAST_DATA_END"); old <- Sys.getenv(vars, unset = NA)
  Sys.setenv(BLAST_RUN_DATE = "2026-10-05", BLAST_UTC_DATE = "2026-10-05"); Sys.unsetenv("BLAST_DATA_END")
  ok("a Monday-UTC run models to the previous Monday, a day later than a Sunday-UTC one",
     blast_data_end(blast_run_date()) - DAY_CUT_LAG_DAYS == as.Date("2026-09-28"),
     sprintf("(got %s)", format(blast_data_end(blast_run_date()) - DAY_CUT_LAG_DAYS)))
  for (v in vars) if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, as.list(setNames(old[[v]], v)))
}

cat("\n25. The weather cache is not versioned on main\n")
# Every run committed a new 6.5 MB cache to main, and git keeps every version of
# a committed file: after 38 runs the old caches came to 104 MB of the 110 MB
# repository, and three runs a week would add about a gigabyte a year. No model
# reads weather older than the crop window. The cache now lives on `cache-data`
# as ONE parentless commit that each run replaces.
{
  for (wf in c("weekly_blast.yml", "midweek_topup.yml")) {
    y <- paste(readLines(file.path(".github/workflows", wf), warn = FALSE), collapse = "\n")
    pos <- function(s) regexpr(s, y, fixed = TRUE)
    ok(sprintf("%s restores the cache before the grid script and saves it straight after", wf),
       pos("name: Restore weather cache") > 0 &&
       pos("name: Restore weather cache") < pos("run: Rscript run_blast_grid.R") &&
       pos("run: Rscript run_blast_grid.R") < pos("name: Save weather cache"))
    ok(sprintf("%s refuses to run on an empty cache", wf),
       grepl("Refusing to run on an empty cache", y, fixed = TRUE))
    ok(sprintf("%s writes a parentless commit and force-pushes cache-data, nothing else", wf),
       grepl("commit-tree", y, fixed = TRUE) &&
       lengths(regmatches(y, gregexpr("push --force", y, fixed = TRUE))) == 1L &&
       grepl("push --force origin \"$commit:refs/heads/cache-data\"", y, fixed = TRUE))
    ok(sprintf("%s will not replace the cache with one under half its size", wf),
       grepl("Refusing to replace it", y, fixed = TRUE))
    commit_step <- sub("\n      - name:.*$", "", sub("^.*\n      - name: Commit ", "", y))
    ok(sprintf("%s no longer adds the cache to main", wf),
       !grepl("weather_cache", commit_step, fixed = TRUE) && !grepl("cache_version", commit_step, fixed = TRUE))
  }
  gi <- readLines(".gitignore", warn = FALSE)
  ok("the cache and its schema marker are ignored on main",
     !any(grepl("^!blast_outputs/(weather_cache|cache_version)", gi)))
  tracked <- tryCatch(suppressWarnings(system2("git", c("ls-files", "blast_outputs"), stdout = TRUE, stderr = FALSE)),
                      error = function(e) character(0))
  if (length(tracked) == 0L) cat("  SKIP  not a git checkout; cannot check what is tracked\n") else
    ok("and no longer tracked there", !any(grepl("weather_cache|cache_version", tracked)),
       sprintf("(tracked: %s)", paste(grep("cache", tracked, value = TRUE), collapse = ", ")))
  need <- CROP_AGE_DAYS + 1L + GRID_WINDOW_MAX_LAG_DAYS
  ok("the cache keeps what the models need, with room for a fallback window",
     CACHE_HISTORY_DAYS >= need + 14L, sprintf("(%d days kept, %d needed)", CACHE_HISTORY_DAYS, need))
  ok("and no more than about three months", CACHE_HISTORY_DAYS <= 93L)
}

cat("\n26. The run asks the archive what it has, rather than assuming\n")
# Regression: the last day fetched was arithmetic, "UTC date minus 6". Open-Meteo
# publishes that day at about 00:30 to 01:00 UTC, so earlier in a UTC day it is
# not there, and every cell west of 120 E, which needs three hours of it, lands a
# day short (test 16). It was still missing at 00:08, 00:11 and 00:16 UTC on
# 2 October 2026; a run that began at 00:26 UTC on 25 August caught it part way.
{
  mkh <- function(last_full, next_day_hours = integer(0)) {
    tt <- seq(as.POSIXct("2026-09-20 00:00", tz = "UTC"),
              as.POSIXct(paste(as.Date(last_full) + 1L, "23:00"), tz = "UTC"), by = "hour")
    d <- data.table(dt = tt, temp = 20, rh = 60, rain = 0)
    nxt  <- as.Date(format(d$dt, "%Y-%m-%d", tz = "UTC")) == as.Date(last_full) + 1L
    keep <- as.integer(format(d$dt, "%H", tz = "UTC")) %in% next_day_hours
    d[nxt & !keep, `:=`(temp = NA_real_, rh = NA_real_, rain = NA_real_)]
    d
  }
  ok("an unpublished last day is not counted",
     om_archive_edge_from_hourly(mkh("2026-09-25"), 3L) == as.Date("2026-09-25"))
  ok("a day counts once the hours the westernmost cell needs are there",
     om_archive_edge_from_hourly(mkh("2026-09-25", 0:2), 3L) == as.Date("2026-09-26"))
  ok("two of the three hours is not enough",
     om_archive_edge_from_hourly(mkh("2026-09-25", 0:1), 3L) == as.Date("2026-09-25"))
  ok("no data at all gives NA, not an error", is.na(om_archive_edge_from_hourly(NULL, 3L)))
  ok("the hours needed follow from the day cut and the western edge of the grid",
     ceiling(BLASTAM_DAY_CUT_HOUR - GRID_EXTENT[1] / 15) == 3)
  vars <- c("BLAST_RUN_DATE", "BLAST_UTC_DATE", "BLAST_DATA_END"); old <- Sys.getenv(vars, unset = NA)
  Sys.setenv(BLAST_RUN_DATE = "2026-10-05", BLAST_UTC_DATE = "2026-10-05")   # arithmetic edge 2026-09-29
  Sys.setenv(BLAST_DATA_END = "2026-09-28")
  ok("a probed edge behind the arithmetic one is used", blast_data_end(blast_run_date()) == as.Date("2026-09-28"))
  Sys.setenv(BLAST_DATA_END = "2026-10-02")
  ok("a probed edge ahead of it is ignored", blast_data_end(blast_run_date()) == as.Date("2026-09-29"))
  Sys.setenv(BLAST_DATA_END = "rubbish")
  ok("and so is one that is not a date", blast_data_end(blast_run_date()) == as.Date("2026-09-29"))
  for (v in vars) if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, as.list(setNames(old[[v]], v)))
  wk <- paste(readLines(".github/workflows/weekly_blast.yml", warn = FALSE), collapse = "\n")
  tu <- paste(readLines(".github/workflows/midweek_topup.yml", warn = FALSE), collapse = "\n")
  for (y in list(wk, tu))
    ok("the workflow probes after the tests and before the grid script",
       regexpr("name: Offline tests", y, fixed = TRUE) < regexpr("name: Resolve archive edge", y, fixed = TRUE) &&
       regexpr("name: Resolve archive edge", y, fixed = TRUE) < regexpr("run: Rscript run_blast_grid.R", y, fixed = TRUE))
  ok("the weekly run waits for the daily update; a top-up takes what is there",
     !grepl("BLAST_EDGE_WAIT_MIN", wk, fixed = TRUE) && grepl("BLAST_EDGE_WAIT_MIN: \"0\"", tu, fixed = TRUE))
  ok("and the weekly job's timeout allows for that wait",
     as.integer(sub(".*timeout-minutes: ([0-9]+).*", "\\1", wk)) >=
       ARCHIVE_EDGE_WAIT_MAX_MIN + GRID_MAX_MINUTES + GRID_RESERVE_MINUTES + TOWN_QUOTA_WAIT_MAX_MIN)
}

cat(sprintf("\n%d tests, %d failures\n", n, fails))
quit(status = if (fails > 0L) 1L else 0L)

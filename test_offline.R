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

cat("\n13. Grid lattice covers the requested extent\n")
fin <- GRID_RES_FINEST
lat_top <- GRID_EXTENT[3] + ceiling((GRID_EXTENT[4] - GRID_EXTENT[3]) / fin) * fin
ok("the northernmost row is inside the lattice", lat_top >= GRID_EXTENT[4],
   sprintf("(lattice reaches %.2f, extent asks for %.2f)", lat_top, GRID_EXTENT[4]))
ok("the old seq() would have dropped it",
   max(seq(GRID_EXTENT[3], GRID_EXTENT[4], by = fin)) < GRID_EXTENT[4] ||
   isTRUE(all.equal((GRID_EXTENT[4] - GRID_EXTENT[3]) %% fin, 0)))

cat("\n14. Overlay artefacts and label declutter (terra)\n")
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

cat(sprintf("\n%d tests, %d failures\n", n, fails))
quit(status = if (fails > 0L) 1L else 0L)

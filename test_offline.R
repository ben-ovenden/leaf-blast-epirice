#!/usr/bin/env Rscript
################################################################################
# tests/test_offline.R
#
# Offline regression tests. No network, no API quota. Run from the repo root:
#
#   Rscript tests/test_offline.R
#
# Each test guards a bug that was actually shipped at some point, so a failure
# here means a real regression rather than a style complaint.
################################################################################
suppressPackageStartupMessages({library(data.table)})
SCRIPT_DIR <- normalizePath(".", winslash = "/")
source("blast_config.R"); source("epirice_model.R"); source("blastam_model.R")

fails <- 0L; n <- 0L
ok <- function(label, cond, extra = "") {
  n <<- n + 1L
  if (isTRUE(cond)) cat(sprintf("  PASS  %s\n", label))
  else { fails <<- fails + 1L; cat(sprintf("  FAIL  %s %s\n", label, extra)) }
}

# Synthetic hourly UTC series whose humidity follows a LOCAL diurnal cycle.
mkseries <- function(lon, days = 20, wet_from = 19, wet_to = 8, tbase = 24) {
  dt <- as.POSIXct("2026-06-01 00:00", tz = "UTC") + (0:(days * 24 - 1)) * 3600
  loc <- (as.numeric(format(dt, "%H")) + lon / 15) %% 24
  data.table(dt = dt,
             temp = tbase + 4 * cos((loc - 15) / 24 * 2 * pi),
             rh   = ifelse(loc >= wet_from | loc < wet_to, 95, 55),
             rain = 0)
}

cat("\n1. BLASTAM night window is LOCAL SOLAR, not UTC\n")
# Regression: applying the 15:00-09:00 window to UTC timestamps put it at roughly
# 01:00-19:00 local over Australia, returning 8 wet hours instead of 13 and zero
# favourable nights.
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

cat("\n2. Unjudgeable nights are NA, not 0\n")
# Regression: returning 0 asserted "not favourable" when the truth was "not yet
# known", and a short refresh then wrote those zeros over good cached values.
d <- blastam_daily_from_hourly(mkseries(145), lon = 145)
ok("first 4 days NA (no preceding 5-day mean)", all(is.na(head(d$infect, 4))))
ok("lead-in constant covers it", BLASTAM_LEADIN_DAYS >= 5L)
ok("short wetness gives 0, not NA", {
  s <- blastam_daily_from_hourly(mkseries(145, wet_from = 23, wet_to = 3), lon = 145)
  sj <- s[!is.na(infect)]; nrow(sj) > 0 && all(sj$infect == 0L)
})

cat("\n3. Partial local days are dropped\n")
ok("no day carries a truncated aggregate",
   all(diff(as.integer(d$date)) == 1L) && nrow(d) == 19L,
   sprintf("(got %d rows)", nrow(d)))

cat("\n4. EPIRICE date alignment\n")
# Regression: filtering the modelling window with `>` instead of `>=` gave
# CROP_AGE_DAYS rows starting one day AFTER emergence, and SEIR's alignment check
# then threw for every point, silently emptying the EPIRICE map.
end_date <- as.Date("2026-07-22"); emergence <- end_date - CROP_AGE_DAYS
w <- data.table(YYYYMMDD = seq(emergence, end_date, by = "day"))
w[, `:=`(DOY = as.integer(format(YYYYMMDD, "%j")), TEMP = 24, RHUM = 92,
         RAIN = 1, LAT = -25, LON = 148)]
ok("inclusive window models without error",
   !is.null(tryCatch(predict_leaf_blast(w, emergence, nrow(w)), error = function(e) NULL)))
ok("exclusive window is the failure mode this guards",
   is.null(tryCatch(predict_leaf_blast(w[-1], emergence, nrow(w) - 1L),
                    error = function(e) NULL)))

cat("\n5. Weighted cost model\n")
source("openmeteo_batch.R")
# Guard rather than let a missing symbol halt the suite: a swapped-out or stubbed
# openmeteo_batch.R should report a clear failure, not an "Execution halted" with
# no summary line.
have <- function(f) exists(f, mode = "function")
if (!have("om_weight_per_location")) {
  ok("om_weight_per_location() is defined", FALSE, "(openmeteo_batch.R stubbed?)")
} else {
  ok("14 day floor: 1 day costs the same as 14",
     om_weight_per_location(1, 3) == om_weight_per_location(14, 3))
  ok("61 days is ~4.36", abs(om_weight_per_location(61, 3) - 61/14) < 1e-9)
  ok("11+ variables trigger the multiplier", om_weight_per_location(14, 11) > 1)
}

cat("\n6. Pacer holds the sustained rate\n")
if (!have("om_pacer")) {
  ok("om_pacer() is defined", FALSE, "(openmeteo_batch.R stubbed?)")
} else {
  p <- om_pacer(6000); p(6000)
  t0 <- Sys.time(); for (i in 1:5) p(100)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # Loose bound: shared CI runners are not real-time systems.
  ok("500 weighted at 6000/min takes ~5 s", abs(el - 5) < 2.0, sprintf("(got %.2fs)", el))
}

cat("\n7. Cache gz is really gzipped\n")
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

cat(sprintf("\n%d tests, %d failures\n", n, fails))
quit(status = if (fails > 0L) 1L else 0L)

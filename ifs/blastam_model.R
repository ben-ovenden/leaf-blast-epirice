################################################################################
# blastam_model.R  -  BLASTAM leaf blast infection-warning model
#
# WHAT IT IS
# BLASTAM (Koshimizu 1988) is the Japanese infection-warning model run on the
# AMeDAS automated weather network. For each night it judges whether the
# leaf-wetness period and temperatures were favourable for a Magnaporthe oryzae
# (leaf blast) infection event. It does NOT simulate the epidemic; it counts
# favourable infection days over a rolling window. This is fundamentally
# different from EPIRICE, which accumulates a disease progress curve.
#
# blastam_daily_from_hourly() serves BOTH models from one hourly fetch, so it
# also produces the daily TEMP, RHUM and RAIN that EPIRICE reads.
#
# EVERY PARAMETER IS DEFINED IN blast_config.R. The assignments below are
# fallbacks for the case where this file is sourced on its own. They are guarded
# with if (!exists(...)), because previously they were unconditional and, since
# blast_config.R is sourced FIRST, they silently overwrote anything set there.
# The README's advice to set BLASTAM_TWET_MAX in blast_config.R did nothing.
#
# ==========================================================================
# ORIGINAL CRITERIA (Koshimizu 1988, Hayashi & Koshimizu 1988)
# ==========================================================================
# A night is FAVOURABLE when ALL THREE hold:
#   1. leaf wetness duration >= threshold (see below)
#   2. mean air temperature DURING the wetness period: 15-25 C
#   3. mean air temperature of the PRECEDING 5 DAYS:  20-25 C
# If wetness threshold is met but exactly ONE temperature condition fails, the
# night is SEMI-FAVOURABLE.
#
# Night window: 15:00 to 09:00 the following day, local time (Japan Standard
# Time for AMeDAS). Leaf wetness estimated from an energy balance on AMeDAS
# sunshine, wind and rainfall, with hours of >= 4 mm/h rain excluded (Yoshino
# 1988) because heavy rain washes conidia off the leaf.
#
# ==========================================================================
# DEVIATIONS USED HERE, AND RATIONALE
# ==========================================================================
#
# DEVIATION 1. Temperature-dependent wetness threshold (replaces the fixed 10 h)
#   Barksdale & Jones (1965) minimum infection curve, lower 95% CI: 12.2 h at
#   15.6 C falling to 7.7 h at 26.7 C. The fixed 10 h was calibrated for
#   temperate Tohoku, where wetness-period temperatures are 15 to 20 C. In
#   tropical Australia it over-counts; in the cool dry season it under-counts.
#   Restore with BLASTAM_USE_BJ_THRESHOLD <- FALSE.
#
# DEVIATION 2. Upper temperature bounds raised for tropical Australia
#   Wetness-period bound 15-25 C becomes 15-32 C; preceding 5-day bound 20-25 C
#   becomes 20-30 C. Lower bounds unchanged. The shortest wetness duration for
#   infection occurs near 25 to 28 C, exactly where the original cap sits.
#   Restore with BLASTAM_TWET_MAX <- 25 and BLASTAM_PREV5_MAX <- 25.
#
# DEVIATION 3. Heavy rain exclusion stated explicitly
#   Hours with precipitation >= BLASTAM_RAIN_HEAVY are excluded from the
#   infection-conducive wet count. Disable with BLASTAM_RAIN_HEAVY <- Inf.
#
# DEVIATION 4. Leaf wetness estimated from hourly RH
#   An hour counts as wet when RH >= BLASTAM_RH_WET or rain >= BLASTAM_RAIN_WET.
#   ERA5 has no leaf-wetness variable and the energy-balance reconstruction is
#   impractical from three variables. Absolute wet-hour counts are provisional
#   until calibrated against field leaf-wetness data.
#
# DEVIATION 5. Local solar time rather than political timezone
#   The 15:00-09:00 window is applied in local solar time (lon / 15 h), because
#   political zone boundaries would put a step discontinuity in the wet-hour
#   count across state lines on the continental map.
#
# ==========================================================================
# THREE FIXES IN SCHEMA VERSION 3
# ==========================================================================
#
# FIX A. The model day is cut at BLASTAM_DAY_CUT_HOUR (10:00 local solar), not
#   at local midnight. Schema 2 correctly moved the BLASTAM night window onto
#   local solar time but moved the DAILY aggregates with it, onto midnight days.
#   EPIRICE's rainlim gate is a daily SUM, so cutting at midnight splits a
#   nocturnal rain event across two days and halves the peak daily total. On a
#   synthetic series with 6 mm falling between 22:00 and 03:00 local every third
#   night, days reaching the 5 mm gate went from 24 to 0 and final intensity from
#   0.0616% to 0.0000%, on identical rainfall. The measured symptom was Malanda
#   falling from 0.374% to 0.006% between the 2026-07-28 and 2026-07-29 runs.
#   A model day labelled 23 July now runs 10:00 on 23 July to 09:59 on 24 July,
#   local solar, so the whole nocturnal wet and rain period sits in one day and
#   the night window is a subset of it by construction.
#
# FIX B. prev5 is now genuinely the PRECEDING five days. It was
#   frollmean(TEMP, 5, align = "right"), which covers days i-4 to i and so
#   includes the night's own day. It is now lagged one day, and it is computed on
#   a COMPLETE date sequence so a missing day yields NA rather than a mean that
#   silently spans six or more calendar days.
#
# FIX C. Night completeness. A night used to be judged on a single evening hour
#   plus a single morning hour, and an hour with missing humidity counted as dry,
#   so a mostly empty night scored "not favourable" rather than "not judged". A
#   night now needs BLASTAM_MIN_EVE_HOURS and BLASTAM_MIN_MORN_HOURS of usable
#   data and no more than BLASTAM_MAX_NA_FRAC unusable hours.
#
# ==========================================================================
# REFERENCES
# ==========================================================================
#   Koshimizu, Y. (1988) Bull. Tohoku Natl. Agric. Exp. Stn. 78: 67-121 [Jpn]
#   Hayashi, T. and Koshimizu, Y. (1988) ibid. 78: 123-138 [Jpn]
#   Barksdale, T.H. and Jones, M.W. (1965) Phytopathology 55: 1037-1040.
#   Kato, H. (1974) Epidemiology of rice blast disease. Rev. Plant Prot. Res.
#     7: 1-20.
#   Kato, H. and Kozaka, T. (1974) Effect of temperature on lesion enlargement
#     and sporulation of Pyricularia oryzae in rice leaves. Phytopathology 64:
#     828-830. doi:10.1094/Phyto-64-828.
#   Maehara, H. and Yamada, M. (2025) Ann. Rep. Soc. Pl. Prot. North Japan 76: 41-46.
################################################################################

suppressPackageStartupMessages(library(data.table))

# ---- Parameters: fallbacks only, blast_config.R is authoritative ------------
.def <- function(nm, value) if (!exists(nm, inherits = TRUE)) assign(nm, value, envir = globalenv())

.def("BLASTAM_TWET_MIN",  15)
.def("BLASTAM_TWET_MAX",  32)
.def("BLASTAM_PREV5_MIN", 20)
.def("BLASTAM_PREV5_MAX", 30)
.def("BLASTAM_RH_WET",     90)
.def("BLASTAM_RAIN_WET",   0.2)
.def("BLASTAM_RAIN_HEAVY", 4.0)
.def("BLASTAM_USE_BJ_THRESHOLD", TRUE)
.def("BLASTAM_WET_HOURS_FIXED",  10L)
.def("BLASTAM_NIGHT_START", 15L)
.def("BLASTAM_NIGHT_END",    9L)
.def("BLASTAM_DAY_CUT_HOUR", 10L)
.def("BLASTAM_LEADIN_DAYS",   6L)
.def("BLASTAM_MIN_EVE_HOURS",  7L)
.def("BLASTAM_MIN_MORN_HOURS", 7L)
.def("BLASTAM_MAX_NA_FRAC",  0.10)
.def("BLASTAM_WINDOW_DAYS", 21L)
.def("BLASTAM_RECENT_DAYS",  7L)

# Barksdale & Jones (1965) minimum infection curve (lower 95% CI): hours of leaf
# wetness required for infection at each temperature. Fahrenheit original
# converted to Celsius. Flat extrapolation outside the measured range.
.BJ_TEMP <- c(15.56, 18.33, 21.11, 23.89, 26.67)   # (60, 65, 70, 75, 80 F)
.BJ_HRS  <- c(12.2,  10.9,   9.7,   8.6,   7.7)
blastam_bj_min_hours <- approxfun(.BJ_TEMP, .BJ_HRS, rule = 2)

# ---- Helpers ---------------------------------------------------------------

blastam_solar_offset_h <- function(lon) lon / 15

# Required wet hours for a given wetness-period temperature. Vectorised, and NOW
# THE ONLY PLACE the rule is written down: the same logic used to be duplicated
# inline inside blastam_daily_from_hourly(), so a change to one drifted from the
# other.
blastam_min_hours <- function(temp_wet) {
  fixed <- as.numeric(BLASTAM_WET_HOURS_FIXED)
  if (!isTRUE(BLASTAM_USE_BJ_THRESHOLD)) return(rep(fixed, length(temp_wet)))
  ifelse(is.na(temp_wet), fixed, blastam_bj_min_hours(temp_wet))
}

# Is the tail + lead-in arithmetic still exactly the API's 14 day minimum charge?
blastam_check_fetch_arithmetic <- function() {
  if (!exists("REFRESH_TAIL_DAYS") || !exists("DAY_CUT_LAG_DAYS")) return(invisible(NA))
  total <- as.integer(REFRESH_TAIL_DAYS) + as.integer(BLASTAM_LEADIN_DAYS) +
           as.integer(DAY_CUT_LAG_DAYS)
  if (total != 14L)
    warning(sprintf(paste0("REFRESH_TAIL_DAYS + BLASTAM_LEADIN_DAYS + DAY_CUT_LAG_DAYS ",
                           "= %d, not 14. A refresh will cost %.2f weighted calls per ",
                           "point instead of 1.00."), total, max(1, total / 14)),
            call. = FALSE)
  invisible(total)
}

################################################################################
# blastam_daily_from_hourly()
#
# Convert one point's hourly UTC series into per-day rows for BOTH models:
#   EPIRICE inputs: TEMP (daily mean C), RHUM (daily mean %), RAIN (daily mm)
#   BLASTAM inputs: wet_hours, temp_wet, infect (0/1/NA), semi (0/1/NA)
#
#   hourly  data.table(dt = POSIXct UTC, temp, rh, rain)
#   lon     longitude degrees east. Pass this so the model day and the night
#           window are in local solar time. If NULL a warning is issued and the
#           timestamps are treated as local.
#
# `date` labels a model day that STARTS at BLASTAM_DAY_CUT_HOUR local solar and
# runs 24 hours. The night of that day (15:00 to 09:00) is contained within it.
################################################################################
blastam_daily_from_hourly <- function(hourly, lon = NULL) {
  if (is.null(hourly) || nrow(hourly) == 0L) return(NULL)
  h <- data.table::copy(as.data.table(hourly))
  setorder(h, dt)

  if (is.null(lon)) {
    warning("blastam_daily_from_hourly(): no lon supplied; treating timestamps ",
            "as local. Pass lon to apply the model day in local solar time.",
            call. = FALSE)
    off <- 0
  } else {
    off <- blastam_solar_offset_h(lon)
  }
  # Local solar clock. dt carries tzone "UTC" from the fetchers; as.Date() and
  # format() both then read UTC, so the two agree. Do not remove the tz on dt.
  h[, sdt   := dt + off * 3600]
  h[, shour := as.integer(format(sdt, "%H"))]
  # FIX A: the model day starts at BLASTAM_DAY_CUT_HOUR, so subtract the cut
  # before taking the date. Hours before the cut belong to the previous day.
  h[, mday := as.Date(sdt - as.numeric(BLASTAM_DAY_CUT_HOUR) * 3600)]

  # An hour is USABLE only when all three variables are present. Previously a
  # missing humidity value was coerced to "not wet", which reported dry weather
  # where the truth was no data.
  h[, usable := !is.na(temp) & !is.na(rh) & !is.na(rain)]

  # Daily aggregates over COMPLETE model days only: a partial day biases the
  # daily mean and truncates the rain sum, which is what the EPIRICE gate reads.
  hrs <- h[, .(nh = .N, n_bad = sum(!usable)), by = mday]
  full_days <- hrs[nh >= 24L & n_bad <= BLASTAM_MAX_NA_FRAC * nh, mday]
  if (length(full_days) == 0L) return(NULL)
  hf <- h[mday %in% full_days & usable]

  daily <- hf[, .(TEMP = mean(temp), RHUM = mean(rh), RAIN = sum(rain)),
              by = .(date = mday)]
  setorder(daily, date)

  # FIX B: the PRECEDING five days, lagged one day, computed on a complete date
  # sequence so a gap gives NA instead of a mean spanning more than five days.
  span <- data.table(date = seq(min(daily$date), max(daily$date), by = "day"))
  span <- merge(span, daily[, .(date, TEMP)], by = "date", all.x = TRUE)
  setorder(span, date)
  # FULLY QUALIFIED DELIBERATELY. run_blast_grid.R loads terra after data.table,
  # and terra::shift is an S4 generic that masks data.table::shift, so the bare
  # call failed with "unable to find an inherited method for function 'shift'".
  # tryCatch in the caller turned that into every grid point silently coming back
  # "empty", i.e. an entirely blank map with no error in the log.
  span[, prev5 := data.table::shift(
    data.table::frollmean(TEMP, 5L, align = "right"), 1L)]
  daily <- merge(daily, span[, .(date, prev5)], by = "date", all.x = TRUE)

  # An hour is infection-conducive wet when:
  #   (a) it is usable, AND
  #   (b) RH >= BLASTAM_RH_WET or rain >= BLASTAM_RAIN_WET, AND
  #   (c) rain < BLASTAM_RAIN_HEAVY   [Yoshino / Deviation 3]
  h[, infection_wet := usable &
      (rh >= BLASTAM_RH_WET | rain >= BLASTAM_RAIN_WET) &
      (rain < BLASTAM_RAIN_HEAVY)]

  # The night window is a subset of the model day, so night == mday wherever the
  # hour falls inside 15:00-09:00. No separate day arithmetic is needed.
  h[, in_night := shour >= BLASTAM_NIGHT_START | shour < BLASTAM_NIGHT_END]

  nh <- h[in_night == TRUE, .(
    wet_hours = sum(infection_wet),
    temp_wet  = if (any(infection_wet)) mean(temp[infection_wet], na.rm = TRUE)
                else NA_real_,
    n_eve   = sum(shour >= BLASTAM_NIGHT_START & usable),
    n_morn  = sum(shour <  BLASTAM_NIGHT_END   & usable),
    n_hours = .N,
    n_bad   = sum(!usable)
  ), by = .(date = mday)]

  # FIX C: a night must be substantially observed to be judged at all.
  nh[, complete := n_eve >= BLASTAM_MIN_EVE_HOURS &
                   n_morn >= BLASTAM_MIN_MORN_HOURS &
                   n_bad <= BLASTAM_MAX_NA_FRAC * n_hours]

  out <- merge(daily, nh[, .(date, wet_hours, temp_wet, complete)],
               by = "date", all.x = TRUE)
  out[is.na(wet_hours), wet_hours := 0]
  out[is.na(complete),  complete  := FALSE]

  ok <- out$complete & !is.na(out$prev5)
  out[, `:=`(infect = NA_integer_, semi = NA_integer_)]

  min_h <- blastam_min_hours(out$temp_wet)

  out[ok & wet_hours < min_h, `:=`(infect = 0L, semi = 0L)]

  jud <- ok & out$wet_hours >= min_h
  if (any(jud)) {
    c_twet <- !is.na(out$temp_wet) &
              out$temp_wet >= BLASTAM_TWET_MIN &
              out$temp_wet <= BLASTAM_TWET_MAX
    c_prev <- out$prev5 >= BLASTAM_PREV5_MIN &
              out$prev5 <= BLASTAM_PREV5_MAX
    n_met  <- as.integer(c_twet) + as.integer(c_prev)
    out[jud, `:=`(infect = as.integer(n_met[jud] == 2L),
                  semi   = as.integer(n_met[jud] == 1L))]
  }

  setorder(out, date)
  out[, .(date, TEMP, RHUM, RAIN, wet_hours, temp_wet, infect, semi)]
}

################################################################################
# blastam_score()
#
# Count favourable infection days in the rolling window ending at end_date.
#
# THE WINDOW IS NOW BOUNDED AT BOTH ENDS. It was `dates > (end_date - window)`
# with no upper bound, so any row later than end_date was also counted. Because
# run_blast.R passed (end_date - BLASTAM_END_LAG_DAYS) while still holding rows to
# the archive edge, the town table reported a 22 day count and an 8 day "7d"
# count, and disagreed with the map, which used the correct
# `> (bend - win) & <= bend` form. Both now go through this function.
#
# `unjudged` tells you how many days in the window could not be assessed, so a
# low count from missing data is distinguishable from one from dry weather.
################################################################################
blastam_score <- function(infect, semi, dates, end_date,
                          window = BLASTAM_WINDOW_DAYS,
                          recent = BLASTAM_RECENT_DAYS) {
  dates    <- as.Date(dates)
  end_date <- as.Date(end_date)
  window   <- as.integer(window)
  recent   <- as.integer(recent)
  inwin <- dates > (end_date - window) & dates <= end_date
  inrec <- dates > (end_date - recent) & dates <= end_date
  list(events   = sum(infect[inwin], na.rm = TRUE),
       combined = sum(infect[inwin], na.rm = TRUE) +
                  sum(semi[inwin],   na.rm = TRUE),
       recent   = sum(infect[inrec], na.rm = TRUE),
       unjudged = sum(is.na(infect[inwin])),
       n_days   = sum(inwin))
}

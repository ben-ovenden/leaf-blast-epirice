################################################################################
# blastam_model.R  -  BLASTAM leaf blast infection-warning model
#
# WHAT IT IS
# BLASTAM (Koshimizu 1988) is the infection-warning model operated across Japan
# on AMeDAS weather data. For each day it judges whether that day's leaf-wetness
# period and temperatures were favourable for a Magnaporthe oryzae (leaf blast)
# infection event. It does NOT simulate the epidemic; it counts favourable
# infection days. This is fundamentally different from EPIRICE, which simulates
# disease build-up (intensity).
#
# INFECTION CRITERIA. A day is FAVOURABLE when ALL THREE hold:
#   1. leaf wetness duration >= 10 hours (including the night);
#   2. mean air temperature DURING the wetness period is within bounds;
#   3. mean air temperature of the preceding 5 days is within bounds.
# If wetness is >= 10 h but exactly one temperature condition fails, the day is
# SEMI-FAVOURABLE. Otherwise it is not favourable.
#
# Original Japanese BLASTAM bounds (Koshimizu 1988): wetness-period temperature
# 15-25 C; preceding 5-day mean 20-25 C.
#
# LOCAL ADAPTATION (northern Australia) -- USED HERE. The upper temperature
# bounds are raised to: wetness-period 15-32 C; preceding 5-day mean 20-30 C.
# Rationale: rice blast infects efficiently from about 25 up to ~32 C when leaf
# wetness is adequate (the required wetness duration is in fact shortest near
# 25-28 C), so the strict Japanese 25 C caps wrongly exclude warm, humid
# wet-season nights in the tropics. This is a DELIBERATE deviation from Koshimizu
# (1988); the lower bounds are unchanged. To reproduce the original Japanese
# model, set BLASTAM_TWET_MAX <- 25 and BLASTAM_PREV5_MAX <- 25 below.
#
# TIME BASE  ** this was wrong and is the most consequential fix in this file **
# The wetness window runs 15:00 to 09:00, which only means "evening to morning"
# if the clock is LOCAL. The weather is fetched in UTC, and the previous version
# applied the 15:00-09:00 window directly to UTC timestamps. At longitude 145
# (north Queensland) that window is 01:00 to 19:00 local: it misses the entire
# evening dew onset and includes the middle of the day. A test series built with
# 13 wet hours every local night returned 8 wet hours and zero favourable days.
#
# Hours are therefore converted to LOCAL SOLAR time from longitude
# (offset = lon / 15 hours) before the day and night windows are applied. Solar
# time rather than the political timezone, because political zones jump at state
# borders and would put a discontinuity through the middle of a continental map,
# and because dew formation follows the sun, not the clock.
#
# COMPLETENESS
# A night needs hours on both sides of midnight, and the preceding 5-day mean
# needs five earlier days. Where either is unavailable, infect and semi are
# returned as NA rather than 0. Returning 0 asserts "conditions were not
# favourable" when the truth is "not yet known", which biased the most recent
# day of every fetch to zero and, once short refresh windows were introduced,
# would write spurious zeros over good cached values. Callers should fetch
# BLASTAM_LEADIN_DAYS of lead-in and discard it.
#
# LEAF WETNESS
# The original BLASTAM estimates leaf wetness from AMeDAS temperature, rainfall,
# sunshine duration and wind (AMeDAS has no humidity sensor). Here we have hourly
# relative humidity directly from ERA5, so leaf wetness is taken as RH >= 90% or
# rain, the standard humidity-based leaf-wetness proxy. Treat the absolute counts
# as provisional and calibrate against local observation.
#
# REFERENCES
#  Koshimizu, Y. (1988) A forecasting method for occurrence of rice leaf blast
#    with AMeDAS data. Bulletin of the Tohoku National Agricultural Experiment
#    Station 78: 67-121. [in Japanese]
#  Hayashi, T. and Koshimizu, Y. (1988) Computer program BLASTAM for forecasting
#    occurrence of rice leaf blast. Bull. Tohoku Natl. Agric. Exp. Stn. 78:
#    123-138. [in Japanese]
#  Maehara, H. and Yamada, M. (2025) Annual changes in the timing and frequency
#    of favorable conditions for rice leaf blast infection estimated by BLASTAM
#    in Fukushima Prefecture. Ann. Rep. Soc. Pl. Prot. North Japan 76: 41-46.
################################################################################

suppressPackageStartupMessages(library(data.table))

# ---- Parameters -----------------------------------------------------------
BLASTAM_WET_HOURS_MIN <- 10    # leaf wetness duration threshold (hours)
BLASTAM_TWET_MIN      <- 15    # mean temp during wetness, lower bound (C)
BLASTAM_TWET_MAX      <- 32    # mean temp during wetness, upper bound (C)  [Japan: 25]
BLASTAM_PREV5_MIN     <- 20    # preceding 5-day mean temp, lower bound (C)
BLASTAM_PREV5_MAX     <- 30    # preceding 5-day mean temp, upper bound (C) [Japan: 25]
BLASTAM_RH_WET        <- 90    # RH (%) at or above = leaf wet
BLASTAM_RAIN_WET      <- 0.2   # precipitation (mm) at or above = leaf wet

# Night window, in LOCAL SOLAR hours.
BLASTAM_NIGHT_START   <- 15L   # hours from this one onward belong to that night
BLASTAM_NIGHT_END     <- 9L    # hours before this one belong to the previous night

# Lead-in a caller must fetch and then discard: 1 day for the solar shift (the
# first local day is partial) plus 5 for the preceding 5-day mean.
BLASTAM_LEADIN_DAYS   <- 6L

# ---- Helpers ---------------------------------------------------------------

# Solar local offset in hours from longitude. Smooth across the continent, so no
# discontinuity at state borders.
blastam_solar_offset_h <- function(lon) lon / 15

# Judge one night from its hourly vectors. Retained for direct use and tests;
# blastam_daily_from_hourly() does the same test vectorised.
blastam_night <- function(temp, rh, rain, prev5_mean, complete = TRUE) {
  wet <- (rh >= BLASTAM_RH_WET) | (rain >= BLASTAM_RAIN_WET)
  wet[is.na(wet)] <- FALSE
  wh <- sum(wet)
  tw <- if (wh > 0) mean(temp[wet], na.rm = TRUE) else NA_real_

  if (!isTRUE(complete) || is.na(prev5_mean))
    return(data.table(infect = NA_integer_, semi = NA_integer_,
                      wet_hours = wh, temp_wet = tw))
  if (wh < BLASTAM_WET_HOURS_MIN)
    return(data.table(infect = 0L, semi = 0L, wet_hours = wh, temp_wet = tw))

  c_twet  <- !is.na(tw) && tw >= BLASTAM_TWET_MIN && tw <= BLASTAM_TWET_MAX
  c_prev5 <- prev5_mean >= BLASTAM_PREV5_MIN && prev5_mean <= BLASTAM_PREV5_MAX
  n_met <- c_twet + c_prev5
  data.table(infect = if (n_met == 2L) 1L else 0L,
             semi   = if (n_met == 1L) 1L else 0L,
             wet_hours = wh, temp_wet = tw)
}

################################################################################
# Process one point's hourly series into per-day rows carrying BOTH models'
# inputs: EPIRICE daily aggregates (TEMP, RHUM, RAIN) and the BLASTAM night
# judgement (wet_hours, temp_wet, infect, semi).
#
#   hourly  data.table(dt = POSIXct UTC, temp, rh, rain)
#   lon     longitude in degrees east, for the solar time shift. If NULL, dt is
#           assumed local and a warning is issued, because silently treating UTC
#           as local is exactly the bug this argument exists to prevent.
#
# Night d spans BLASTAM_NIGHT_START on day d to BLASTAM_NIGHT_END on day d+1, in
# local solar time. Partial local days at either end are dropped.
################################################################################
blastam_daily_from_hourly <- function(hourly, lon = NULL) {
  if (is.null(hourly) || nrow(hourly) == 0L) return(NULL)
  h <- data.table::copy(as.data.table(hourly))
  setorder(h, dt)

  if (is.null(lon)) {
    warning("blastam_daily_from_hourly(): no lon given, treating timestamps as ",
            "local. Pass lon so the night window is local solar time.",
            call. = FALSE)
    off <- 0
  } else {
    off <- blastam_solar_offset_h(lon)
  }
  h[, sdt := dt + off * 3600]
  h[, `:=`(sday = as.Date(sdt), shour = as.integer(format(sdt, "%H")))]

  # Daily aggregates over COMPLETE local days only. A partial day at either end
  # biases the means and truncates the rainfall sum.
  hrs <- h[, .(nh = .N), by = sday]
  full_days <- hrs[nh >= 24L, sday]
  if (length(full_days) == 0L) return(NULL)
  hf <- h[sday %in% full_days]

  daily <- hf[, .(TEMP = mean(temp, na.rm = TRUE),
                  RHUM = mean(rh,   na.rm = TRUE),
                  RAIN = sum(rain,  na.rm = TRUE)), by = .(date = sday)]
  setorder(daily, date)
  # Trailing 5-day mean of daily mean temperature. NA for the first four days,
  # which is why callers fetch BLASTAM_LEADIN_DAYS of lead-in.
  daily[, prev5 := frollmean(TEMP, 5, align = "right")]

  # Vectorised night assignment. The previous version scanned the whole hourly
  # series once per day, which is O(days * hours) per point.
  h[, wet := ((rh >= BLASTAM_RH_WET) | (rain >= BLASTAM_RAIN_WET)) %in% TRUE]
  h[, night := fifelse(shour >= BLASTAM_NIGHT_START, sday,
                fifelse(shour <  BLASTAM_NIGHT_END,  sday - 1L, as.Date(NA)))]
  nh <- h[!is.na(night), .(
    wet_hours = sum(wet),
    temp_wet  = if (any(wet)) mean(temp[wet], na.rm = TRUE) else NA_real_,
    n_eve     = sum(shour >= BLASTAM_NIGHT_START),
    n_morn    = sum(shour <  BLASTAM_NIGHT_END)
  ), by = night]
  # A night is complete only with hours on both sides of midnight.
  nh[, complete := n_eve > 0L & n_morn > 0L]

  out <- merge(daily, nh[, .(date = night, wet_hours, temp_wet, complete)],
               by = "date", all.x = TRUE)
  out[is.na(wet_hours), wet_hours := 0]
  out[is.na(complete), complete := FALSE]

  ok  <- out$complete & !is.na(out$prev5)
  out[, `:=`(infect = NA_integer_, semi = NA_integer_)]
  out[ok & wet_hours <  BLASTAM_WET_HOURS_MIN, `:=`(infect = 0L, semi = 0L)]

  jud <- ok & out$wet_hours >= BLASTAM_WET_HOURS_MIN
  if (any(jud)) {
    c_twet <- !is.na(out$temp_wet) & out$temp_wet >= BLASTAM_TWET_MIN &
              out$temp_wet <= BLASTAM_TWET_MAX
    c_prev <- out$prev5 >= BLASTAM_PREV5_MIN & out$prev5 <= BLASTAM_PREV5_MAX
    n_met  <- as.integer(c_twet) + as.integer(c_prev)
    out[jud, `:=`(infect = as.integer(n_met[jud] == 2L),
                  semi   = as.integer(n_met[jud] == 1L))]
  }

  setorder(out, date)
  out[, .(date, TEMP, RHUM, RAIN, wet_hours, temp_wet, infect, semi)]
}

# BLASTAM score over the most recent `window` days (default 21): the count of
# FAVOURABLE infection days (and favourable plus semi-favourable), indicating
# where infection pressure is building NOW rather than a whole-season total.
# `recent` is a shorter sub-count for the trend. Days that could not be judged
# are excluded from the counts and reported in `unjudged`, so a low count caused
# by missing data is distinguishable from one caused by unfavourable weather.
blastam_score <- function(infect, semi, dates, end_date, window = 21L, recent = 7L) {
  inwin <- dates > (end_date - window)
  list(events   = sum(infect[inwin], na.rm = TRUE),
       combined = sum(infect[inwin], na.rm = TRUE) + sum(semi[inwin], na.rm = TRUE),
       recent   = sum(infect[dates > (end_date - recent)], na.rm = TRUE),
       unjudged = sum(is.na(infect[inwin])))
}

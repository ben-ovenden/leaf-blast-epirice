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
# Original Japanese BLASTAM bounds (Koshimizu 1988, as operated by the Fukushima,
# Kochi, Saga, Ehime, Miyagi, Kyoto and Hokkaido plant-protection stations):
# wetness-period temperature 15-25 C; preceding 5-day mean 20-25 C.
#
# LOCAL ADAPTATION (northern Australia) -- USED HERE. The upper temperature
# bounds are raised to: wetness-period 15-32 C; preceding 5-day mean 20-30 C.
# Rationale: rice blast infects efficiently from about 25 up to ~32 C when leaf
# wetness is adequate (the required wetness duration is in fact shortest near
# 25-28 C), so the strict Japanese 25 C caps wrongly exclude warm, humid
# wet-season nights in the tropics. This is a DELIBERATE deviation from Koshimizu
# (1988), calibrated for warmer conditions; the lower bounds are unchanged. To
# reproduce the original Japanese model, set BLASTAM_TWET_MAX <- 25 and
# BLASTAM_PREV5_MAX <- 25 below.
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

# ---- Parameters -----------------------------------------------------------
# Upper temperature bounds are RAISED from the original Japanese BLASTAM to suit
# warmer northern-Australian conditions (see LOCAL ADAPTATION in the header). For
# the strict Japanese model, set BLASTAM_TWET_MAX <- 25 and BLASTAM_PREV5_MAX <- 25.
BLASTAM_WET_HOURS_MIN <- 10    # leaf wetness duration threshold (hours)
BLASTAM_TWET_MIN      <- 15    # mean temp during wetness, lower bound (C)
BLASTAM_TWET_MAX      <- 32    # mean temp during wetness, upper bound (C)  [Japan: 25]
BLASTAM_PREV5_MIN     <- 20    # preceding 5-day mean temp, lower bound (C)
BLASTAM_PREV5_MAX     <- 30    # preceding 5-day mean temp, upper bound (C) [Japan: 25]
BLASTAM_RH_WET        <- 90    # RH (%) at or above = leaf wet
BLASTAM_RAIN_WET      <- 0.2   # precipitation (mm) at or above = leaf wet

# Judge one day's infection from hourly vectors over the overnight wetness
# window, plus the preceding 5-day mean air temperature. Returns a one-row
# data.table: infect (1 = favourable), semi (1 = semi-favourable), wet_hours,
# temp_wet.
blastam_night <- function(temp, rh, rain, prev5_mean) {
  none <- data.table::data.table(infect = 0L, semi = 0L, wet_hours = 0,
                                 temp_wet = NA_real_)
  wet <- (rh >= BLASTAM_RH_WET) | (rain >= BLASTAM_RAIN_WET)
  wet[is.na(wet)] <- FALSE
  wh <- sum(wet)
  if (wh < BLASTAM_WET_HOURS_MIN) { none$wet_hours <- wh; return(none) }
  tw <- mean(temp[wet], na.rm = TRUE)
  c_twet  <- !is.na(tw) && tw >= BLASTAM_TWET_MIN && tw <= BLASTAM_TWET_MAX
  c_prev5 <- !is.na(prev5_mean) && prev5_mean >= BLASTAM_PREV5_MIN &&
             prev5_mean <= BLASTAM_PREV5_MAX
  n_met <- c_twet + c_prev5
  inf <- if (n_met == 2L) 1L else 0L
  sem <- if (n_met == 1L) 1L else 0L
  data.table::data.table(infect = inf, semi = sem, wet_hours = wh, temp_wet = tw)
}

# Process one point's hourly series into per-day rows carrying BOTH models'
# inputs: EPIRICE daily aggregates (TEMP, RHUM, RAIN) and the BLASTAM day
# judgement (wet_hours, temp_wet, infect, semi). The BLASTAM day d spans 15:00 on
# d to 09:00 on d+1 (the evening-to-morning wetness period).
#
# hourly: data.table(dt = POSIXct UTC, temp, rh, rain)
blastam_daily_from_hourly <- function(hourly) {
  h <- data.table::copy(hourly)
  data.table::setorder(h, dt)
  h[, date := as.Date(dt)]
  daily <- h[, .(TEMP = mean(temp, na.rm = TRUE),
                 RHUM = mean(rh,  na.rm = TRUE),
                 RAIN = sum(rain, na.rm = TRUE)), by = date]
  data.table::setorder(daily, date)
  # preceding 5-day mean air temperature (trailing 5-day mean of daily mean temp)
  daily[, prev5 := data.table::frollmean(TEMP, 5, align = "right")]
  hh <- as.integer(format(h$dt, "%H"))
  res <- vector("list", nrow(daily))
  for (i in seq_len(nrow(daily))) {
    d <- daily$date[i]
    sel <- (h$date == d & hh >= 15) | (h$date == d + 1 & hh < 9)
    res[[i]] <- if (!any(sel)) blastam_night(numeric(0), numeric(0), numeric(0), NA)
                else blastam_night(h$temp[sel], h$rh[sel], h$rain[sel], daily$prev5[i])
  }
  nb <- data.table::rbindlist(res)
  cbind(daily[, .(date, TEMP, RHUM, RAIN)],
        nb[, .(wet_hours, temp_wet, infect, semi)])
}

# BLASTAM score over the most recent `window` days (default 21): the count of
# FAVOURABLE infection days (and favourable + semi-favourable), indicating where
# infection pressure is building NOW rather than a whole-season total. `recent`
# is a shorter sub-count for the trend.
blastam_score <- function(infect, semi, dates, end_date, window = 21L, recent = 7L) {
  inwin <- dates > (end_date - window)
  fav  <- sum(infect[inwin], na.rm = TRUE)
  comb <- fav + sum(semi[inwin], na.rm = TRUE)
  rec  <- sum(infect[dates > (end_date - recent)], na.rm = TRUE)
  list(events = fav, combined = comb, recent = rec)
}

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
# ORIGINAL WETNESS THRESHOLD
# Koshimizu (1988) used a fixed threshold of 10 hours, derived empirically from
# AMeDAS data in temperate Tohoku. The underlying biology shows the required
# duration decreases with temperature (Barksdale & Jones 1965; Kato 1974):
# roughly 12.2 h at 15.6 C and 7.7 h at 26.7 C. The fixed 10 h is a reasonable
# approximation for Japanese temperate conditions, where wetness-period
# temperatures are typically 15-20 C; it over-counts cool nights in Australia.
#
# ORIGINAL LEAF WETNESS ESTIMATION (Koshimizu 1988)
# AMeDAS has no humidity sensor. Koshimizu estimated leaf wetness from sunshine
# duration, precipitation and wind speed using a leaf energy balance. Hours with
# precipitation >= 4 mm/h were excluded from the infection-conducive count
# (following Yoshino 1988) because heavy rain washes conidia off the leaf.
#
# NIGHT WINDOW (Koshimizu 1988)
# 15:00 to 09:00 the following day, LOCAL time (Japan Standard Time for AMeDAS).
# This code uses local solar time (longitude / 15 h) instead of political
# timezone, to avoid a step discontinuity across the continental map.
#
# ==========================================================================
# DEVIATIONS USED HERE, AND RATIONALE
# ==========================================================================
#
# DEVIATION 1. Temperature-dependent wetness threshold (Deviation from fixed 10 h)
#   WHAT: the required wetness duration is computed from the Barksdale & Jones
#     (1965) minimum infection curve (lower 95% CI): 12.2 h at 15.6 C,
#     decreasing to 7.7 h at 26.7 C.
#   WHY: the fixed 10 h threshold was calibrated for temperate Japan (wetness-
#     period temps typically 15-20 C, where the B&J curve also gives ~10-12 h).
#     In tropical Australia, wetness-period temps are typically 22-28 C, where
#     the B&J curve gives ~8-9 h. Using 10 h in the tropics over-counts
#     infection-conducive nights by ~41% relative to the B&J threshold (measured
#     on the 2026-07-29 cache). In the cool dry season it under-counts, since
#     17 C nights need 11.5 h but are scored against the 10 h bar.
#   HOW TO RESTORE: set BLASTAM_USE_BJ_THRESHOLD <- FALSE below.
#   Ref: Barksdale, T.H. and Jones, M.W. (1965) Phytopathology 55: 1037-1040.
#        Kato, H. (1974) Rev. Plant Prot. Res. 7: 1-20.
#
# DEVIATION 2. Upper temperature bounds raised for tropical Australia
#   WHAT: wetness-period bound raised from 15-25 C to 15-32 C; preceding 5-day
#     bound raised from 20-25 C to 20-30 C. Lower bounds unchanged.
#   WHY: blast infects efficiently from ~15 C to ~32 C when leaf wetness is
#     adequate. The shortest wetness duration for infection occurs near 25-28 C
#     (Barksdale & Jones 1965; Kato 1974), exactly the temperatures excluded by
#     the original 25 C cap. Koshimizu's bounds reflect Tohoku conditions; they
#     wrongly classify warm tropical wet-season nights as non-favourable.
#   HOW TO RESTORE: set BLASTAM_TWET_MAX <- 25 and BLASTAM_PREV5_MAX <- 25.
#
# DEVIATION 3. Heavy rain exclusion rate matched to Yoshino (1988)
#   WHAT: hours with precipitation >= BLASTAM_RAIN_HEAVY (4 mm/h) are excluded
#     from the infection-conducive wet count.
#   WHY: Yoshino (1988, cited in Hayashi & Koshimizu 1988) excluded hours with
#     rain >= 4 mm/h because heavy rain washes conidia off leaves. The original
#     BLASTAM handled this implicitly through the energy-balance estimator;
#     our RH proxy requires it to be stated explicitly.
#   HOW TO DISABLE: set BLASTAM_RAIN_HEAVY <- Inf.
#
# DEVIATION 4. Leaf wetness estimated from hourly RH (not energy balance)
#   WHAT: a leaf hour is counted as wet when RH >= 90% or rain >= 0.2 mm/h.
#   WHY: ERA5 provides no measured leaf-wetness variable; the energy-balance
#     reconstruction (sunshine, wind, rain) is impractical from the three
#     variables fetched. RH >= 90% is the standard humidity proxy for leaf
#     wetness (Kim et al. 2002; Dalla Marta et al. 2005). The 0.2 mm/h rain
#     threshold picks up drizzle that the RH sensor may lag.
#   NOTE: the mean temp during the wet period is then the mean during RH >= 90%
#     or light rain hours (excluding heavy rain); this differs from the AMeDAS
#     energy-balance approach and absolute counts should be treated as
#     provisional until calibrated against field leaf-wetness data.
#
# DEVIATION 5. Local solar time rather than political timezone
#   WHAT: the 15:00-09:00 window is applied in local solar time (lon / 15 h).
#   WHY: political timezone boundaries would put a discontinuity in the wet-hour
#     count across state lines on the continental map. Solar time follows the
#     sun, which drives dew formation and leaf drying.
#
# ==========================================================================
# REFERENCES
# ==========================================================================
#   Koshimizu, Y. (1988) A forecasting method for occurrence of rice leaf blast
#     with AMeDAS data. Bull. Tohoku Natl. Agric. Exp. Stn. 78: 67-121 [Japanese]
#   Hayashi, T. and Koshimizu, Y. (1988) Computer program BLASTAM for forecasting
#     occurrence of rice leaf blast. Bull. Tohoku Natl. Agric. Exp. Stn. 78:
#     123-138 [Japanese]
#   Barksdale, T.H. and Jones, M.W. (1965) Minimum conditions of temperature and
#     dew period for infection of rice by Piricularia oryzae. Phytopathology 55:
#     1037-1040.
#   Kato, H. (1974) Epidemiology of blast. Rev. Plant Prot. Res. 7: 1-20.
#   Maehara, H. and Yamada, M. (2025) Annual changes in the timing and frequency
#     of favorable conditions for rice leaf blast infection estimated by BLASTAM
#     in Fukushima Prefecture. Ann. Rep. Soc. Pl. Prot. North Japan 76: 41-46.
################################################################################

suppressPackageStartupMessages(library(data.table))

# ---- Parameters -----------------------------------------------------------
# Temperature bounds (Koshimizu 1988 originals in brackets)
BLASTAM_TWET_MIN  <- 15    # wetness-period mean lower bound (C)  [Japan: 15]
BLASTAM_TWET_MAX  <- 32    # wetness-period mean upper bound (C)  [Japan: 25]
BLASTAM_PREV5_MIN <- 20    # preceding 5-day mean lower bound (C) [Japan: 20]
BLASTAM_PREV5_MAX <- 30    # preceding 5-day mean upper bound (C) [Japan: 25]

# Leaf wetness proxy (Deviation 4)
BLASTAM_RH_WET       <- 90    # RH (%) >= this = leaf wet
BLASTAM_RAIN_WET     <- 0.2   # hourly rain (mm) >= this = leaf wet
BLASTAM_RAIN_HEAVY   <- 4.0   # hourly rain (mm) >= this = spores washed off;
                               # those hours excluded from infection-wet count
                               # (Yoshino 1988 via Hayashi & Koshimizu 1988).
                               # Set to Inf to disable.

# Wetness duration threshold (Deviation 1)
# TRUE:  Barksdale & Jones (1965) temperature-dependent minimum
# FALSE: original Koshimizu (1988) fixed 10 h
BLASTAM_USE_BJ_THRESHOLD <- TRUE
BLASTAM_WET_HOURS_FIXED  <- 10L   # used when BLASTAM_USE_BJ_THRESHOLD is FALSE

# Barksdale & Jones (1965) minimum infection curve (lower 95% CI):
# hours of leaf wetness required for infection at each temperature.
# Fahrenheit original converted to Celsius. Flat extrapolation outside range.
.BJ_TEMP <- c(15.56, 18.33, 21.11, 23.89, 26.67)   # (60, 65, 70, 75, 80 F)
.BJ_HRS  <- c(12.2,  10.9,   9.7,   8.6,   7.7)
blastam_bj_min_hours <- approxfun(.BJ_TEMP, .BJ_HRS, rule = 2)

# Night window in LOCAL SOLAR hours (Deviation 5)
BLASTAM_NIGHT_START <- 15L
BLASTAM_NIGHT_END   <- 9L

# Lead-in days callers must fetch and discard: 1 day for the solar shift plus
# 5 days for the preceding 5-day mean.
BLASTAM_LEADIN_DAYS <- 6L

# ---- Helpers ---------------------------------------------------------------

blastam_solar_offset_h <- function(lon) lon / 15

# Required wet hours for a given wetness-period temperature.
blastam_min_hours <- function(temp_wet) {
  if (isTRUE(BLASTAM_USE_BJ_THRESHOLD))
    blastam_bj_min_hours(temp_wet)
  else
    as.numeric(BLASTAM_WET_HOURS_FIXED)
}

################################################################################
# blastam_daily_from_hourly()
#
# Convert one point's hourly UTC series into per-day rows for BOTH models:
#   EPIRICE inputs: TEMP (daily mean C), RHUM (daily mean %), RAIN (daily mm)
#   BLASTAM inputs: wet_hours, temp_wet, infect (0/1/NA), semi (0/1/NA)
#
#   hourly  data.table(dt = POSIXct UTC, temp, rh, rain)
#   lon     longitude degrees east. Pass this so the night window is in local
#           solar time; if NULL a warning is issued and UTC is used as-is.
################################################################################
blastam_daily_from_hourly <- function(hourly, lon = NULL) {
  if (is.null(hourly) || nrow(hourly) == 0L) return(NULL)
  h <- data.table::copy(as.data.table(hourly))
  setorder(h, dt)

  if (is.null(lon)) {
    warning("blastam_daily_from_hourly(): no lon supplied; treating timestamps ",
            "as local. Pass lon to apply the night window in local solar time.",
            call. = FALSE)
    off <- 0
  } else {
    off <- blastam_solar_offset_h(lon)
  }
  h[, sdt   := dt + off * 3600]
  h[, sday  := as.Date(sdt)]
  h[, shour := as.integer(format(sdt, "%H"))]

  # Daily aggregates over COMPLETE local days only: a partial day biases the
  # daily mean and truncates the rain sum.
  hrs       <- h[, .(nh = .N), by = sday]
  full_days <- hrs[nh >= 24L, sday]
  if (length(full_days) == 0L) return(NULL)
  hf <- h[sday %in% full_days]

  daily <- hf[, .(TEMP = mean(temp, na.rm = TRUE),
                  RHUM = mean(rh,   na.rm = TRUE),
                  RAIN = sum(rain,  na.rm = TRUE)), by = .(date = sday)]
  setorder(daily, date)
  daily[, prev5 := frollmean(TEMP, 5, align = "right")]

  # Vectorised night assignment.
  # An hour is infection-conducive wet when:
  #   (a) RH >= BLASTAM_RH_WET or rain >= BLASTAM_RAIN_WET, AND
  #   (b) rain < BLASTAM_RAIN_HEAVY  [Yoshino / Deviation 3]
  h[, infection_wet := (
    (rh >= BLASTAM_RH_WET | rain >= BLASTAM_RAIN_WET) &
    (rain < BLASTAM_RAIN_HEAVY)
  ) %in% TRUE]

  h[, night := fifelse(shour >= BLASTAM_NIGHT_START, sday,
                fifelse(shour <  BLASTAM_NIGHT_END,  sday - 1L, as.Date(NA)))]

  nh <- h[!is.na(night), .(
    wet_hours = sum(infection_wet),
    temp_wet  = if (any(infection_wet)) mean(temp[infection_wet], na.rm = TRUE)
                else NA_real_,
    n_eve  = sum(shour >= BLASTAM_NIGHT_START),
    n_morn = sum(shour <  BLASTAM_NIGHT_END)
  ), by = night]
  nh[, complete := n_eve > 0L & n_morn > 0L]

  out <- merge(daily, nh[, .(date = night, wet_hours, temp_wet, complete)],
               by = "date", all.x = TRUE)
  out[is.na(wet_hours), wet_hours := 0]
  out[is.na(complete),  complete  := FALSE]

  ok <- out$complete & !is.na(out$prev5)
  out[, `:=`(infect = NA_integer_, semi = NA_integer_)]

  # Minimum required hours: temperature-dependent (B&J) or fixed.
  min_h <- if (isTRUE(BLASTAM_USE_BJ_THRESHOLD))
    ifelse(!is.na(out$temp_wet), blastam_bj_min_hours(out$temp_wet),
           as.numeric(BLASTAM_WET_HOURS_FIXED))
  else
    rep(as.numeric(BLASTAM_WET_HOURS_FIXED), nrow(out))

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

# BLASTAM score: count favourable infection days in the rolling window.
# `unjudged` tells you how many days could not be assessed, so a low count
# from missing data is distinguishable from one from unfavourable weather.
blastam_score <- function(infect, semi, dates, end_date,
                          window = 21L, recent = 7L) {
  inwin <- dates > (end_date - window)
  list(events   = sum(infect[inwin], na.rm = TRUE),
       combined = sum(infect[inwin], na.rm = TRUE) +
                  sum(semi[inwin],   na.rm = TRUE),
       recent   = sum(infect[dates > (end_date - recent)], na.rm = TRUE),
       unjudged = sum(is.na(infect[inwin])))
}

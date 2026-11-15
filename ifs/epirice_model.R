################################################################################
# epirice_model.R
#
# The EPIRICE leaf blast model, vendored as a standalone file.
#
# This is the SEIR engine and the leaf blast parameterisation from the epicrop
# R package by Adam H. Sparks and colleagues, which implements the EPIRICE model
# of Savary et al. (2012), Crop Protection 34:6-17,
# doi:10.1016/j.cropro.2011.11.009.
#
#   Original model: Serge Savary & Rene Pangga (IRRI)
#   Original R implementation: Robert J. Hijmans, Rene Pangga, Jorrel Aunario
#   epicrop package + leaf blast parameters: Adam H. Sparks (DPIRD, WA)
#   SEIR framework: Zadoks (1971), doi:10.1094/Phyto-61-600
#
# The SEIR() and .fn_Rc() functions below are reproduced from the epicrop source
# so the model runs without installing the package. Please retain this
# attribution and cite Savary et al. (2012) if you publish or distribute
# results. Check the epicrop licence before redistributing the code itself.
#
# TWO DELIBERATE DEPARTURES FROM THE PACKAGE
#
#  1. .calculate_audpc() uses the standard trapezoidal definition (Madden, Hughes
#     & van den Bosch 2007). It is an internal helper and affects only the AUDPC
#     output column, not the disease dynamics.
#
#     Temperature response evidence behind the default optimum:
#       Kato, H. and Kozaka, T. (1974) Effect of temperature on lesion
#         enlargement and sporulation of Pyricularia oryzae in rice leaves.
#         Phytopathology 64: 828-830. doi:10.1094/Phyto-64-828.
#       Hashioka, Y. (1965) Effects of environmental factors on development of
#         causal fungus, infection, disease development, and epidemiology in
#         rice blast disease. In: The Rice Blast Disease. J Hopkins Press,
#         pp. 153-161.
#
#  2. The RcT infection-optimum temperature is SELECTABLE, via EPIRICE_RCT_PEAK
#     in blast_config.R, and defaults to the published 25 C rather than epicrop's
#     20 C. See section 4b of blast_config.R for the reasoning and for the
#     magnitude of the difference. Previously the README argued at length for the
#     25 C peak while this file shipped the 20 C curve and carried a comment
#     asserting 20 C, so the documented model and the running model disagreed by
#     roughly a factor of two at northern Australian temperatures.
#
# NOTE ON THE SEIR BODY. It is reproduced from upstream, including three
# constructs that look wrong but are load-bearing and match the package:
#   * removed[d] is read on the line above the one that assigns it, so it reads
#     the pre-allocated zero;
#   * sum(infectious) sums the whole pre-allocated vector, giving the cumulative
#     total of infectious increments;
#   * removed_today reads `infday` carried over from the previous iteration.
# Diff against epicrop before changing any of them.
################################################################################

suppressPackageStartupMessages(library(data.table))

# Standard area under the disease progress curve (daily step = 1 day)
.calculate_audpc <- function(intensity) {
  n <- length(intensity)
  if (n < 2L) return(0)
  sum((intensity[-1L] + intensity[-n]) / 2)
}

# Linear interpolation of an Rc modifier curve (age or temperature)
.fn_Rc <- function(.Rc, .xout) {
  stats::approx(
    x = .Rc[, 1],
    y = .Rc[, 2],
    method = "linear",
    xout = .xout,
    yleft = 0,
    yright = 0
  )$y
}

# ---- Temperature response curves -------------------------------------------
# Relative infection rate against daily mean air temperature.
#
#  peak 25: as published, Table 2 of Savary et al. 2012.
#  peak 20: as implemented in epicrop, whose source asserts the table has a typo.
EPIRICE_RCT_CURVES <- list(
  `25` = cbind(c(10L, 15L, 20L, 25L, 30L, 35L, 40L, 45L),
               c(0,   0.5, 0.6, 1.0, 0.6, 0.2, 0.05, 0)),
  `20` = cbind(c(10L, 15L, 20L, 25L, 30L, 35L, 40L, 45L),
               c(0,   0.5, 1.0, 0.6, 0.2, 0.05, 0.01, 0))
)

epirice_rct <- function(peak = if (exists("EPIRICE_RCT_PEAK", inherits = TRUE))
                                EPIRICE_RCT_PEAK else 25L) {
  key <- as.character(as.integer(peak))
  if (!key %in% names(EPIRICE_RCT_CURVES))
    stop("EPIRICE_RCT_PEAK must be 25 or 20, not ", peak, call. = FALSE)
  EPIRICE_RCT_CURVES[[key]]
}

# ---- SEIR engine (reproduced from epicrop::SEIR) --------------------------
SEIR <-
  function(wth,
           emergence,
           onset,
           duration,
           rhlim,
           rainlim,
           H0,
           I0,
           RcA,
           RcT,
           RcOpt,
           p,
           i,
           Sx,
           a,
           RRS,
           RRG) {

    if (!inherits(wth, "data.table")) {
      wth <- as.data.table(wth)
    }

    if (a < 1L) {
      stop(call. = FALSE,
           "`a` cannot be set to less than 1.")
    }
    if (H0 < 0) {
      stop(call. = FALSE, "H0 cannot be < 0.")
    }
    if (I0 < 0) {
      stop(call. = FALSE, "I0 cannot be < 0.")
    }

    # Set and check dates
    emergence <- as.Date(emergence)
    harvest <- emergence + sum(duration, -1)
    dates <- seq(from = emergence, to = harvest, by = "day")

    if (!(emergence >= wth[1, YYYYMMDD]) ||
        (max(dates) > max(wth[, YYYYMMDD]))) {
      stop(call. = FALSE,
           "Incomplete weather data or dates do not align")
    }

    if (nrow(wth) > duration) {
      wth <-
        wth[YYYYMMDD %between% c(emergence, emergence + sum(duration, -1))]
    }

    # The daily loop indexes the weather BY POSITION, so a missing calendar day
    # would shift every later day against crop age and be modelled silently.
    # Both callers screen for this, but a direct caller may not.
    if (nrow(wth) > 1L &&
        any(as.integer(diff(sort(wth$YYYYMMDD))) != 1L)) {
      stop(call. = FALSE,
           "Weather series has gaps; SEIR indexes by position and cannot be trusted")
    }

    # Reference vectors
    wth_rain <- wth$RAIN
    wth_rhum <- wth$RHUM
    Rc_temp <- .fn_Rc(.Rc = RcT, .xout = wth$TEMP)
    Rc_age <- .fn_Rc(.Rc = RcA, .xout = 0:duration)

    cofr <-
      rc <-
      RcW <- latency <- infectious <- intensity <- rsenesced <-
      rgrowth <-
      rtransfer <- infection <- diseased <- senesced <- removed <-
      now_infectious <- now_latent <- sites <- total_sites <-
      vector(mode = "double", length = duration)

    for (d in 1:duration) {
      d_1 <- sum(d, -1)

      if (d == 1) {
        sites[d] <- H0
        rsenesced[d] <- RRS * sites[d]
      } else {
        if (d > i) {
          removed_today <- infectious[infday + 1]
        } else {
          removed_today <- 0
        }

        sites[[d]] <-
          sum(sites[d_1], rgrowth[d_1], -infection[d_1], -rsenesced[d_1])
        rsenesced[[d]] <- sum(removed_today, RRS * sites[d])
        senesced[[d]] <- sum(senesced[d_1], rsenesced[d_1])

        latency[[d]] <- infection[d_1]
        latday <- sum(d, -p)
        latday <- max(1, latday)
        now_latent[[d]] <- sum(latency[latday:d])

        infectious[[d]] <- rtransfer[d_1]
        infday <- sum(d, -i)
        infday <- max(1, infday)
        now_infectious[[d]] <- sum(infectious[infday:d])
      }

      if (wth_rhum[[d]] >= rhlim || wth_rain[[d]] >= rainlim) {
        RcW[[d]] <- 1
      }

      rc[[d]] <- RcOpt * Rc_age[d] * Rc_temp[d] * RcW[d]

      diseased[[d]] <- sum(sum(infectious), now_latent[d], removed[d])
      removed[[d]] <- sum(sum(infectious), -now_infectious[d])
      cofr[[d]] <- 1 - diseased[d] / sum(sites[d], diseased[d])

      if (d == onset) {
        infection[[d]] <- I0
      } else if (d > onset) {
        infection[[d]] <- now_infectious[d] * rc[d] * (cofr[d] ^ a)
      } else {
        infection[[d]] <- 0
      }

      if (d >= p) {
        rtransfer[[d]] <- latency[latday]
      } else {
        rtransfer[[d]] <- 0
      }

      total_sites[[d]] <- sum(diseased[d], sites[d])
      rgrowth[[d]] <- RRG * sites[d] * sum(1, -(total_sites[d] / Sx))
      intensity[[d]] <- sum(diseased[d], -removed[d]) /
        sum(total_sites[d], -removed[d])
    } # end loop

    simday <- seq_len(duration)
    AUDPC <- .calculate_audpc(intensity)

    out <-
      setDT(
        list(
          "simday" = simday,
          "dates" = dates[seq_len(duration)],
          "sites" = sites,
          "latent" = now_latent,
          "infectious" = now_infectious,
          "removed" = removed,
          "senesced" = senesced,
          "rateinf" = infection,
          "rtransfer" = rtransfer,
          "rgrowth" = rgrowth,
          "rsenesced" = rsenesced,
          "diseased" = diseased,
          "intensity" = intensity,
          "AUDPC" = rep_len(AUDPC, duration)
        )
      )

    if (all(c("LAT", "LON") %in% names(wth))) {
      out[, lat := rep_len(wth[, LAT], .N)]
      out[, lon := rep_len(wth[, LON], .N)]
    }

    return(out[])
  }

# ---- Leaf blast parameterisation (from epicrop::predict_leaf_blast) ---------
# Published parameters, Savary et al. (2012) Table 2:
#   onset 15 d, duration 120 d, rhlim 90%, rainlim 5 mm, H0 600, I0 1,
#   RcOpt 1.14, p 5 d, i 20 d, a 1, Sx 30000, RRS 0.01, RRG 0.1.
#
# `duration` is exposed so an in-season weekly run can simulate only as far as
# today's weather. The runners pass CROP_AGE_DAYS + 1 rows, which gives a final
# crop age of CROP_AGE_DAYS.
#
# NOTE ON THE WETNESS GATE. RcW is switched on by a daily MEAN RH of 90% or a
# daily rain SUM of 5 mm. A 24 hour mean of 90% needs an all-day saturated air
# mass: a series with 95% nights and 58% afternoons averages about 75% and never
# opens the RH branch. In practice this configuration is driven almost entirely
# by the 5 mm rain branch, which is why BLASTAM_DAY_CUT_HOUR matters so much.
predict_leaf_blast <- function(wth, emergence, duration = 120L,
                               rct = epirice_rct()) {
  SEIR(
    wth = wth,
    emergence = emergence,
    onset = 15L,
    duration = as.integer(duration),
    rhlim = 90L,
    rainlim = 5L,
    H0 = 600L,
    I0 = 1L,
    RcA = cbind(
      c(0L, 5L, 10L, 15L, 20L, 25L, 30L, 35L, 40L, 45L, 50L, 55L, 60L,
        65L, 70L, 75L, 80L, 85L, 90L, 95L, 100L, 105L, 110L, 115L, 120L),
      c(1, 1, 1, 0.9, 0.8, 0.7, 0.64, 0.59, 0.53, 0.43, 0.32, 0.22, 0.16,
        0.09, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
    ),
    RcT = rct,
    RcOpt = 1.14,
    p = 5L,
    i = 20L,
    a = 1L,
    Sx = 30000L,
    RRS = 0.01,
    RRG = 0.1
  )
}

predict_lb <- predict_leaf_blast

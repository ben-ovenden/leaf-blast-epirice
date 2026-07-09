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
# The SEIR() and .fn_Rc() functions below are reproduced from the epicrop
# source so the model runs without installing the package. predict_leaf_blast()
# carries the exact leaf blast parameters. Please retain this attribution and
# cite Savary et al. (2012) if you publish or distribute results. Check the
# epicrop licence before redistributing the code itself.
#
# Only change from the package: .calculate_audpc() is implemented here with the
# standard trapezoidal definition (Madden, Hughes & van den Bosch 2007), since
# it is an internal helper. It affects only the AUDPC output column, not the
# disease dynamics.
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
          "dates" = dates[1:d],
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

# ---- Leaf blast parameterisation (exact, from epicrop::predict_leaf_blast) --
# Note: the optimum temperature is 20 C. Table 2 of Savary et al. 2012 shows a
# typo of 25 C; the epicrop implementation uses the corrected 20 C, as below in
# the RcT curve peak.
# duration defaults to the published season length of 120 days. It is exposed
# so an in-season weekly run can simulate only as far as today's weather.
predict_leaf_blast <- function(wth, emergence, duration = 120L) {
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
    RcT = cbind(
      c(10L, 15L, 20L, 25L, 30L, 35L, 40L, 45L),
      c(0, 0.5, 1, 0.6, 0.2, 0.05, 0.01, 0)
    ),
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

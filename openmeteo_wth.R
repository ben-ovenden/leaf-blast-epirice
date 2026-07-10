################################################################################
# openmeteo_wth.R
#
# Replacement for epicrop::get_wth(). Fetches daily weather from the Open-Meteo
# historical archive (ERA5 reanalysis, free, no key) and returns it in the
# schema the SEIR model expects:
#
#   YYYYMMDD (Date), DOY (int), TEMP (mean C), RHUM (mean %), RAIN (mm),
#   LAT, LON
#
# The archive has a lag of roughly 5 days behind real time, which is fine for a
# weekly retrospective risk signal. Open-Meteo docs:
# https://open-meteo.com/en/docs/historical-weather-api
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

OPENMETEO_ARCHIVE_URL <- "https://archive-api.open-meteo.com/v1/archive"

# Parse an Open-Meteo daily JSON list into the wth schema. Separated from the
# HTTP call so it can be tested offline.
.parse_openmeteo_daily <- function(daily, lat, lon) {

  if (is.null(daily) || is.null(daily$time) || length(daily$time) == 0) {
    return(NULL)
  }

  dt <- data.table(YYYYMMDD = as.Date(daily$time))
  dt[, DOY := as.integer(format(YYYYMMDD, "%j"))]

  # TEMP: prefer daily mean, fall back to (max + min) / 2
  temp_mean <- daily$temperature_2m_mean
  if (is.null(temp_mean) || all(is.na(temp_mean))) {
    temp_mean <- (as.numeric(daily$temperature_2m_max) +
                    as.numeric(daily$temperature_2m_min)) / 2
  }
  dt[, TEMP := as.numeric(temp_mean)]

  # RHUM: prefer daily mean, fall back to (max + min) / 2
  rh_mean <- daily$relative_humidity_2m_mean
  if (is.null(rh_mean) || all(is.na(rh_mean))) {
    rh_mean <- (as.numeric(daily$relative_humidity_2m_max) +
                  as.numeric(daily$relative_humidity_2m_min)) / 2
  }
  dt[, RHUM := as.numeric(rh_mean)]

  dt[, RAIN := as.numeric(daily$precipitation_sum)]
  dt[, `:=`(LAT = lat, LON = lon)]

  setorder(dt, YYYYMMDD)
  dt[]
}

# Fill short internal gaps so the SEIR daily loop stays continuous. Carries the
# last observation forward, then back-fills any leading NA. Stops if too much is
# missing, since that means the fetch was incomplete.
.fill_gaps <- function(dt, max_missing_frac = 0.1) {
  for (col in c("TEMP", "RHUM", "RAIN")) {
    v <- dt[[col]]
    miss <- sum(is.na(v))
    if (miss > 0 && miss / length(v) > max_missing_frac) {
      stop(sprintf("Too many missing values in %s (%d of %d)",
                   col, miss, length(v)), call. = FALSE)
    }
    if (miss > 0) {
      # carry forward
      for (k in seq_along(v)[-1]) if (is.na(v[k])) v[k] <- v[k - 1]
      # back-fill any leading NA
      for (k in rev(seq_along(v))[-1]) if (is.na(v[k])) v[k] <- v[k + 1]
      dt[[col]] <- v
    }
  }
  dt
}

# Main adapter: fetch daily weather for one location and date range.
get_openmeteo_wth <- function(lat, lon, start_date, end_date,
                              timeout_seconds = 30) {

  url <- sprintf(
    paste0("%s?latitude=%.4f&longitude=%.4f&start_date=%s&end_date=%s",
           "&daily=temperature_2m_mean,temperature_2m_max,temperature_2m_min,",
           "relative_humidity_2m_mean,relative_humidity_2m_max,",
           "relative_humidity_2m_min,precipitation_sum&timezone=auto"),
    OPENMETEO_ARCHIVE_URL, lat, lon,
    format(as.Date(start_date), "%Y-%m-%d"),
    format(as.Date(end_date), "%Y-%m-%d")
  )

  resp <- tryCatch(
    withCallingHandlers(
      jsonlite::fromJSON(url),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) {
      warning("Open-Meteo fetch failed for (", lat, ", ", lon, "): ",
              conditionMessage(e), call. = FALSE)
      NULL
    }
  )

  if (is.null(resp) || is.null(resp$daily)) return(NULL)

  dt <- .parse_openmeteo_daily(resp$daily, lat, lon)
  if (is.null(dt)) return(NULL)

  dt <- .fill_gaps(dt)
  dt
}

################################################################################
# Gridded fetch: many locations per request (Open-Meteo allows up to 1000
# coordinates per call). Returns a list of per-point data.tables in the SEIR
# schema, indexed to the input order. Used to build a risk heatmap.
################################################################################

# Parse one location element (from a simplifyVector = FALSE response) to schema
.parse_openmeteo_point <- function(el, lat, lon) {
  d <- el$daily
  if (is.null(d) || is.null(d$time) || length(d$time) == 0) return(NULL)
  gv <- function(x) if (is.null(x)) NA_real_ else as.numeric(unlist(x))
  tmean <- gv(d$temperature_2m_mean)
  if (all(is.na(tmean))) tmean <- (gv(d$temperature_2m_max) + gv(d$temperature_2m_min)) / 2
  rmean <- gv(d$relative_humidity_2m_mean)
  if (all(is.na(rmean))) rmean <- (gv(d$relative_humidity_2m_max) + gv(d$relative_humidity_2m_min)) / 2
  dates <- as.Date(unlist(d$time))
  dt <- data.table(
    YYYYMMDD = dates,
    DOY = as.integer(format(dates, "%j")),
    TEMP = tmean, RHUM = rmean, RAIN = gv(d$precipitation_sum),
    LAT = lat, LON = lon
  )
  setorder(dt, YYYYMMDD)
  dt[]
}

# Small wrapper around the JSON fetch so it can be stubbed in tests.
.om_fetch <- function(url) {
  tryCatch(jsonlite::fromJSON(url, simplifyVector = FALSE),
           error = function(e) {
             warning("grid fetch failed: ", conditionMessage(e), call. = FALSE)
             NULL
           })
}

# Fetch weather for a set of grid points. lats/lons are equal-length vectors.
#
# Strategy: request coordinates in SMALL batches (keeps the URL short; large
# multi-coordinate URLs are commonly rejected). Any cell that still has no data
# afterwards is retried with a single-point request, which is the proven path
# used elsewhere. A summary line reports how many cells returned data.
get_openmeteo_grid <- function(lats, lons, start_date, end_date,
                               batch = 100L, pause = 1) {
  stopifnot(length(lats) == length(lons))
  n <- length(lats)
  out <- vector("list", n)
  vars <- paste0("temperature_2m_mean,temperature_2m_max,temperature_2m_min,",
                 "relative_humidity_2m_mean,relative_humidity_2m_max,",
                 "relative_humidity_2m_min,precipitation_sum")
  sd <- format(as.Date(start_date), "%Y-%m-%d")
  ed <- format(as.Date(end_date), "%Y-%m-%d")

  idx <- split(seq_len(n), ceiling(seq_len(n) / batch))
  for (g in idx) {
    url <- sprintf(
      "%s?latitude=%s&longitude=%s&start_date=%s&end_date=%s&daily=%s&timezone=UTC",
      OPENMETEO_ARCHIVE_URL,
      paste(sprintf("%.4f", lats[g]), collapse = ","),
      paste(sprintf("%.4f", lons[g]), collapse = ","), sd, ed, vars)
    resp <- .om_fetch(url)
    if (!is.null(resp)) {
      # multi-point response is a top-level array; single point is one object
      elems <- if (!is.null(resp$daily)) list(resp) else resp
      for (j in seq_along(g)) {
        el <- tryCatch(elems[[j]], error = function(e) NULL)
        if (!is.null(el)) {
          p <- tryCatch(.parse_openmeteo_point(el, lats[g[j]], lons[g[j]]),
                        error = function(e) NULL)
          if (!is.null(p)) p <- tryCatch(.fill_gaps(p), error = function(e) NULL)
          out[[g[j]]] <- p
        }
      }
    }
    if (pause > 0) Sys.sleep(pause)
  }

  # Fall back to single-point requests for any cell still missing.
  missing <- which(vapply(out, is.null, logical(1)))
  if (length(missing) > 0) {
    cat(sprintf("  %d of %d cells missing after batched fetch; retrying singly\n",
                length(missing), n))
    for (k in missing) {
      out[[k]] <- tryCatch(
        get_openmeteo_wth(lats[k], lons[k], start_date, end_date),
        error = function(e) NULL)
      if (pause > 0) Sys.sleep(0.15)
    }
  }

  ndata <- sum(!vapply(out, is.null, logical(1)))
  cat(sprintf("  grid weather: %d of %d cells returned data\n", ndata, n))
  out
}

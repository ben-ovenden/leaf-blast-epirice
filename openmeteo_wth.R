################################################################################
# openmeteo_wth.R
#
# Single point weather adapters for the Open-Meteo historical archive (ERA5).
# The GRID path no longer lives here: see openmeteo_batch.R, which fetches many
# locations per request. What remains is used by run_blast.R for the town table
# and by the offline tests.
#
# Returns the schema the SEIR model expects:
#   YYYYMMDD (Date), DOY (int), TEMP (mean C), RHUM (mean %), RAIN (mm), LAT, LON
#
# The archive lags real time by roughly 5 days, which is why ARCHIVE_LAG_DAYS
# exists. Docs: https://open-meteo.com/en/docs/historical-weather-api
#
# TIMEZONE. Requests are made in UTC, deliberately and consistently: timezone=auto
# returns different zones for adjacent cells near state borders, which would make
# the aggregation window discontinuous across the continental map.
#
# The UTC timestamps are then converted to LOCAL SOLAR time inside
# blastam_daily_from_hourly(), which cuts the 24 hour model day at
# BLASTAM_DAY_CUT_HOUR (10:00 local solar) so the nocturnal wet and rain period
# stays inside one day. An earlier version of this note claimed the UTC day was
# itself the mechanism keeping the dew period intact, which was true only for
# eastern longitudes and stopped being true when schema 2 moved the aggregates
# onto local midnight days. Do not read the raw UTC day as the model day.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# On Windows behind a corporate proxy or SSL inspection, route web requests
# through the Windows system proxy and certificate store. No-op on Linux.
if (.Platform$OS.type == "windows") {
  options(url.method = "wininet", download.file.method = "wininet")
}

OPENMETEO_ARCHIVE_URL <- "https://archive-api.open-meteo.com/v1/archive"

.wth_cfg <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default

# Parse an Open-Meteo daily JSON list into the wth schema. Separated from the
# HTTP call so it can be tested offline.
.parse_openmeteo_daily <- function(daily, lat, lon) {

  if (is.null(daily) || is.null(daily$time) || length(daily$time) == 0) return(NULL)

  dt <- data.table(YYYYMMDD = as.Date(daily$time))
  dt[, DOY := as.integer(format(YYYYMMDD, "%j"))]
  n <- nrow(dt)

  # Recycle-safe getter. as.numeric(NULL) is numeric(0), and assigning a
  # zero-length column into a data.table errors rather than returning NULL,
  # which is how a missing variable used to crash the parse.
  gv <- function(x) {
    v <- suppressWarnings(as.numeric(x))
    if (length(v) == 0L) rep(NA_real_, n) else v
  }

  # TEMP: prefer daily mean, fall back to (max + min) / 2
  temp_mean <- gv(daily$temperature_2m_mean)
  if (all(is.na(temp_mean)))
    temp_mean <- (gv(daily$temperature_2m_max) + gv(daily$temperature_2m_min)) / 2
  if (length(temp_mean) != n || all(is.na(temp_mean))) return(NULL)
  dt[, TEMP := temp_mean]

  # RHUM: prefer daily mean, fall back to (max + min) / 2
  rh_mean <- gv(daily$relative_humidity_2m_mean)
  if (all(is.na(rh_mean)))
    rh_mean <- (gv(daily$relative_humidity_2m_max) + gv(daily$relative_humidity_2m_min)) / 2
  dt[, RHUM := rh_mean]

  dt[, RAIN := gv(daily$precipitation_sum)]
  dt[, `:=`(LAT = lat, LON = lon)]

  setorder(dt, YYYYMMDD)
  dt[]
}

# Fill short internal gaps so the SEIR daily loop stays continuous. Carries the
# last observation forward, then back-fills any leading NA. Stops if too much is
# missing, since that means the fetch was incomplete.
#
# The run length guard is new: total missingness under 10% could previously still
# be one unbroken 6 day block carried forward from a single value, which is not a
# short gap in any useful sense.
.fill_gaps <- function(dt, max_missing_frac = 0.1, max_run = 3L) {
  for (col in c("TEMP", "RHUM", "RAIN")) {
    v <- dt[[col]]
    miss <- sum(is.na(v))
    if (miss == 0L) next
    if (miss / length(v) > max_missing_frac)
      stop(sprintf("Too many missing values in %s (%d of %d)", col, miss, length(v)),
           call. = FALSE)
    r <- rle(is.na(v))
    longest <- if (any(r$values)) max(r$lengths[r$values]) else 0L
    if (longest > max_run)
      stop(sprintf("Gap of %d consecutive days in %s exceeds max_run=%d",
                   longest, col, max_run), call. = FALSE)
    for (k in seq_along(v)[-1]) if (is.na(v[k])) v[k] <- v[k - 1]
    for (k in rev(seq_along(v))[-1]) if (is.na(v[k])) v[k] <- v[k + 1]
    set(dt, j = col, value = v)
  }
  dt
}

# Daily adapter for one location and date range.
get_openmeteo_wth <- function(lat, lon, start_date, end_date,
                              timeout_seconds = .wth_cfg("OM_TIMEOUT_S", 60L)) {
  old <- options(timeout = timeout_seconds)
  on.exit(options(old), add = TRUE)

  url <- sprintf(
    paste0("%s?latitude=%.4f&longitude=%.4f&start_date=%s&end_date=%s",
           "&daily=temperature_2m_mean,temperature_2m_max,temperature_2m_min,",
           "relative_humidity_2m_mean,relative_humidity_2m_max,",
           "relative_humidity_2m_min,precipitation_sum&timezone=UTC&models=%s"),
    OPENMETEO_ARCHIVE_URL, lat, lon,
    format(as.Date(start_date), "%Y-%m-%d"),
    format(as.Date(end_date), "%Y-%m-%d"),
    .wth_cfg("OPENMETEO_MODEL", "era5")
  )

  resp <- tryCatch(
    withCallingHandlers(jsonlite::fromJSON(url),
                        warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) {
      warning("Open-Meteo fetch failed for (", lat, ", ", lon, "): ",
              conditionMessage(e), call. = FALSE)
      NULL
    })

  if (is.null(resp) || is.null(resp$daily)) return(NULL)
  dt <- .parse_openmeteo_daily(resp$daily, lat, lon)
  if (is.null(dt)) return(NULL)
  tryCatch(.fill_gaps(dt), error = function(e) NULL)
}

# Parse one location element from a simplifyVector = FALSE response. Retained
# because openmeteo_batch.R has its own hourly equivalent and this is the daily
# one used by the offline tests.
.parse_openmeteo_point <- function(el, lat, lon) {
  d <- el$daily
  if (is.null(d) || is.null(d$time) || length(d$time) == 0) return(NULL)
  n <- length(d$time)
  dates <- as.Date(vapply(d$time, function(e) as.character(e)[1], character(1)))
  # NULL-safe: unlist() drops JSON nulls, which shortens the vector and
  # mis-aligns every later value against the date axis.
  gv <- function(x) {
    if (is.null(x)) return(rep(NA_real_, n))
    v <- suppressWarnings(vapply(x, function(e)
      if (is.null(e) || length(e) == 0L) NA_real_ else as.numeric(e)[1], numeric(1)))
    if (length(v) != n) rep(NA_real_, n) else v
  }
  tmean <- gv(d$temperature_2m_mean)
  if (all(is.na(tmean))) tmean <- (gv(d$temperature_2m_max) + gv(d$temperature_2m_min)) / 2
  rmean <- gv(d$relative_humidity_2m_mean)
  if (all(is.na(rmean))) rmean <- (gv(d$relative_humidity_2m_max) + gv(d$relative_humidity_2m_min)) / 2
  dt <- data.table(YYYYMMDD = dates, DOY = as.integer(format(dates, "%j")),
                   TEMP = tmean, RHUM = rmean, RAIN = gv(d$precipitation_sum),
                   LAT = lat, LON = lon)
  setorder(dt, YYYYMMDD)
  dt[]
}

################################################################################
# Hourly fetch, one point. Feeds BLASTAM its leaf wetness hours and is aggregated
# to daily values for EPIRICE, so a single fetch serves both models.
#
# NOTE. For more than a handful of points use fetch_points_batched() in
# openmeteo_batch.R instead. Multi location requests do not reduce the weighted
# quota, because weight scales with locations, but the quota was never what
# limited wall clock here: request latency was, and batching removes it.
################################################################################
get_openmeteo_hourly <- function(lat, lon, start_date, end_date,
                                 timeout_seconds = .wth_cfg("OM_TIMEOUT_S", 60L)) {
  old <- options(timeout = timeout_seconds)
  on.exit(options(old), add = TRUE)
  vars <- .wth_cfg("OM_HOURLY_VARS",
                   c("temperature_2m", "relative_humidity_2m", "precipitation"))
  url <- sprintf(
    "%s?latitude=%.4f&longitude=%.4f&start_date=%s&end_date=%s&hourly=%s&timezone=UTC&models=%s",
    OPENMETEO_ARCHIVE_URL, lat, lon,
    format(as.Date(start_date), "%Y-%m-%d"), format(as.Date(end_date), "%Y-%m-%d"),
    paste(vars, collapse = ","), .wth_cfg("OPENMETEO_MODEL", "era5"))
  resp <- tryCatch(
    withCallingHandlers(jsonlite::fromJSON(url),
                        warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) NULL)
  if (is.null(resp) || is.null(resp$hourly) || is.null(resp$hourly$time)) return(NULL)
  h <- resp$hourly
  nh <- length(h$time)
  # Length-safe. jsonlite may return a short vector if the response is ragged,
  # and recycling a short column against the time axis mis-aligns every later
  # value, which is corruption rather than a visible failure.
  gv <- function(x) {
    if (is.null(x)) return(rep(NA_real_, nh))
    v <- suppressWarnings(as.numeric(x))
    if (length(v) != nh) rep(NA_real_, nh) else v
  }
  dt <- as.POSIXct(h$time, format = "%Y-%m-%dT%H:%M", tz = "UTC")
  out <- data.table(dt = dt, temp = gv(h$temperature_2m),
                    rh = gv(h$relative_humidity_2m), rain = gv(h$precipitation))
  out <- out[!is.na(dt)]
  if (nrow(out) == 0L || all(is.na(out$temp)) || all(is.na(out$rh))) return(NULL)
  out
}

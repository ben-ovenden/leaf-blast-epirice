################################################################################
# openmeteo_batch.R
#
# Batched fetch path for run_blast_grid.R. Replaces the per point mclapply loop.
#
# Three changes, each fixing a measured bug:
#
#  1. MANY LOCATIONS PER REQUEST. Open-Meteo accepts comma separated latitude and
#     longitude lists. Weighted cost is unchanged (it scales with locations), but
#     HTTP round trips fall by a factor of OM_BATCH_SIZE. That removes the two
#     things that actually governed wall clock: per request latency, and the
#     mclapply chunk barrier, where one slow request held up its whole chunk.
#
#  2. REAL HTTP STATUS CODES. jsonlite::fromJSON(url) collapses a 429 quota
#     rejection, a 5xx and a socket timeout into one anonymous error, so the run
#     could not tell "we are out of quota" from "the network hiccuped". curl
#     exposes status and body, so quota exhaustion stops the run in seconds
#     instead of burning three hours confirming itself.
#
#  3. TOKEN BUCKET PACING ON WEIGHTED CALLS. The old pacer forced every chunk to
#     last (conc * cost) / rate seconds. Points per chunk is conc, so points per
#     second was rate / (60 * cost) with the conc cancelling exactly: throughput
#     was algebraically independent of concurrency, while the failure rate was
#     not. Concurrency is therefore gone entirely.
#
# Three more, after the September 2026 outage:
#
#  4. A MINUTELY OR HOURLY 429 IS WAITED OUT, NOT TREATED AS FATAL. The 2026-09-07
#     grid run hit the hourly ceiling at exactly 5,000 weighted, abandoned the
#     remaining 42% of the grid with 140 minutes of deadline unused, and the town
#     fetch in the same hour was refused too. Only the DAILY ceiling ends a run.
#
#  5. RETRIES ARE PACED. A retried request was sent straight away, so the server
#     saw more calls per hour than the token bucket accounted for.
#
#  6. COHORT REFRESH PLANNING (om_plan_refresh). Each cached cell is fetched from
#     the day after its newest cached row, so a cell 16 days behind costs 1.64
#     weighted rather than the 4.86 of a full crop-window refetch.
#
# WEIGHTED COST MODEL
#   weight_per_location = max(1, n_variables / 10) * max(1, n_days / 14)
# The 14 day floor is why a one day top up costs the same as a fourteen day one.
# Verify against the weight note the Open-Meteo docs UI prints for a given
# request URL before trusting the numbers in a paper.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(curl)
})

OPENMETEO_ARCHIVE_URL <- "https://archive-api.open-meteo.com/v1/archive"

# Defined here rather than after om_request(), which uses it. R resolves the
# symbol at call time so the old order worked, but only if the whole file was
# sourced before the first call.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Defined here rather than after om_request(), which used it. R resolves the
# symbol at call time so it worked, but only if the whole file was sourced.

.cfg <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default

# ---- weighted cost ---------------------------------------------------------
# Vectorised over n_days, so a planner can price many cohorts in one call.
om_weight_per_location <- function(n_days, n_vars = length(.cfg("OM_HOURLY_VARS", c("a","b","c")))) {
  pmax(1, n_vars / 10) * pmax(1, n_days / 14)
}

# ---- token bucket ----------------------------------------------------------
# Holds the sustained weighted rate below the hourly ceiling. The per minute
# ceiling (600) is roughly seven times the sustained hourly one (5,000/h =
# 83/min), so the hourly limit is what binds and the bucket is sized to it.
om_pacer <- function(rate_per_min = 80) {
  last <- Sys.time()
  allowance <- rate_per_min
  function(weight) {
    now <- Sys.time()
    allowance <<- min(rate_per_min,
                      allowance + as.numeric(difftime(now, last, units = "mins")) * rate_per_min)
    last <<- now
    if (allowance < weight) {
      Sys.sleep((weight - allowance) / rate_per_min * 60)
      allowance <<- 0
      last <<- Sys.time()
    } else {
      allowance <<- allowance - weight
    }
    invisible(NULL)
  }
}

# ---- which ceiling did a 429 hit? -------------------------------------------
# Open-Meteo names the ceiling in the body: "Minutely API request limit
# exceeded...", "Hourly API request limit exceeded. Please try again in the next
# hour.", "Daily API request limit exceeded. Please try again tomorrow." An
# unrecognised body is treated like an hourly one: a BOUNDED wait costs wall
# clock, which the run has, while treating it as daily abandons the run, which is
# the failure this replaces. The body is logged either way.
om_quota_kind <- function(msg) {
  m <- tolower(paste(msg %||% "", collapse = " "))
  if (grepl("minute", m, fixed = TRUE)) return("minute")
  if (grepl("hour", m, fixed = TRUE))   return("hour")
  if (grepl("daily|\\bday\\b|tomorrow", m)) return("day")
  "unknown"
}

################################################################################
# Archive edge: which day can the archive actually serve?
#
# om_archive_edge_from_hourly() is the pure part: the newest UTC day whose first
# `need_hours` hours are all present. Those are the hours of the last fetched day
# that the westernmost cell needs to complete its final model day; ERA5 arrives a
# whole UTC day at a time, so in practice this is "the newest published day".
#
# om_probe_archive_edge() makes the one request. See probe_archive_edge.R for
# why: the arithmetic edge names a day that is not published yet for the first
# part of every UTC day, and the western cells then land a day short.
################################################################################
om_archive_edge_from_hourly <- function(hw, need_hours = 3L) {
  if (is.null(hw) || nrow(hw) == 0L) return(as.Date(NA))
  h <- data.table::copy(as.data.table(hw))
  h[, `:=`(day = as.Date(format(dt, "%Y-%m-%d", tz = "UTC")),
           hr  = as.integer(format(dt, "%H", tz = "UTC")))]
  ok <- h[hr < need_hours,
          .(n = sum(!is.na(temp) & !is.na(rh) & !is.na(rain))), by = day][n >= need_hours, day]
  if (length(ok) == 0L) as.Date(NA) else max(ok)
}

om_probe_archive_edge <- function(edge, lon, lat, need_hours = 3L, look_back = 6L) {
  edge <- as.Date(edge)
  res <- om_request(lat, lon, edge - as.integer(look_back), edge)
  if (!identical(res$status, "ok"))
    return(list(date = as.Date(NA), status = res$status, msg = res$msg %||% ""))
  el <- if (!is.null(res$body$hourly)) res$body else res$body[[1]]
  hw <- tryCatch(.om_hourly_dt(el), error = function(e) NULL)
  list(date = om_archive_edge_from_hourly(hw, need_hours), status = "ok", msg = "")
}

################################################################################
# Refresh planner: one cohort per "how far behind", each paying only for the
# days it is missing.
#
#   last_by       data.table(pid, last): each cached cell's newest cached day
#   end_date      last day to be MODELLED this run
#   data_end      last day to be FETCHED (end_date + DAY_CUT_LAG_DAYS)
#   add_keep_from first day a brand new cell keeps; no cohort starts earlier, so
#                 a cell can never cost more than a new one
#   lead          BLASTAM_LEADIN_DAYS, fetched and discarded
#   stale_days    REFRESH_MIN_STALE_DAYS
#   budget, max_n weighted and count caps; `take` marks the affordable prefix
#
# Returns one row per eligible cell, CHEAPEST FIRST (newest `last` first within
# a cost), with keep_from, fetch_from, n_days, missing, cost and take.
#
# WHY. The previous planner refreshed cells with a FIXED 7 day tail and refetched
# anything further behind over the full crop window at 4.86 weighted. After the
# 2026-09-07 run was cut short, 42% of the grid was more than a tail behind, the
# weekly budget bought only ~800 of those cells per run at 4.86, and the map sat
# on the same 29 August window for three weeks. Fetched from the day after its
# newest row, a cell 16 days behind costs 1.64, and that whole backlog would have
# fitted inside one run's budget.
#
# Cheapest first, because GRID_WINDOW_MIN_COVERAGE counts cells: a weighted call
# buys the most coverage when spent on the cells missing the fewest days. What is
# left over is the oldest cohort, and it costs 0.5 more next week, so a chronic
# shortfall wants a midweek top-up run rather than a different order.
################################################################################
om_plan_refresh <- function(last_by, end_date, data_end, add_keep_from,
                            lead = .cfg("BLASTAM_LEADIN_DAYS", 6L),
                            n_vars = length(.cfg("OM_HOURLY_VARS", c("a", "b", "c"))),
                            stale_days = 0L, budget = Inf, max_n = Inf) {
  end_date <- as.Date(end_date); data_end <- as.Date(data_end)
  add_keep_from <- as.Date(add_keep_from)
  el <- as.data.table(last_by)[, .(pid = as.character(pid), last = as.Date(last))]
  el <- el[last < end_date & last <= (end_date - as.integer(stale_days))]
  el[, keep_from := last + 1L]
  el[keep_from < add_keep_from, keep_from := add_keep_from]
  el[, fetch_from := keep_from - as.integer(lead)]
  el[, n_days  := as.integer(data_end - fetch_from) + 1L]
  el[, missing := as.integer(end_date - last)]
  el[, cost    := om_weight_per_location(n_days, n_vars)]
  setorder(el, cost, -last)
  el[, take := cumsum(cost) <= budget & seq_len(.N) <= max_n]
  el[]
}

# ---- one batched request ---------------------------------------------------
# Returns list(status, code, retry_after, body, msg) where status is one of
# "ok", "quota", "http", "transport".
om_request <- function(lats, lons, start_date, end_date, timeout_s = 60) {
  vars  <- .cfg("OM_HOURLY_VARS", c("temperature_2m", "relative_humidity_2m", "precipitation"))
  model <- .cfg("OPENMETEO_MODEL", "era5")
  url <- sprintf(
    "%s?latitude=%s&longitude=%s&start_date=%s&end_date=%s&hourly=%s&timezone=UTC&models=%s",
    OPENMETEO_ARCHIVE_URL,
    paste(sprintf("%.4f", lats), collapse = ","),
    paste(sprintf("%.4f", lons), collapse = ","),
    format(as.Date(start_date), "%Y-%m-%d"),
    format(as.Date(end_date), "%Y-%m-%d"),
    paste(vars, collapse = ","),
    model)

  h <- curl::new_handle(timeout = timeout_s, connecttimeout = 20,
                        accept_encoding = "gzip",
                        useragent = "leaf-blast-epirice (non-commercial research)")
  res <- tryCatch(curl::curl_fetch_memory(url, handle = h),
                  error = function(e) structure(list(msg = conditionMessage(e)),
                                                class = "om_transport_error"))
  if (inherits(res, "om_transport_error"))
    return(list(status = "transport", code = NA_integer_, retry_after = NA_real_,
                body = NULL, msg = res$msg))

  short <- function(x) {
    s <- tryCatch(rawToChar(utils::head(x, 160L)), error = function(e) "")
    substr(gsub("\\s+", " ", trimws(s)), 1L, 120L)
  }
  if (identical(as.integer(res$status_code), 429L)) {
    hdr <- tryCatch(curl::parse_headers_list(res$headers), error = function(e) list())
    ra <- suppressWarnings(as.numeric(hdr[["retry-after"]]))
    return(list(status = "quota", code = 429L,
                retry_after = if (isTRUE(is.finite(ra))) ra else NA_real_,
                body = NULL, msg = short(res$content)))
  }
  if (res$status_code >= 400L)
    return(list(status = "http", code = as.integer(res$status_code),
                retry_after = NA_real_, body = NULL, msg = short(res$content)))

  parsed <- tryCatch(jsonlite::fromJSON(rawToChar(res$content), simplifyVector = FALSE),
                     error = function(e) NULL)
  if (is.null(parsed))
    return(list(status = "http", code = as.integer(res$status_code),
                retry_after = NA_real_, body = NULL, msg = "unparseable JSON"))
  # Open-Meteo signals some argument errors with 200 plus an error object.
  if (isTRUE(parsed$error))
    return(list(status = "http", code = 200L, retry_after = NA_real_, body = NULL,
                msg = as.character(parsed$reason %||% "error flag set")))
  list(status = "ok", code = as.integer(res$status_code), retry_after = NA_real_,
       body = parsed, msg = "")
}


# ---- hourly JSON for one location -> hourly data.table ---------------------
# Indexed BY POSITION, never by the returned latitude and longitude: Open-Meteo
# echoes the grid CELL CENTRE, not the requested coordinate, so matching on
# coordinates would silently mis-assign points. Response order follows the input.
.om_hourly_dt <- function(el) {
  h <- el$hourly
  if (is.null(h) || is.null(h$time) || length(h$time) == 0L) return(NULL)
  nh <- length(h$time)
  # NULL-SAFE. Open-Meteo emits JSON null for a missing hour, and unlist() DROPS
  # nulls rather than preserving them as NA. Using unlist() therefore returns a
  # short vector and silently mis-aligns every later value against the time axis,
  # which is corruption rather than a visible failure. Map element by element.
  gv <- function(x) {
    if (is.null(x)) return(rep(NA_real_, nh))
    v <- suppressWarnings(vapply(x, function(e)
      if (is.null(e) || length(e) == 0L) NA_real_ else as.numeric(e)[1], numeric(1)))
    if (length(v) != nh) rep(NA_real_, nh) else v
  }
  tt <- as.POSIXct(vapply(h$time, function(e) as.character(e)[1], character(1)),
                   format = "%Y-%m-%dT%H:%M", tz = "UTC")
  out <- data.table(dt = tt, temp = gv(h$temperature_2m),
                    rh = gv(h$relative_humidity_2m), rain = gv(h$precipitation))
  out <- out[!is.na(dt)]
  if (nrow(out) == 0L) return(NULL)
  # A location can be structurally valid but entirely empty. Treat that as no
  # data rather than a fetch failure, so it is not retried forever.
  #
  # HUMIDITY IS CHECKED TOO. Only temperature used to be checked, so a response
  # with temperature but a null relative_humidity_2m column passed through and
  # BLASTAM scored every night as not favourable. On the map that is
  # indistinguishable from genuinely dry weather. Rain may legitimately be all
  # zero, so it is not tested for presence.
  if (all(is.na(out$temp)) || all(is.na(out$rh))) return(NULL)
  out
}

################################################################################
# Main entry point.
#
#   pts        data.table with columns pid, lon, lat
#   start_date single Date shared by the whole call. The grid runner groups its
#              cells into cohorts by their newest cached day and makes one call
#              per cohort (om_plan_refresh); the town runner has one window.
#   end_date   single Date
#   on_point   function(pid, lon, lat, hourly_dt) -> rows to append, or NULL
#   deadline   POSIXct wall clock stop, or NA
#   budget     remaining weighted calls for this run
#   pacer      an om_pacer() closure shared across calls, so several cohort calls
#              in one run draw on ONE token bucket instead of each opening with a
#              full one. NULL makes a private pacer.
#   quota_wait_max_s  total seconds this call may spend waiting out minutely or
#              hourly HTTP 429s before giving up. A daily 429 is never waited out.
#   batch_size locations per request. The town runner's second pass uses 1, so a
#              town that failed inside a batch is retried on its own through the
#              same budget, pacer and ledger rather than through a side path.
#
# Returns list(rows, ledger, spent, stopped, n_ok, waited_s).
#   stopped is "" | "deadline" | "budget" | "quota"
#   ledger  is pid, status ("ok"|"empty"|"http"|"transport"|"quota"), code
################################################################################
fetch_points_batched <- function(pts, start_date, end_date, on_point,
                                 deadline = as.POSIXct(NA),
                                 budget = Inf,
                                 label = "fetch",
                                 pacer = NULL,
                                 quota_wait_max_s = .cfg("OM_QUOTA_WAIT_MAX_MIN", 60) * 60,
                                 batch_size = .cfg("OM_BATCH_SIZE", 25L)) {
  n <- nrow(pts)
  if (n == 0L)
    return(list(rows = list(), ledger = data.table(pid = character(), status = character(),
                                                   code = integer()),
                spent = 0, stopped = "", n_ok = 0L, waited_s = 0))

  batch_size <- max(1L, as.integer(batch_size))
  rate       <- .cfg("GRID_TARGET_PER_MIN", 80)
  timeout_s  <- as.integer(.cfg("OM_TIMEOUT_S", 60L))
  max_att    <- as.integer(.cfg("OM_MAX_ATTEMPTS", 3L))
  wait_s_cfg <- as.numeric(.cfg("OM_QUOTA_WAIT_S", 300))

  n_days <- as.integer(as.Date(end_date) - as.Date(start_date)) + 1L
  w_loc  <- om_weight_per_location(n_days)
  pace   <- if (is.function(pacer)) pacer else om_pacer(rate)

  rows <- vector("list", n)
  n_rows <- 0L
  spent <- 0; stopped <- ""; waited <- 0; n_waits <- 0L
  ledger <- data.table(pid = pts$pid, status = NA_character_, code = NA_integer_)
  t0 <- Sys.time(); n_ok <- 0L; n_req <- 0L; n_logged <- 0L
  max_log <- 5L    # a wholesale outage otherwise fills the log with one line per batch
  tally <- c(ok = 0L, empty = 0L, quota = 0L, http = 0L, transport = 0L)

  cat(sprintf("  %s: %d points, %d days, %.2f weighted each, batches of %d, pacing %.0f weighted/min\n",
              label, n, n_days, w_loc, batch_size, rate))

  i <- 1L
  while (i <= n) {
    if (!is.na(deadline) && Sys.time() > deadline) { stopped <- "deadline"; break }
    j <- min(i + batch_size - 1L, n)
    idx <- i:j
    w_batch <- length(idx) * w_loc
    if (spent + w_batch > budget) { stopped <- "budget"; break }

    attempt <- 1L
    repeat {
      # EVERY attempt goes through the pacer. Retries used to be sent straight
      # away, so a flaky stretch pushed the server-side rate above the hourly
      # ceiling while the local bucket still believed it was under it.
      pace(w_batch)
      n_req <- n_req + 1L
      res <- om_request(pts$lat[idx], pts$lon[idx], start_date, end_date, timeout_s)
      # EVERY attempt is charged, not just the successful one. A 429 is a
      # rejection rather than work performed, so it alone is free; a 400 or a
      # timeout may well have cost quota server side, and assuming it did is the
      # safe direction. Charging only the last attempt let a flaky run make three
      # times the requests it accounted for and walk into a 429.
      if (res$status != "quota") spent <- spent + w_batch
      if (res$status == "ok") break
      if (res$status == "quota") break              # decided below: wait, or stop
      if (attempt >= max_att) break
      if (spent + w_batch > budget) break           # retries must not overrun it
      # Jittered exponential backoff, transport errors and 5xx only.
      Sys.sleep(min(30, 2^attempt) * runif(1, 0.5, 1.5))
      attempt <- attempt + 1L
    }

    if (res$status == "quota") {
      # WHICH ceiling? Only the daily one ends the run. A minutely or hourly
      # rejection is waited out, then the SAME batch is sent again (i is not
      # advanced), within this call's wait allowance and the run's deadline. The
      # 2026-09-07 run stopped here at exactly 5,000 weighted, one hour in, with
      # 140 minutes of deadline unused and 42% of the grid unrefreshed.
      kind <- om_quota_kind(res$msg)
      wait_s <- if (isTRUE(is.finite(res$retry_after)) && res$retry_after > 0) res$retry_after
                else if (kind == "minute") 65 else wait_s_cfg
      past_deadline <- !is.na(deadline) && (Sys.time() + wait_s) >= deadline
      can_wait <- kind != "day" && (waited + wait_s) <= quota_wait_max_s && !past_deadline
      if (can_wait) {
        n_waits <- n_waits + 1L
        cat(sprintf("  %s: HTTP 429 (%s): %s ceiling. Waiting %.0f s, then resuming at batch %d-%d (waited %.0f of at most %.0f s this call).\n",
                    label, res$msg, kind, wait_s, i, j, waited, quota_wait_max_s))
        Sys.sleep(wait_s)
        waited <- waited + wait_s
        next
      }
      why <- if (kind == "day") "the DAILY ceiling, which waiting cannot fix"
             else if (past_deadline) "waiting would pass the fetch deadline"
             else sprintf("already waited %.0f s of the %.0f s allowed", waited, quota_wait_max_s)
      cat(sprintf("  %s: HTTP 429 (%s): %s. Stopping fetches for this run.\n",
                  label, res$msg, why))
      tally["quota"] <- tally["quota"] + 1L
      ledger[idx, `:=`(status = "quota", code = 429L)]
      stopped <- "quota"
      break
    }

    if (res$status != "ok") {
      tally[res$status] <- tally[res$status] + 1L
      ledger[idx, `:=`(status = res$status, code = res$code)]
      n_logged <- n_logged + 1L
      if (n_logged <= max_log)
        cat(sprintf("  %s: batch %d-%d failed after %d attempts (%s, code %s): %s\n",
                    label, i, j, attempt, res$status, as.character(res$code), res$msg))
      else if (n_logged == max_log + 1L)
        cat(sprintf("  %s: further batch failures suppressed; see the status tally below.\n", label))
      i <- j + 1L
      next
    }

    body <- res$body
    # One location returns an object; several return an array. Normalise.
    els <- if (!is.null(body$hourly)) list(body) else body
    if (length(els) != length(idx))
      cat(sprintf("  %s: WARNING batch %d-%d asked for %d locations, got %d; the surplus are marked empty.\n",
                  label, i, j, length(idx), length(els)))
    for (k in seq_along(idx)) {
      m <- idx[k]
      el <- if (k <= length(els)) els[[k]] else NULL
      hw <- if (is.null(el)) NULL else .om_hourly_dt(el)
      r  <- if (is.null(hw)) NULL else
        tryCatch(on_point(pts$pid[m], pts$lon[m], pts$lat[m], hw), error = function(e) NULL)
      if (!is.null(r) && nrow(r) > 0) {
        n_rows <- n_rows + 1L
        rows[[n_rows]] <- r
        ledger[m, `:=`(status = "ok", code = 200L)]
        n_ok <- n_ok + 1L; tally["ok"] <- tally["ok"] + 1L
      } else {
        ledger[m, `:=`(status = "empty", code = 200L)]
        tally["empty"] <- tally["empty"] + 1L
      }
    }
    if (j %% 500 < batch_size || j == n)
      cat(sprintf("  %s %d/%d (%d ok, %.0f weighted spent)\n", label, j, n, n_ok, spent))
    i <- j + 1L
  }

  length(rows) <- n_rows
  el_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("  %s done: %d ok of %d attempted in %.1f min (%.1f pts/min) over %d requests; ~%.0f weighted%s%s\n",
              label, n_ok, min(i - 1L, n), el_min,
              if (el_min > 0) n_ok / el_min else NA_real_, n_req, spent,
              if (n_waits > 0) sprintf("; waited out %d HTTP 429(s) for %.0f s", n_waits, waited) else "",
              if (nzchar(stopped)) sprintf("; stopped early: %s", stopped) else ""))
  cat(sprintf("  %s status: %s\n", label,
              paste(sprintf("%s=%d", names(tally), tally), collapse = " ")))
  list(rows = rows, ledger = ledger[!is.na(status)], spent = spent,
       stopped = stopped, n_ok = n_ok, waited_s = waited)
}

################################################################################
# Shared weighted-spend ledger
#
# DAILY_WEIGHTED_CAP was applied per RUN, not per day, and run_blast.R fetched
# its 31 towns with budget = Inf on top of whatever the grid run had already
# spent. The 2026-07-30 grid run reported 8,667 weighted calls against a 9,000
# cap, and the town run then added roughly 150 more, unbudgeted.
#
# Both scripts now append their spend here, keyed on the UTC day the quota resets
# on, and read it back before deciding their own budget.
################################################################################
# The ledger day is the run's PINNED UTC date (BLAST_UTC_DATE, via blast_utc_date()
# when blast_config.R is loaded), not the clock at the moment of writing. The
# Sunday-evening grid run crosses 00:00 UTC, so the wall clock booked its whole
# 8,550 to Monday: the 2026-09-21 evening rerun then read "8,700 already spent
# today" and was capped to 800, for spend the API had largely counted against
# Sunday. Keyed on the day the run started, a rerun on the same pinned day still
# shares the budget, which is the property the ledger exists for.
om_spend_utc_day <- function() {
  if (exists("blast_utc_date", mode = "function")) return(as.Date(blast_utc_date()))
  as.Date(format(Sys.time(), "%Y-%m-%d", tz = "UTC"))
}

om_spend_read <- function(file, day = om_spend_utc_day()) {
  if (!file.exists(file)) return(0)
  d <- tryCatch(fread(file), error = function(e) NULL)
  if (is.null(d) || !all(c("utc_day", "spent") %in% names(d))) return(0)
  s <- suppressWarnings(sum(d[as.Date(utc_day) == as.Date(day), spent], na.rm = TRUE))
  if (!is.finite(s)) 0 else s
}

om_spend_add <- function(file, spent, label, day = om_spend_utc_day()) {
  row <- data.table(utc_day = format(as.Date(day)), stamp = format(Sys.time(), tz = "UTC"),
                    label = label, spent = round(as.numeric(spent), 1))
  old <- if (file.exists(file))
    tryCatch(fread(file, colClasses = list(character = c("utc_day", "stamp", "label"))),
             error = function(e) NULL) else NULL
  d <- if (is.null(old) || !all(c("utc_day", "spent") %in% names(old))) row else
    rbind(old, row, fill = TRUE)
  # Keep a fortnight; the ledger only has to survive one quota window.
  d <- d[as.Date(utc_day) >= (as.Date(day) - 14L)]
  fwrite(d, file)
  invisible(d)
}

#!/usr/bin/env Rscript
################################################################################
# run_blast_grid.R
#
# Two continental heatmaps from a shared, cached weather record with DYNAMIC
# resolution:
#   1) EPIRICE leaf blast intensity (mechanistic epidemic, % leaf area diseased)
#   2) BLASTAM infection days (count of days favourable for new infection)
#
# Both are derived from ONE hourly Open-Meteo fetch per point: hourly data gives
# BLASTAM its leaf-wetness days and is aggregated to daily values for EPIRICE.
# The cache holds daily rows (pruned to the window); each run only fetches the
# newest days for existing points and adds new points, so the maps sharpen over
# successive weeks. See blastam_model.R for BLASTAM parameters and caveats.
################################################################################

# Dates in Australian Eastern time (AEST/AEDT), matching run_blast.R, so the
# run tag and map filenames agree with the Monday-morning email.
Sys.setenv(TZ = "Australia/Sydney")

# Start of the run: the fetch deadline is measured from here, so time spent loading
# the cache and building the target grid counts against the budget too.
RUN_T0 <- Sys.time()

SCRIPT_DIR <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile), winslash = "/"),
  error = function(e) normalizePath(getwd(), winslash = "/"))

source(file.path(SCRIPT_DIR, "epirice_model.R"))
source(file.path(SCRIPT_DIR, "blastam_model.R"))
source(file.path(SCRIPT_DIR, "openmeteo_wth.R"))
source(file.path(SCRIPT_DIR, "blast_config.R"))

suppressPackageStartupMessages({library(data.table); library(terra)})

# One hourly fetch -> daily rows carrying both models' inputs.
if (!exists("POINT_FETCH_FN")) POINT_FETCH_FN <- function(lat, lon, start_date, end_date) {
  to <- if (exists("GRID_FETCH_TIMEOUT_S")) GRID_FETCH_TIMEOUT_S else 45
  hw <- get_openmeteo_hourly(lat, lon, start_date, end_date, timeout_seconds = to)
  if (is.null(hw) || nrow(hw) == 0) return(NULL)
  blastam_daily_from_hourly(hw)   # date, TEMP, RHUM, RAIN, wet_hours, temp_wet, infect
}

ext <- GRID_EXTENT
key_of <- function(lon, lat) paste0(sprintf("%.4f", lon), "_", sprintf("%.4f", lat))
CACHE_COLS <- c("pid","lon","lat","date","TEMP","RHUM","RAIN","wet_hours","temp_wet","infect","semi")

# ---- Land polygon (bundled, terra only) -----------------------------------
get_land <- function() {
  f <- file.path(SCRIPT_DIR, "australia_land.geojson")
  v <- tryCatch(if (file.exists(f)) terra::vect(f) else NULL, error = function(e) NULL)
  if (!is.null(v) && nrow(v) > 0) return(v)
  tryCatch({
    if (requireNamespace("rnaturalearth", quietly = TRUE))
      terra::vect(rnaturalearth::ne_countries(country = "Australia", returnclass = "sf"))
    else NULL
  }, error = function(e) NULL)
}
land_poly <- if (isTRUE(LAND_ONLY)) get_land() else NULL
if (isTRUE(LAND_ONLY) && is.null(land_poly))
  stop("No Australia polygon (australia_land.geojson missing?).", call. = FALSE)
if (!is.null(land_poly)) land_poly <- terra::crop(land_poly, terra::ext(ext))

# ---- Target lattice, coarse-to-fine ---------------------------------------
build_targets <- function() {
  fin <- GRID_RES_FINEST
  g <- as.data.table(expand.grid(lon = seq(ext[1], ext[2], by = fin),
                                 lat = seq(ext[3], ext[4], by = fin)))
  if (!is.null(land_poly)) {
    pts <- terra::vect(as.matrix(g[, .(lon, lat)]), type = "points", crs = "EPSG:4326")
    g <- g[!is.na(terra::extract(land_poly, pts)[, 2])]
  }
  lv <- rep(NA_integer_, nrow(g))
  for (li in seq_along(GRID_RES_LEVELS)) {
    R <- GRID_RES_LEVELS[li]
    on <- (abs((g$lon - ext[1]) / R - round((g$lon - ext[1]) / R)) < 1e-6) &
          (abs((g$lat - ext[3]) / R - round((g$lat - ext[3]) / R)) < 1e-6)
    lv[is.na(lv) & on] <- li
  }
  lv[is.na(lv)] <- length(GRID_RES_LEVELS) + 1L
  g[, lvl := lv]; setorder(g, lvl, lat, lon); g[, pid := key_of(lon, lat)]; g[]
}
targets <- build_targets()

# ---- Window ----------------------------------------------------------------
end_date  <- Sys.Date() - ARCHIVE_LAG_DAYS
emergence <- end_date - CROP_AGE_DAYS
win_dates <- seq(emergence, end_date, by = "day")
cat(sprintf("Window %s to %s (%d days); target %d land points at %.2f deg\n",
            emergence, end_date, length(win_dates), nrow(targets), GRID_RES_FINEST))

# ---- Load cache (rebuild if it predates the BLASTAM columns) ---------------
OUT <- file.path(SCRIPT_DIR, OUTPUT_DIR)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
gz_file  <- file.path(OUT, WEATHER_CACHE_GZ)
csv_file <- file.path(OUT, WEATHER_CACHE_CSV)

# Read the cache preferring the small gz; fall back to the plain CSV if the gz is
# missing or unreadable (also covers the one-time migration from a CSV-only repo).
#
# NOTE: fread() cannot read .gz directly unless the R.utils package is installed
# (it is not in the runner container), so gz files are decompressed with base R's
# gzfile connection into a temp CSV first. That keeps the cache dependency-free.
read_gz_dt <- function(f) {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  inp <- gzfile(f, "rb"); outp <- file(tmp, "wb")
  repeat {
    chunk <- readBin(inp, "raw", 1048576L)
    if (length(chunk) == 0L) break
    writeBin(chunk, outp)
  }
  close(inp); close(outp)
  fread(tmp, colClasses = list(character = "pid"))
}
read_cache <- function() {
  need <- c("pid","lon","lat","date","TEMP","RHUM","RAIN","wet_hours","temp_wet","infect","semi")
  read_one <- function(f) {
    if (!file.exists(f)) return(NULL)
    is_gz <- grepl("\\.gz$", f)
    d <- tryCatch(if (is_gz) read_gz_dt(f) else fread(f, colClasses = list(character = "pid")),
                  error = function(e) {
                    cat(sprintf("Cache read FAILED (%s): %s\n", basename(f), conditionMessage(e))); NULL })
    # A damaged file can read without error yet be empty or missing columns (a
    # non-gzip file opened through gzfile reads as plain text, for instance).
    if (is.null(d) || nrow(d) == 0 || !all(need %in% names(d))) {
      if (!is.null(d))
        cat(sprintf("Cache file %s is empty or malformed (%d rows); ignoring it.\n",
                    basename(f), nrow(d)))
      return(NULL)
    }
    d
  }
  gz  <- read_one(gz_file)
  csv <- read_one(csv_file)
  # A truncated gzip can read part-way without error, so when both copies exist
  # (WEATHER_CACHE_KEEP_CSV) take whichever is more complete and say so. They are
  # written in the same run, so any disagreement means one copy is damaged.
  if (!is.null(gz) && !is.null(csv) && nrow(csv) > nrow(gz))
    cat(sprintf("gz cache looks truncated (%d rows vs %d in the CSV); using the CSV.\n",
                nrow(gz), nrow(csv)))
  pick <- if (is.null(gz)) csv
          else if (!is.null(csv) && nrow(csv) > nrow(gz)) csv else gz
  if (is.null(pick)) {
    cat("No readable cache found; building fresh.\n")
    return(list(data = NULL, fmt = NA_character_))
  }
  used <- if (!is.null(gz) && identical(pick, gz)) gz_file else csv_file
  cat(sprintf("Cache file: %s (%.0f KB); read %d rows, %d points\n",
              basename(used), file.info(used)$size / 1024,
              nrow(pick), length(unique(pick$pid))))
  list(data = pick, fmt = if (grepl("\\.gz$", used)) "gz" else "csv")
}
cache_in <- read_cache()
cache <- cache_in$data
read_fmt <- if (is.na(cache_in$fmt)) "none" else cache_in$fmt
if (is.null(cache) || !all(c("infect","semi","wet_hours") %in% names(cache))) {
  if (!is.null(cache)) cat("Cache lacks BLASTAM columns; rebuilding from scratch.\n")
  cache <- data.table(pid = character(), lon = numeric(), lat = numeric(),
                      date = as.Date(character()), TEMP = numeric(), RHUM = numeric(),
                      RAIN = numeric(), wet_hours = numeric(), temp_wet = numeric(),
                      infect = integer(), semi = integer())
}
cache[, date := as.Date(date)]

# Drop cached points that are not on the current target grid. This matters after a
# change to GRID_RES_FINEST: points from the old grid do not coincide with the new
# one, so they would never be refreshed again, yet would still be modelled - feeding
# progressively staler weather into the map. Better to drop them than to map them.
if (nrow(cache) > 0) {
  orphans <- setdiff(unique(cache$pid), targets$pid)
  if (length(orphans) > 0) {
    cat(sprintf("Dropping %d cached point(s) not on the current %.2f deg grid (off-grid after a resolution change; they could never be refreshed).\n",
                length(orphans), GRID_RES_FINEST))
    cache <- cache[!pid %in% orphans]
  }
}

# ---- Decide what to fetch --------------------------------------------------
cost_new <- (CROP_AGE_DAYS / 7) * 0.75
cached_pids <- unique(cache$pid)
last_by <- if (nrow(cache) > 0) cache[, .(last = max(date)), by = pid] else data.table(pid = character(), last = as.Date(character()))
# Cap the whole run by NUMBER OF FETCHES so a run always finishes inside the time
# limit (cold fetches are the bottleneck, ~15-17/min). Refresh only points stale by
# >= REFRESH_MIN_STALE_DAYS, most-stale first, up to that cap; fresh/current points
# are skipped. Any leftover fetches go to adding new points, also kept under the
# daily weighted-call ceiling. Across seven daily runs this covers the grid ~weekly.
max_fetch  <- if (exists("GRID_MAX_FETCHES_PER_RUN")) GRID_MAX_FETCHES_PER_RUN else 1000L
wt_cap     <- if (exists("DAILY_WEIGHTED_CAP")) DAILY_WEIGHTED_CAP else TARGET_CALLS_PER_RUN
stale_days <- if (exists("REFRESH_MIN_STALE_DAYS")) REFRESH_MIN_STALE_DAYS else 0L
eligible <- last_by[last < end_date & last <= (end_date - stale_days)][order(last)]
n_refresh <- min(nrow(eligible), max_fetch)
maintain <- targets[pid %in% eligible$pid[seq_len(n_refresh)]]
to_add   <- targets[!pid %in% cached_pids]
maintain_cost <- nrow(maintain) * 1
# adds limited by (a) remaining fetch budget and (b) remaining weighted headroom
n_add <- min(nrow(to_add),
             max_fetch - n_refresh,
             floor(max(0, wt_cap - maintain_cost) / cost_new))
add <- if (n_add > 0) to_add[seq_len(n_add)] else to_add[0]
skipped_fresh <- length(cached_pids) - nrow(eligible)
cat(sprintf("Cache: %d points (%d fresh, skipped). Refresh %d, add %d new (%d/%d fetches, ~%.0f weighted)\n",
            length(cached_pids), skipped_fresh, nrow(maintain), nrow(add),
            nrow(maintain) + nrow(add), max_fetch, maintain_cost + nrow(add) * cost_new))

# ---- Fetch (per point, only missing dates) --------------------------------
fetch_point <- function(lon, lat, start) {
  w <- tryCatch(POINT_FETCH_FN(lat = lat, lon = lon, start_date = start, end_date = end_date),
                error = function(e) NULL)
  if (is.null(w) || nrow(w) == 0) return(NULL)
  data.table(pid = key_of(lon, lat), lon = lon, lat = lat, date = as.Date(w$date),
             TEMP = w$TEMP, RHUM = w$RHUM, RAIN = w$RAIN, wet_hours = w$wet_hours,
             temp_wet = w$temp_wet, infect = as.integer(w$infect),
             semi = as.integer(w$semi))
}
new_rows <- list()
# Points actually attempted this run. The retry below must only revisit points that
# were tried and failed, never points the time budget stopped us from reaching.
attempted_pids <- character(0)
# Parallel, rate-limited fetch. Points are fetched in chunks of GRID_CONC at once
# (mclapply forks on Linux; serial on Windows), and each chunk is held to at least
# a minimum duration so the weighted-call rate stays under GRID_TARGET_PER_MIN
# (< the 600/min limit) whether the archive responds fast or slow. This overlaps
# the slow cold hourly requests instead of doing them one at a time.
GRID_CONC <- if (exists("GRID_CONC")) GRID_CONC else 6L
GRID_TARGET_PER_MIN <- if (exists("GRID_TARGET_PER_MIN")) GRID_TARGET_PER_MIN else 450
pace_mult <- if (exists("GRID_PACE")) GRID_PACE else 1

# Coarser batch pause: once a batch's worth of weighted calls is made, wait out
# the hour so the next batch starts with a clear hourly quota. This lets one run
# do GRID_BATCHES * GRID_BATCH_CALLS calls (more points -> higher resolution)
# without breaching the 5000/hour limit.
n_batches   <- if (exists("GRID_BATCHES")) GRID_BATCHES else 1L
batch_calls <- if (exists("GRID_BATCH_CALLS")) GRID_BATCH_CALLS else TARGET_CALLS_PER_RUN
batch_wait  <- if (exists("GRID_BATCH_WAIT_S")) GRID_BATCH_WAIT_S else 3660
planned_total <- nrow(maintain) * 1 + nrow(add) * cost_new
calls_made <- 0
next_pause <- batch_calls
account_and_maybe_pause <- function(weighted) {
  calls_made <<- calls_made + weighted
  if (n_batches > 1 && calls_made >= next_pause && calls_made < planned_total && pace_mult > 0) {
    cat(sprintf("  batch complete (~%.0f weighted calls); pausing %.0f min for the hourly limit...\n",
                calls_made, batch_wait / 60))
    Sys.sleep(batch_wait)
    next_pause <<- next_pause + batch_calls
  }
}

do_fetch <- function(tab, start_fun, label, cost, deadline = fetch_deadline) {
  n <- nrow(tab); got <- 0L
  # Concurrency adapts to whatever actually delivers the most SUCCESSFUL fetches per
  # minute, which is not the same as the highest request rate. Raising concurrency
  # raises the request rate but also the failure rate, and past some point the extra
  # failures cost more than the extra throughput buys. So rather than chasing a gross
  # rate, this measures net (successful) throughput over a window of chunks and hill
  # climbs: keep moving concurrency the way that improved it, reverse when it drops.
  # The pacer below still caps the request rate, so this cannot breach API limits.
  conc     <- max(1L, as.integer(GRID_CONC))
  conc_min <- if (exists("GRID_CONC_MIN")) as.integer(GRID_CONC_MIN) else 4L
  conc_max <- if (exists("GRID_CONC_MAX")) as.integer(GRID_CONC_MAX) else 32L
  fail_lim <- if (exists("GRID_FAIL_BACKOFF")) GRID_FAIL_BACKOFF else 0.5
  probe_every <- if (exists("GRID_PROBE_CHUNKS")) as.integer(GRID_PROBE_CHUNKS) else 6L
  win_ok <- 0L; win_secs <- 0; win_chunks <- 0L; win_tried <- 0L
  last_net <- -1; step_dir <- 1L; conc_start <- conc
  tried <- 0L; got <- 0L
  i <- 1L
  t_start <- Sys.time()
  while (i <= n) {
    # Wall-clock guard: if the archive is slow, stop starting new chunks once the
    # fetch budget is spent, so the run still models, saves and commits what it has
    # instead of being killed by the job timeout (which would lose all progress).
    if (!is.na(deadline) && Sys.time() > deadline) {
      cat(sprintf("  %s: time budget reached at %d/%d; stopping this phase and moving on.\n",
                  label, i - 1L, n))
      break
    }
    j <- min(i + conc - 1L, n)
    t0 <- Sys.time()
    idx <- i:j
    chunk <- tryCatch(
      parallel::mclapply(idx, function(m)
        fetch_point(tab$lon[m], tab$lat[m], start_fun(tab$pid[m])), mc.cores = conc),
      error = function(e) lapply(idx, function(m)
        fetch_point(tab$lon[m], tab$lat[m], start_fun(tab$pid[m]))))
    ok_n <- 0L
    for (r in chunk) if (is.data.frame(r) && nrow(r) > 0) {
      new_rows[[length(new_rows) + 1L]] <<- r; ok_n <- ok_n + 1L
    }
    got <- got + ok_n; tried <- tried + length(idx)
    attempted_pids <<- c(attempted_pids, tab$pid[idx])
    if (j %% 100 < conc || j == n) cat(sprintf("  %s %d/%d (%d ok)\n", label, j, n, got))
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    win_ok <- win_ok + ok_n; win_tried <- win_tried + length(idx)
    win_secs <- win_secs + max(el, 0.01); win_chunks <- win_chunks + 1L
    if (win_chunks >= probe_every) {
      net <- win_ok / win_secs                       # successful fetches per second
      wfail <- 1 - win_ok / max(1L, win_tried)
      if (wfail > fail_lim) {
        # Most requests are failing: back off regardless of throughput, both to stop
        # hammering a struggling archive and because the work is largely wasted.
        step_dir <- -1L
      } else if (last_net >= 0 && net < last_net) {
        step_dir <- -step_dir                        # that move made things worse
      }
      step <- max(1L, as.integer(round(conc * 0.25)))
      conc <- min(conc_max, max(conc_min, conc + step_dir * step))
      last_net <- net
      win_ok <- 0L; win_tried <- 0L; win_secs <- 0; win_chunks <- 0L
    }
    min_chunk_s <- (conc * cost) / (GRID_TARGET_PER_MIN / 60)   # recomputed: conc varies
    if (pace_mult > 0 && el < min_chunk_s) Sys.sleep((min_chunk_s - el) * pace_mult)
    account_and_maybe_pause(length(idx) * cost)
    i <- j + 1L
  }
  # Gross and NET rates plus where concurrency settled. Net is the number that
  # matters: it is the cells actually added to the map per minute.
  el_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  if (tried > 0 && el_min > 0)
    cat(sprintf("  %s done: %d/%d ok (%.1f%% failed) in %.1f min; gross %.0f/min, NET %.1f/min; concurrency %d -> %d\n",
                label, got, tried, 100 * (tried - got) / tried, el_min,
                tried / el_min, got / el_min, conc_start, conc))
}
# Deadline for all fetching this run (refresh + add), measured from RUN_T0 (script
# start) so cache loading and grid building count against it, not just fetching.
# Adds are fetched after refreshes, so a slow archive spends the budget keeping
# existing cells current first, then grows the grid with whatever time is left.
# GRID_RESERVE_MINUTES is the time still needed after fetching for modelling,
# rendering, saving and committing. The workflow's timeout-minutes must exceed
# GRID_MAX_MINUTES + GRID_RESERVE_MINUTES, or the job is killed before this
# script's own graceful stop can save anything.
GRID_MAX_MINUTES <- if (exists("GRID_MAX_MINUTES")) GRID_MAX_MINUTES else 70L
reserve_min <- if (exists("GRID_RESERVE_MINUTES")) GRID_RESERVE_MINUTES else 25L
fetch_deadline <- as.POSIXct(NA)
if (!is.na(GRID_MAX_MINUTES) && GRID_MAX_MINUTES > 0)
  fetch_deadline <- RUN_T0 + GRID_MAX_MINUTES * 60
cat(sprintf("Fetch budget %.0f min from run start (%.0f min used so far); reserve %.0f min after fetching, so the workflow timeout should exceed %.0f min.\n",
            as.numeric(GRID_MAX_MINUTES), as.numeric(difftime(Sys.time(), RUN_T0, units = "mins")),
            as.numeric(reserve_min), as.numeric(GRID_MAX_MINUTES) + as.numeric(reserve_min)))
# Refresh gets its own, earlier deadline so it cannot consume the entire budget.
# Without this the grid stalls: once refreshing every cached cell fills the budget,
# no time is ever left to add new ones and the map can never grow. Reserving a share
# for adds costs a few cells on this run's map (unrefreshed cells are absent under
# the "latest" window) and buys permanent growth.
add_frac <- if (exists("GRID_ADD_RESERVE_FRAC")) GRID_ADD_RESERVE_FRAC else 0.2
refresh_deadline <- fetch_deadline
if (!is.na(fetch_deadline) && nrow(add) > 0 && add_frac > 0)
  refresh_deadline <- RUN_T0 + GRID_MAX_MINUTES * 60 * (1 - add_frac)
if (nrow(maintain) > 0) {
  last_by <- cache[, .(last = max(date)), by = pid]; setkey(last_by, pid)
  ms <- function(k) { l <- last_by[.(k), last]; if (length(l) == 0 || is.na(l)) emergence else max(emergence, l + 1) }
  do_fetch(maintain, ms, "refresh", 1, deadline = refresh_deadline)
}
if (nrow(add) > 0) do_fetch(add, function(k) emergence, "add", cost_new, deadline = fetch_deadline)

# Retry points that were ATTEMPTED and failed, in a second pass at modest
# concurrency. Both refreshes and adds are retried: with the latest-date window a
# point that misses end_date is dropped from this run's map, so a failed refresh
# costs a visible cell, not just a delay.
# Two bounds matter. Only attempted points are eligible, so points the time budget
# stopped us reaching are left for next run rather than fetched here. And the retry
# honours the same deadline, so it cannot run past the budget into the job timeout.
fetched_pids <- if (length(new_rows) > 0) unique(rbindlist(new_rows)$pid) else character(0)
tried_pids <- unique(attempted_pids)
failed_ref <- if (nrow(maintain) > 0) maintain[pid %in% tried_pids & !pid %in% fetched_pids] else maintain[0]
failed_add <- add[pid %in% tried_pids & !pid %in% fetched_pids]
retry_one <- function(tab, start_fun, label) {
  if (nrow(tab) == 0) return(invisible())
  if (!is.na(fetch_deadline) && Sys.time() > fetch_deadline) {
    cat(sprintf("  skipping retry of %d failed %s point(s): fetch budget spent.\n",
                nrow(tab), label))
    return(invisible())
  }
  cat(sprintf("  parallel retry for %d failed %s point(s)\n", nrow(tab), label))
  rc <- max(1L, min(as.integer(GRID_CONC), 8L))
  done <- 0L
  i <- 1L
  while (i <= nrow(tab)) {
    if (!is.na(fetch_deadline) && Sys.time() > fetch_deadline) {
      cat(sprintf("  %s retry stopped at %d/%d: fetch budget spent.\n", label, i - 1L, nrow(tab)))
      break
    }
    j <- min(i + rc - 1L, nrow(tab)); idx <- i:j
    res <- tryCatch(
      parallel::mclapply(idx, function(k)
        fetch_point(tab$lon[k], tab$lat[k], start_fun(tab$pid[k])), mc.cores = rc),
      error = function(e) lapply(idx, function(k)
        fetch_point(tab$lon[k], tab$lat[k], start_fun(tab$pid[k]))))
    for (r in res) if (is.data.frame(r) && nrow(r) > 0) {
      new_rows[[length(new_rows) + 1L]] <<- r; done <- done + 1L
    }
    i <- j + 1L
  }
  cat(sprintf("  %s retry recovered %d point(s)\n", label, done))
}
if (nrow(failed_ref) > 0) {
  last_by_r <- cache[, .(last = max(date)), by = pid]; setkey(last_by_r, pid)
  msr <- function(k) { l <- last_by_r[.(k), last]; if (length(l) == 0 || is.na(l)) emergence else max(emergence, l + 1) }
  retry_one(failed_ref, msr, "refresh")
}
retry_one(failed_add, function(k) emergence, "add")
if (length(new_rows) > 0)
  cache <- rbindlist(list(cache, rbindlist(new_rows)), use.names = TRUE)

cache <- cache[date >= emergence & date <= end_date]
setorder(cache, pid, date)
cache <- unique(cache, by = c("pid", "date"), fromLast = TRUE)

# ---- Model both per point --------------------------------------------------
run_tag <- format(Sys.Date(), "%Y-%m-%d")
writeLines(run_tag, file.path(OUT, "run_date.txt"))  # so the email attaches the exact files
model_pt <- function(dt) {
  epi <- NA_real_
  if (nrow(dt) >= 20) {
    w <- data.table(YYYYMMDD = dt$date, DOY = as.integer(format(dt$date, "%j")),
                    TEMP = dt$TEMP, RHUM = dt$RHUM, RAIN = dt$RAIN, LAT = dt$lat[1], LON = dt$lon[1])
    setorder(w, YYYYMMDD)
    lb <- tryCatch(predict_leaf_blast(w, emergence = emergence,
                                      duration = as.integer(min(120L, nrow(w)))),
                   error = function(e) NULL)
    if (!is.null(lb) && nrow(lb) > 0) epi <- lb$intensity[nrow(lb)]
  }
  win <- if (exists("BLASTAM_WINDOW_DAYS")) BLASTAM_WINDOW_DAYS else 21L
  list(intensity = epi,
       events = sum(dt$infect[dt$date > (model_end - win)], na.rm = TRUE))
}
# COMMON WINDOW, REAL WEATHER ONLY. With a rolling refresh, cells differ in how
# current they are, so the map's window ends at the newest date that essentially
# every cell has actually reached, and each cell's series is TRUNCATED to it. That
# gives one shared calendar window using only real weather: nothing is padded or
# extrapolated. The map therefore lags the newest data by about the refresh cycle,
# which is the price of covering more cells than one run can fetch.
# GRID_WINDOW_COVERAGE is the fraction of cells that must reach the window end;
# setting it below 1 lets a few stragglers be dropped so the window stays fresher.
pt_end <- cache[, .(mx = max(date)), by = pid]
n_cache_pts <- nrow(pt_end)
wmode <- if (exists("GRID_WINDOW_MODE")) GRID_WINDOW_MODE else "coverage"
if (identical(wmode, "latest")) {
  # LATEST: the window always ends at the archive edge. Cells that did not reach it
  # this run (fetch failure, or the time budget stopping the run) are simply absent
  # from the map. Freshness is fixed; the cell count is what varies.
  model_end <- end_date
} else {
  # COVERAGE: pull the window back to the newest date GRID_WINDOW_COVERAGE of cells
  # reached, so nearly every cell is mapped at the cost of a staler window.
  cover <- if (exists("GRID_WINDOW_COVERAGE")) GRID_WINDOW_COVERAGE else 1
  mx_sorted <- sort(pt_end$mx)
  idx <- max(1L, ceiling((1 - cover) * length(mx_sorted)))
  model_end <- min(mx_sorted[idx], end_date)
}
current_pids <- pt_end[mx >= model_end, pid]
held_out <- n_cache_pts - length(current_pids)
model_cache <- cache[pid %in% current_pids & date <= model_end]
cat(sprintf("Window [%s] ends %s: %d of %d cells reach it, %d absent from this map (%d days behind the archive edge).\n",
            wmode, format(model_end), length(current_pids), n_cache_pts, held_out,
            as.integer(end_date - model_end)))

pm <- model_cache[, { m <- model_pt(.SD); .(lon = lon[1], lat = lat[1],
              intensity = m$intensity, events = m$events) }, by = pid]
pm_epi <- pm[!is.na(intensity)]
map_spacing <- if (nrow(pm) > 0) sqrt(prod(ext[c(2,4)] - ext[c(1,3)]) / nrow(pm)) else NA
cat(sprintf("Modelled %d points (all over the same window to %s); effective spacing ~%.2f deg\n",
            nrow(pm), format(model_end), map_spacing))

# ---- Render helper (IDW surface -> PNG + GeoTIFF) -------------------------
load_bundled <- function(fname) {
  f <- file.path(SCRIPT_DIR, fname)
  v <- tryCatch(if (file.exists(f)) terra::vect(f) else NULL, error = function(e) NULL)
  if (is.null(v) || nrow(v) == 0) return(NULL)
  tryCatch(terra::crop(v, terra::ext(ext)), error = function(e) v)
}
rivers <- if (isTRUE(SHOW_RIVERS)) load_bundled("australia_rivers.geojson") else NULL
roads  <- if (isTRUE(SHOW_ROADS))  load_bundled("australia_roads.geojson") else NULL

render_map <- function(pts, valcol, title, colours, hmax, legend, base) {
  if (nrow(pts) < 3) { cat("Too few points to render", base, "\n"); return(invisible()) }
  r0 <- terra::rast(xmin = ext[1], xmax = ext[2], ymin = ext[3], ymax = ext[4],
                    resolution = GRID_RES_FINEST / 2, crs = "EPSG:4326")
  v <- terra::vect(as.matrix(pts[, .(lon, lat)]), type = "points", crs = "EPSG:4326")
  v$val <- pts[[valcol]]
  r <- tryCatch(terra::interpIDW(r0, v, field = "val", radius = 6, power = 2, maxPoints = 12),
                error = function(e) terra::rasterize(v, r0, field = "val", fun = "mean"))
  names(r) <- valcol
  if (!is.null(land_poly)) r <- terra::mask(r, land_poly)
  if (isTRUE(WRITE_GEOTIFF))
    terra::writeRaster(r, file.path(OUT, sprintf("%s_%s.tif", base, run_tag)), overwrite = TRUE)
  rmax <- if (!is.null(hmax)) hmax else { m <- terra::global(r, "max", na.rm = TRUE)[1,1]; if (!is.finite(m) || m <= 0) 1 else m }
  png_file <- file.path(OUT, sprintf("%s_%s.png", base, run_tag))
  png(png_file, width = 1000, height = 900, res = 120)
  op <- par(mar = c(4, 4, 3, 5))
  ramp <- grDevices::colorRampPalette(colours)(100)
  terra::plot(r, col = ramp, range = c(0, rmax), xlab = "Longitude", ylab = "Latitude",
              main = sprintf("%s  %s", title, run_tag), plg = list(title = legend))
  if (isTRUE(SHOW_RIVERS) && !is.null(rivers)) try(terra::lines(rivers, col = COL_RIVER, lwd = 0.6), silent = TRUE)
  if (isTRUE(SHOW_ROADS)  && !is.null(roads))  try(terra::lines(roads,  col = COL_ROAD,  lwd = 0.5), silent = TRUE)
  if (isTRUE(SHOW_COAST)  && !is.null(land_poly)) try(terra::lines(land_poly, col = COL_COAST, lwd = 1), silent = TRUE)
  if (isTRUE(SHOW_TOWNS) && exists("MONITOR_TOWNS") && nrow(MONITOR_TOWNS) > 0) {
    points(MONITOR_TOWNS$lon, MONITOR_TOWNS$lat, pch = 21, bg = "white", col = COL_TOWN, cex = 1.1, lwd = 1.4)
    text(MONITOR_TOWNS$lon, MONITOR_TOWNS$lat, MONITOR_TOWNS$name, pos = 4, offset = 0.3, cex = 0.5, col = COL_TOWN)
  }
  par(op); dev.off()
  file.copy(png_file, file.path(OUT, sprintf("%s_latest.png", base)), overwrite = TRUE)
  cat("Heatmap:", png_file, "\n")
}

# EPIRICE intensity (%) and BLASTAM infection days
pm_epi[, intensity_pct := intensity * 100]
render_map(pm_epi, "intensity_pct", "EPIRICE potential risk (%)", HEAT_COLOURS, HEAT_MAX,
           "intensity %", "epirice_heatmap")
render_map(pm, "events", sprintf("BLASTAM infection days (last %d)", BLASTAM_WINDOW_DAYS),
           BLASTAM_HEAT_COLOURS,
           BLASTAM_HEAT_MAX, "days", "blastam_heatmap")

# ---- Save cache -----------------------------------------------------------
# Round weather values before saving to keep the cache compact (no meaningful
# model impact). lon/lat/pid/date/flags are left exact.
csvdt <- copy(cache)
csvdt[, `:=`(TEMP = round(TEMP, 1), RHUM = round(RHUM, 0), RAIN = round(RAIN, 1),
             temp_wet = round(temp_wet, 1), wet_hours = round(wet_hours, 0),
             lon = round(lon, 4), lat = round(lat, 4))]

# Primary write is gz (small); verify it reads back, else fall back to plain CSV.
write_cache <- function(dt) {
  ok_gz <- tryCatch({
    fwrite(dt, gz_file)
    nrow(read_gz_dt(gz_file)) == nrow(dt)
  }, error = function(e) { cat("gz verify error:", conditionMessage(e), "\n"); FALSE })
  if (ok_gz) {
    if (WEATHER_CACHE_KEEP_CSV) { fwrite(dt, csv_file); fmt <- "gz+csv" }
    else { if (file.exists(csv_file)) file.remove(csv_file); fmt <- "gz" }
    return(list(fmt = fmt, kb = file.info(gz_file)$size / 1024))
  }
  cat("gz cache write/verify failed; falling back to plain CSV.\n")
  fwrite(dt, csv_file)
  if (file.exists(gz_file)) file.remove(gz_file)
  list(fmt = "csv", kb = file.info(csv_file)$size / 1024)
}
wc <- write_cache(csvdt)
cat(sprintf("Cache saved: %d points as %s (%.0f KB)\nDone.\n",
            length(unique(cache$pid)), wc$fmt, wc$kb))

# Stats line for the email: mapped points, spacing, last run's points, target
# spacing, stored format and size, format actually read, and the common window end.
# The count is the number MAPPED (nrow(pm)), matching map_spacing, not the whole
# cache, since cells short of the common window are not on this run's map.
writeLines(sprintf("%d|%.2f|%d|%.2f|%s|%.0f|%s|%s",
                   nrow(pm), map_spacing,
                   length(cached_pids), GRID_RES_FINEST, wc$fmt, wc$kb, read_fmt,
                   format(model_end)),
           file.path(OUT, "map_stats.txt"))

# On the silent midweek fetch-only run, record a status line (committed to the
# repo) so the next full run can confirm in the email that it ran and went well.
if (Sys.getenv("BLAST_MIDWEEK") == "1") {
  added <- length(unique(cache$pid)) - length(cached_pids)
  writeLines(sprintf("%s|%d|%d|%d|%s|%.0f",
                     format(Sys.Date()), length(unique(cache$pid)),
                     length(cached_pids), added, wc$fmt, wc$kb),
             file.path(OUT, "midweek_status.txt"))
  cat(sprintf("Midweek fetch-only run: +%d points, cache now %d.\n",
              added, length(unique(cache$pid))))
}

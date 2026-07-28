#!/usr/bin/env Rscript
################################################################################
# run_blast_grid.R
#
# Two continental heatmaps from a shared, cached weather record with dynamic
# resolution:
#   1) EPIRICE leaf blast intensity (mechanistic epidemic, % leaf area diseased)
#   2) BLASTAM infection days (count of days favourable for new infection)
#
# Both derive from ONE hourly Open-Meteo fetch per point: hourly data gives
# BLASTAM its leaf wetness hours and is aggregated to daily values for EPIRICE.
# The cache holds daily rows; each run refreshes existing points and adds new
# ones, so the map sharpens over successive runs.
#
# CHANGES FROM THE CONCURRENT VERSION
#   * Fetching is batched (openmeteo_batch.R): many locations per request, real
#     HTTP status codes, token bucket pacing on weighted calls. The adaptive
#     concurrency controller is gone; its pacer made throughput algebraically
#     independent of concurrency, so it hill climbed on noise while the failure
#     rate rose.
#   * cost_new now uses the API's actual weight formula. The old
#     (CROP_AGE_DAYS / 7) * 0.75 overestimated by about 47%, so runs stopped
#     roughly 1.5x short of the budget they were allowed to spend.
#   * Refreshes request a fixed REFRESH_TAIL_DAYS tail rather than only the
#     missing days. Same weighted cost because of the API's 14 day floor, and it
#     self heals gaps and absorbs ERA5T revisions.
#   * The retry pass has its own reserved wall clock. Previously the add phase
#     ran to the full deadline, so the retry was unreachable by construction.
#   * Effective spacing is computed over LAND, not the bounding box, and the
#     completed lattice level is reported alongside it.
#   * The IDW search radius follows the achieved spacing instead of a fixed 6
#     degrees (about 660 km), and uncovered cells are masked.
#
# The test hook changed: POINT_FETCH_FN(lat, lon, start, end) is replaced by
# GRID_ON_POINT(pid, lon, lat, hourly_dt) -> cache rows.
################################################################################

Sys.setenv(TZ = "Australia/Sydney")
RUN_T0 <- Sys.time()

SCRIPT_DIR <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile), winslash = "/"),
  error = function(e) normalizePath(getwd(), winslash = "/"))

source(file.path(SCRIPT_DIR, "blast_config.R"))
source(file.path(SCRIPT_DIR, "epirice_model.R"))
source(file.path(SCRIPT_DIR, "blastam_model.R"))
source(file.path(SCRIPT_DIR, "openmeteo_wth.R"))
source(file.path(SCRIPT_DIR, "openmeteo_batch.R"))

suppressPackageStartupMessages({library(data.table); library(terra)})

ext <- GRID_EXTENT
key_of <- function(lon, lat) paste0(sprintf("%.4f", lon), "_", sprintf("%.4f", lat))
CACHE_COLS <- c("pid","lon","lat","date","TEMP","RHUM","RAIN","wet_hours","temp_wet","infect","semi")

# Hourly JSON for one point -> the daily cache rows both models read.
if (!exists("GRID_ON_POINT")) GRID_ON_POINT <- function(pid, lon, lat, hw) {
  w <- tryCatch(blastam_daily_from_hourly(hw), error = function(e) NULL)
  if (is.null(w) || nrow(w) == 0) return(NULL)
  data.table(pid = pid, lon = lon, lat = lat, date = as.Date(w$date),
             TEMP = w$TEMP, RHUM = w$RHUM, RAIN = w$RAIN,
             wet_hours = w$wet_hours, temp_wet = w$temp_wet,
             infect = as.integer(w$infect), semi = as.integer(w$semi))
}

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

# ---- Target lattice, coarse to fine ---------------------------------------
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
  g[, lvl := lv]
  setorder(g, lvl, lat, lon)
  g[, pid := key_of(lon, lat)]
  g[]
}
targets <- build_targets()
LEVEL_RES <- c(GRID_RES_LEVELS, GRID_RES_FINEST)[seq_len(max(targets$lvl))]

# ---- Window ----------------------------------------------------------------
end_date  <- Sys.Date() - ARCHIVE_LAG_DAYS
emergence <- end_date - CROP_AGE_DAYS
win_dates <- seq(emergence, end_date, by = "day")
cat(sprintf("Window %s to %s (%d days); target %d land points at %.2f deg\n",
            emergence, end_date, length(win_dates), nrow(targets), GRID_RES_FINEST))

# ---- Load cache ------------------------------------------------------------
OUT <- file.path(SCRIPT_DIR, OUTPUT_DIR)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
gz_file  <- file.path(OUT, WEATHER_CACHE_GZ)
csv_file <- file.path(OUT, WEATHER_CACHE_CSV)

# fread() cannot read .gz unless R.utils is installed (it is not in the runner
# container), so gz is decompressed with a base R gzfile connection first. That
# keeps the cache dependency free.
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
  read_one <- function(f) {
    if (!file.exists(f)) return(NULL)
    is_gz <- grepl("\\.gz$", f)
    d <- tryCatch(if (is_gz) read_gz_dt(f) else fread(f, colClasses = list(character = "pid")),
                  error = function(e) {
                    cat(sprintf("Cache read FAILED (%s): %s\n", basename(f), conditionMessage(e))); NULL })
    # A damaged file can read without error yet be empty or missing columns (a
    # non-gzip file opened through gzfile reads as plain text, for instance).
    if (is.null(d) || nrow(d) == 0 || !all(CACHE_COLS %in% names(d))) {
      if (!is.null(d))
        cat(sprintf("Cache file %s is empty or malformed (%d rows); ignoring it.\n",
                    basename(f), nrow(d)))
      return(NULL)
    }
    d
  }
  gz  <- read_one(gz_file)
  csv <- read_one(csv_file)
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

# Drop cached points that are not on the current target grid. After a change to
# GRID_RES_FINEST the old points do not coincide with the new lattice, so they
# would never be refreshed again yet would still be modelled, feeding
# progressively staler weather into the map.
if (nrow(cache) > 0) {
  orphans <- setdiff(unique(cache$pid), targets$pid)
  if (length(orphans) > 0) {
    cat(sprintf("Dropping %d cached point(s) not on the current %.2f deg grid.\n",
                length(orphans), GRID_RES_FINEST))
    cache <- cache[!pid %in% orphans]
  }
}

# ---- Failure ledger --------------------------------------------------------
# Points whose fetch failed previously. Without this the planner re-selects the
# same points in the same order every run and fails the same way.
ledger_file <- file.path(OUT, if (exists("FAIL_LEDGER_FILE")) FAIL_LEDGER_FILE else "fetch_failures.csv")
max_strikes <- if (exists("FAIL_LEDGER_MAX_STRIKES")) FAIL_LEDGER_MAX_STRIKES else 4L
fails <- if (file.exists(ledger_file))
  tryCatch(fread(ledger_file, colClasses = list(character = "pid")), error = function(e) NULL) else NULL
if (is.null(fails) || !all(c("pid","strikes","last_try","last_status") %in% names(fails)))
  fails <- data.table(pid = character(), strikes = integer(),
                      last_try = as.Date(character()), last_status = character())
fails[, last_try := as.Date(last_try)]
# Exponential cooloff: a point with s strikes is skipped for 2^s days.
benched <- fails[strikes >= max_strikes & last_try > (Sys.Date() - 2^pmin(strikes, 8L)), pid]
if (length(benched) > 0)
  cat(sprintf("%d point(s) benched after repeated fetch failures; they will be retried later.\n",
              length(benched)))

# ---- Decide what to fetch --------------------------------------------------
# Weighted cost. The API charges max(1, nvars/10) * max(1, ndays/14) per
# location, with a 14 day FLOOR. That floor is why a refresh costs 1 regardless
# of how few days it needs, and why refreshes take a full 14 day tail below.
n_days_add <- as.integer(end_date - emergence) + 1L
cost_new   <- om_weight_per_location(n_days_add, length(OM_HOURLY_VARS))
tail_start <- max(emergence, end_date - (REFRESH_TAIL_DAYS - 1L))
cost_ref   <- om_weight_per_location(as.integer(end_date - tail_start) + 1L,
                                     length(OM_HOURLY_VARS))

cached_pids <- unique(cache$pid)
last_by <- if (nrow(cache) > 0) cache[, .(last = max(date)), by = pid] else
  data.table(pid = character(), last = as.Date(character()))

max_fetch  <- GRID_MAX_FETCHES_PER_RUN
wt_cap     <- DAILY_WEIGHTED_CAP
stale_days <- REFRESH_MIN_STALE_DAYS

eligible <- last_by[last < end_date & last <= (end_date - stale_days)][order(last)]
n_refresh <- min(nrow(eligible), max_fetch)
maintain  <- targets[pid %in% eligible$pid[seq_len(n_refresh)]]
maintain_cost <- nrow(maintain) * cost_ref

to_add <- targets[!pid %in% cached_pids & !pid %in% benched]
n_add <- min(nrow(to_add),
             max_fetch - nrow(maintain),
             floor(max(0, wt_cap - maintain_cost) / cost_new))
add <- if (n_add > 0) to_add[seq_len(n_add)] else to_add[0]

skipped_fresh <- length(cached_pids) - nrow(eligible)
cat(sprintf("Cache: %d points (%d fresh, skipped). Refresh %d @ %.2f, add %d new @ %.2f (%d/%d fetches, ~%.0f of %.0f weighted)\n",
            length(cached_pids), skipped_fresh, nrow(maintain), cost_ref,
            nrow(add), cost_new, nrow(maintain) + nrow(add), max_fetch,
            maintain_cost + nrow(add) * cost_new, wt_cap))

# ---- Deadlines -------------------------------------------------------------
# Measured from RUN_T0 so cache loading and grid building count against the
# budget, not just fetching.
reserve_min   <- GRID_RESERVE_MINUTES
retry_reserve <- if (exists("GRID_RETRY_RESERVE_MINUTES")) GRID_RETRY_RESERVE_MINUTES else 20L
# CLAMP. The retry slice is subtracted from the fetch budget, so if it is larger
# than the budget the add deadline lands in the PAST and the run fetches nothing
# while reporting "deadline". That is exactly what a short dry run would hit.
if (!is.na(GRID_MAX_MINUTES) && GRID_MAX_MINUTES > 0)
  retry_reserve <- min(as.numeric(retry_reserve), as.numeric(GRID_MAX_MINUTES) * 0.25)

fetch_deadline <- if (!is.na(GRID_MAX_MINUTES) && GRID_MAX_MINUTES > 0)
  RUN_T0 + GRID_MAX_MINUTES * 60 else as.POSIXct(NA)
add_deadline <- if (is.na(fetch_deadline)) fetch_deadline else fetch_deadline - retry_reserve * 60
add_frac <- GRID_ADD_RESERVE_FRAC
# Refresh gets a share of the ADD deadline, not of the whole budget, so the three
# phases are always ordered refresh < add < retry whatever the budget is.
refresh_deadline <- if (is.na(add_deadline)) add_deadline else {
  if (nrow(add) > 0 && add_frac > 0)
    RUN_T0 + as.numeric(difftime(add_deadline, RUN_T0, units = "secs")) * (1 - add_frac)
  else add_deadline
}
stopifnot(is.na(add_deadline) || add_deadline > RUN_T0)

cat(sprintf("Fetch budget %.0f min from run start (%.1f used); %.1f min held for retries, %.0f min reserved after fetching. Workflow timeout must exceed %.0f min.\n",
            as.numeric(GRID_MAX_MINUTES),
            as.numeric(difftime(Sys.time(), RUN_T0, units = "mins")),
            as.numeric(retry_reserve), as.numeric(reserve_min),
            as.numeric(GRID_MAX_MINUTES) + as.numeric(reserve_min)))

# ---- Fetch -----------------------------------------------------------------
new_rows <- list()
all_ledger <- list()
spent <- 0
quota_hit <- FALSE

run_phase <- function(tab, start_date, deadline, label) {
  if (nrow(tab) == 0L || isTRUE(quota_hit)) return(invisible())
  r <- fetch_points_batched(tab[, .(pid, lon, lat)], start_date, end_date,
                            on_point = GRID_ON_POINT,
                            deadline = deadline,
                            budget = max(0, wt_cap - spent),
                            label = label)
  if (length(r$rows) > 0) new_rows <<- c(new_rows, r$rows)
  if (nrow(r$ledger) > 0) all_ledger[[length(all_ledger) + 1L]] <<- r$ledger
  spent <<- spent + r$spent
  if (identical(r$stopped, "quota")) quota_hit <<- TRUE
  invisible()
}

# Refresh first: keeping existing cells current matters more than growing, since
# under "latest" window mode a stale cell drops off the map entirely.
run_phase(maintain, tail_start, refresh_deadline, "refresh")
run_phase(add,      emergence,  add_deadline,     "add")

# ---- Retry -----------------------------------------------------------------
# Only points that were ATTEMPTED and failed with a transport or HTTP error.
# "empty" means the API answered and had nothing, so retrying is pointless, and
# points the budget never reached are left for next run. This pass now has its
# own reserved wall clock; previously the add phase consumed the whole budget so
# the retry was skipped every time.
led <- if (length(all_ledger)) rbindlist(all_ledger) else
  data.table(pid = character(), status = character(), code = integer())
retryable <- led[status %in% c("http", "transport"), unique(pid)]
# If most of the run failed, the archive is having a bad day and retrying
# everything simply doubles the wasted quota and wall clock. Retry only when the
# failures look transient.
attempted <- nrow(led)
ok_frac <- if (attempted > 0) nrow(led[status == "ok"]) / attempted else 1
retry_min_ok <- if (exists("GRID_RETRY_MIN_OK_FRAC")) GRID_RETRY_MIN_OK_FRAC else 0.25
if (length(retryable) > 0 && !quota_hit && ok_frac < retry_min_ok) {
  cat(sprintf("Skipping the retry pass: only %.0f%% of %d attempted points succeeded, which looks like a systemic failure rather than transient flakiness.\n",
              100 * ok_frac, attempted))
} else if (length(retryable) > 0 && !quota_hit) {
  rt_ref <- maintain[pid %in% retryable]
  rt_add <- add[pid %in% retryable]
  cat(sprintf("Retrying %d refresh and %d add point(s) that failed with a transport or HTTP error.\n",
              nrow(rt_ref), nrow(rt_add)))
  run_phase(rt_ref, tail_start, fetch_deadline, "refresh-retry")
  run_phase(rt_add, emergence,  fetch_deadline, "add-retry")
}

# ---- Update the failure ledger ---------------------------------------------
led <- if (length(all_ledger)) rbindlist(all_ledger) else led
if (nrow(led) > 0) {
  final <- led[, .(status = last(status)), by = pid]
  ok_now <- final[status == "ok", pid]
  bad_now <- final[status %in% c("http", "transport", "empty"), ]
  fails <- fails[!pid %in% ok_now]
  if (nrow(bad_now) > 0) {
    upd <- merge(bad_now, fails[, .(pid, strikes)], by = "pid", all.x = TRUE)
    upd[is.na(strikes), strikes := 0L]
    upd[, `:=`(strikes = strikes + 1L, last_try = Sys.Date(), last_status = status)]
    fails <- rbind(fails[!pid %in% upd$pid],
                   upd[, .(pid, strikes, last_try, last_status)], fill = TRUE)
  }
  fwrite(fails, ledger_file)
  cat(sprintf("Failure ledger: %d point(s) carrying at least one strike.\n", nrow(fails)))
}

if (length(new_rows) > 0)
  cache <- rbindlist(list(cache, rbindlist(new_rows)), use.names = TRUE)

# Retain history rather than pruning to the modelling window. Weather already
# fetched costs nothing to keep, and discarding it means no past map can be
# reproduced and any retrospective analysis needs the whole grid fetched again.
keep_hist <- if (exists("CACHE_KEEP_HISTORY")) isTRUE(CACHE_KEEP_HISTORY) else FALSE
hist_days <- if (exists("CACHE_HISTORY_DAYS")) CACHE_HISTORY_DAYS else 400L
cache <- if (keep_hist) {
  cache[date >= (end_date - hist_days) & date <= end_date]
} else {
  cache[date >= emergence & date <= end_date]
}
setorder(cache, pid, date)
cache <- unique(cache, by = c("pid", "date"), fromLast = TRUE)

# ---- Model both per point --------------------------------------------------
run_tag <- format(Sys.Date(), "%Y-%m-%d")
writeLines(run_tag, file.path(OUT, "run_date.txt"))

model_pt <- function(dt) {
  epi <- NA_real_
  if (nrow(dt) >= 20) {
    w <- data.table(YYYYMMDD = dt$date, DOY = as.integer(format(dt$date, "%j")),
                    TEMP = dt$TEMP, RHUM = dt$RHUM, RAIN = dt$RAIN,
                    LAT = dt$lat[1], LON = dt$lon[1])
    setorder(w, YYYYMMDD)
    lb <- tryCatch(predict_leaf_blast(w, emergence = emergence,
                                      duration = as.integer(min(120L, nrow(w)))),
                   error = function(e) NULL)
    if (!is.null(lb) && nrow(lb) > 0) epi <- lb$intensity[nrow(lb)]
  }
  win <- BLASTAM_WINDOW_DAYS
  list(intensity = epi,
       events = sum(dt$infect[dt$date > (model_end - win)], na.rm = TRUE))
}

# COMMON WINDOW, REAL WEATHER ONLY. Cells differ in how current they are, so the
# map's window ends at a date essentially every cell has reached, and each cell's
# series is TRUNCATED to it. One shared calendar window, nothing padded or
# extrapolated.
pt_end <- if (nrow(cache) > 0) cache[, .(mx = max(date)), by = pid] else
  data.table(pid = character(), mx = as.Date(character()))
n_cache_pts <- nrow(pt_end)
wmode <- GRID_WINDOW_MODE
if (identical(wmode, "latest")) {
  model_end <- end_date
} else {
  cover <- GRID_WINDOW_COVERAGE
  mx_sorted <- sort(pt_end$mx)
  idx <- max(1L, ceiling((1 - cover) * length(mx_sorted)))
  model_end <- min(mx_sorted[idx], end_date)
}
current_pids <- pt_end[mx >= model_end, pid]
held_out <- n_cache_pts - length(current_pids)
model_cache <- cache[pid %in% current_pids &
                     date > (model_end - CROP_AGE_DAYS) & date <= model_end]
cat(sprintf("Window [%s] ends %s: %d of %d cells reach it, %d absent (%d days behind the archive edge).\n",
            wmode, format(model_end), length(current_pids), n_cache_pts, held_out,
            as.integer(end_date - model_end)))

pm <- model_cache[, { m <- model_pt(.SD); .(lon = lon[1], lat = lat[1],
              intensity = m$intensity, events = m$events) }, by = pid]
pm_epi <- pm[!is.na(intensity)]

# ---- Coverage reporting ----------------------------------------------------
# Mean spacing over LAND, from the target lattice. The old statistic divided the
# BOUNDING BOX (42 x 34 = 1428 sq deg) by the point count, counting the Southern
# Ocean as unsampled land, which understated the achieved resolution by about 40%.
land_area_deg2 <- nrow(targets) * GRID_RES_FINEST^2
map_spacing <- if (nrow(pm) > 0) sqrt(land_area_deg2 / nrow(pm)) else NA_real_

# Mean spacing is misleading while the lattice fills coarse to fine, because
# coverage is a complete coarse level plus a partial finer one. The completed
# level is the defensible claim.
lvl_target <- targets[, .(N_target = .N), by = lvl]
lvl_done   <- targets[pid %in% pm$pid, .(N_done = .N), by = lvl]
lvl_tab <- merge(lvl_target, lvl_done, by = "lvl", all.x = TRUE)[order(lvl)]
lvl_tab[is.na(N_done), N_done := 0L]
lvl_tab[, frac := N_done / N_target]
complete_lvl <- suppressWarnings(max(lvl_tab[frac >= 0.995, lvl]))
res_complete <- if (is.finite(complete_lvl)) LEVEL_RES[complete_lvl] else NA_real_
cat(sprintf("Level coverage: %s\n",
            paste(sprintf("L%d(%.2f deg) %d/%d", lvl_tab$lvl, LEVEL_RES[lvl_tab$lvl],
                          lvl_tab$N_done, lvl_tab$N_target), collapse = "; ")))
cat(sprintf("Modelled %d points to %s; %s; mean land spacing ~%s deg\n",
            nrow(pm), format(model_end),
            if (is.na(res_complete)) "no lattice level complete yet"
            else sprintf("complete to %.2f deg", res_complete),
            if (is.na(map_spacing)) "NA" else sprintf("%.2f", map_spacing)))

# ---- Render ----------------------------------------------------------------
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

  # Search radius tied to the achieved spacing. The old fixed radius of 6 was six
  # DEGREES, roughly 660 km and about fifteen times the point spacing, so it
  # interpolated smoothly straight across unsampled regions and a gap in the
  # lattice was indistinguishable from data.
  idw_rad <- min(if (exists("IDW_RADIUS_MAX")) IDW_RADIUS_MAX else 1.0,
                 max(GRID_RES_FINEST, IDW_RADIUS_MULT * map_spacing))
  r <- tryCatch(terra::interpIDW(r0, v, field = "val", radius = idw_rad,
                                 power = IDW_POWER, maxPoints = IDW_MAX_POINTS),
                error = function(e) terra::rasterize(v, r0, field = "val", fun = "mean"))
  names(r) <- valcol
  if (!is.null(land_poly)) r <- terra::mask(r, land_poly)
  if (isTRUE(MASK_UNCOVERED)) {
    # Blank anywhere with no modelled point within the search radius, so gaps
    # read as gaps. distance() on a lonlat grid returns metres.
    dmask <- tryCatch(terra::distance(r0, v), error = function(e) NULL)
    if (!is.null(dmask)) r[dmask > idw_rad * 111320] <- NA
  }

  if (isTRUE(WRITE_GEOTIFF))
    terra::writeRaster(r, file.path(OUT, sprintf("%s_%s.tif", base, run_tag)),
                       overwrite = TRUE)
  rmax <- if (!is.null(hmax)) hmax else {
    m <- terra::global(r, "max", na.rm = TRUE)[1,1]; if (!is.finite(m) || m <= 0) 1 else m }

  # Filenames stay on the RUN date so send_email.py keeps finding them, but the
  # title and footer carry the DATA date. The two differ by ARCHIVE_LAG_DAYS plus
  # any refresh lag, and a run-dated title implied a currency the map lacks.
  # To move the filename onto the data date, swap run_tag for format(model_end)
  # here and in the .tif above, and update send_email.py to match.
  png_file <- file.path(OUT, sprintf("%s_%s.png", base, run_tag))
  png(png_file, width = 1000, height = 900, res = 120)
  op <- par(mar = c(4, 4, 3, 5))
  ramp <- grDevices::colorRampPalette(colours)(100)
  terra::plot(r, col = ramp, range = c(0, rmax), xlab = "Longitude", ylab = "Latitude",
              main = sprintf("%s  weather to %s", title, format(model_end)),
              plg = list(title = legend))
  mtext(sprintf("run %s | %d cells | complete to %s deg | driver %s",
                run_tag, nrow(pts),
                if (is.na(res_complete)) "no level" else sprintf("%.2f", res_complete),
                OPENMETEO_MODEL),
        side = 1, line = 2.6, cex = 0.55, col = "#6b7378")
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

pm_epi[, intensity_pct := intensity * 100]
render_map(pm_epi, "intensity_pct", "EPIRICE potential risk (%)", HEAT_COLOURS, HEAT_MAX,
           "intensity %", "epirice_heatmap")
render_map(pm, "events", sprintf("BLASTAM infection days (last %d)", BLASTAM_WINDOW_DAYS),
           BLASTAM_HEAT_COLOURS, BLASTAM_HEAT_MAX, "days", "blastam_heatmap")

# ---- Save cache ------------------------------------------------------------
csvdt <- copy(cache)
csvdt[, `:=`(TEMP = round(TEMP, 1), RHUM = round(RHUM, 0), RAIN = round(RAIN, 1),
             temp_wet = round(temp_wet, 1), wet_hours = round(wet_hours, 0),
             lon = round(lon, 4), lat = round(lat, 4))]

# Write to a temp path, verify, then rename. A job cancelled mid write used to
# leave a truncated cache that the next run had to detect and discard.
write_cache <- function(dt) {
  tmp_gz <- paste0(gz_file, ".tmp")
  ok_gz <- tryCatch({
    fwrite(dt, tmp_gz)
    n <- nrow(read_gz_dt(tmp_gz))
    n == nrow(dt)
  }, error = function(e) { cat("gz verify error:", conditionMessage(e), "\n"); FALSE })
  if (ok_gz) {
    file.rename(tmp_gz, gz_file)
    if (WEATHER_CACHE_KEEP_CSV) {
      tmp_csv <- paste0(csv_file, ".tmp"); fwrite(dt, tmp_csv); file.rename(tmp_csv, csv_file)
      fmt <- "gz+csv"
    } else {
      if (file.exists(csv_file)) file.remove(csv_file); fmt <- "gz"
    }
    return(list(fmt = fmt, kb = file.info(gz_file)$size / 1024))
  }
  unlink(tmp_gz)
  cat("gz cache write/verify failed; falling back to plain CSV.\n")
  tmp_csv <- paste0(csv_file, ".tmp"); fwrite(dt, tmp_csv); file.rename(tmp_csv, csv_file)
  if (file.exists(gz_file)) file.remove(gz_file)
  list(fmt = "csv", kb = file.info(csv_file)$size / 1024)
}
wc <- write_cache(csvdt)
cat(sprintf("Cache saved: %d points, %d rows as %s (%.0f KB)\n",
            length(unique(cache$pid)), nrow(cache), wc$fmt, wc$kb))

# ---- Stats line for the email ----------------------------------------------
# Field 3 is now the number MAPPED on the previous run, read back before this
# run overwrites the file. It used to be the cache size at the start of this run,
# which is a different quantity, so "up from N" compared unlike things.
stats_file <- file.path(OUT, "map_stats.txt")
prev_mapped <- if (!file.exists(stats_file)) NA_integer_ else tryCatch({
  s <- strsplit(readLines(stats_file, warn = FALSE)[1], "\\|")[[1]]
  as.integer(s[1])
}, error = function(e) NA_integer_, warning = function(w) NA_integer_)

writeLines(sprintf("%d|%.2f|%d|%.2f|%s|%.0f|%s|%s|%s|%.0f",
                   nrow(pm), map_spacing,
                   if (is.na(prev_mapped)) 0L else prev_mapped,
                   GRID_RES_FINEST, wc$fmt, wc$kb, read_fmt,
                   format(model_end),
                   if (is.na(res_complete)) "" else sprintf("%.2f", res_complete),
                   spent),
           stats_file)

if (Sys.getenv("BLAST_MIDWEEK") == "1") {
  added <- length(unique(cache$pid)) - length(cached_pids)
  writeLines(sprintf("%s|%d|%d|%d|%s|%.0f",
                     format(Sys.Date()), length(unique(cache$pid)),
                     length(cached_pids), added, wc$fmt, wc$kb),
             file.path(OUT, "midweek_status.txt"))
  cat(sprintf("Midweek fetch-only run: +%d points, cache now %d.\n",
              added, length(unique(cache$pid))))
}
cat(sprintf("Done. Spent ~%.0f of %.0f weighted calls this run.\n", spent, wt_cap))

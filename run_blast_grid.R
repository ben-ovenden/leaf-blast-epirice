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
  hw <- get_openmeteo_hourly(lat, lon, start_date, end_date)
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
read_cache <- function() {
  for (f in c(gz_file, csv_file)) {
    if (!file.exists(f)) next
    d <- tryCatch(fread(f, colClasses = list(character = "pid")),
                  error = function(e) {
                    cat(sprintf("Cache read FAILED (%s): %s\n", basename(f), conditionMessage(e))); NULL })
    if (!is.null(d)) {
      cat(sprintf("Cache file: %s (%.0f KB); read %d rows, %d points\n",
                  basename(f), file.info(f)$size / 1024, nrow(d), length(unique(d$pid))))
      return(list(data = d, fmt = if (grepl("\\.gz$", f)) "gz" else "csv"))
    }
  }
  cat("No readable cache found; building fresh.\n")
  list(data = NULL, fmt = NA_character_)
}
cache <- read_cache()$data
if (is.null(cache) || !all(c("infect","semi","wet_hours") %in% names(cache))) {
  if (!is.null(cache)) cat("Cache lacks BLASTAM columns; rebuilding from scratch.\n")
  cache <- data.table(pid = character(), lon = numeric(), lat = numeric(),
                      date = as.Date(character()), TEMP = numeric(), RHUM = numeric(),
                      RAIN = numeric(), wet_hours = numeric(), temp_wet = numeric(),
                      infect = integer(), semi = integer())
}
cache[, date := as.Date(date)]

# ---- Decide what to fetch --------------------------------------------------
cost_new <- (CROP_AGE_DAYS / 7) * 0.75
cached_pids <- unique(cache$pid)
last_by <- if (nrow(cache) > 0) cache[, .(last = max(date)), by = pid] else data.table(pid = character(), last = as.Date(character()))
# Refresh only points stale by >= REFRESH_MIN_STALE_DAYS (most-stale first), capped
# to the per-run budget; skip fresh/current points so their budget adds new points
# instead. Two staggered runs a week then cover the grid ~weekly within the limits.
stale_days <- if (exists("REFRESH_MIN_STALE_DAYS")) REFRESH_MIN_STALE_DAYS else 0L
eligible <- last_by[last < end_date & last <= (end_date - stale_days)][order(last)]
n_refresh <- min(nrow(eligible), TARGET_CALLS_PER_RUN)
maintain <- targets[pid %in% eligible$pid[seq_len(n_refresh)]]
to_add   <- targets[!pid %in% cached_pids]
maintain_cost <- nrow(maintain) * 1
remaining <- max(0, TARGET_CALLS_PER_RUN - maintain_cost)
n_add <- min(nrow(to_add), floor(remaining / cost_new))
add <- if (n_add > 0) to_add[seq_len(n_add)] else to_add[0]
skipped_fresh <- length(cached_pids) - nrow(eligible)
cat(sprintf("Cache: %d points (%d fresh, skipped). Refresh %d, add %d new (~%.0f weighted, budget %d)\n",
            length(cached_pids), skipped_fresh, nrow(maintain), nrow(add),
            maintain_cost + nrow(add) * cost_new, TARGET_CALLS_PER_RUN))

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

do_fetch <- function(tab, start_fun, label, cost) {
  n <- nrow(tab); got <- 0L
  min_chunk_s <- (GRID_CONC * cost) / (GRID_TARGET_PER_MIN / 60)  # hold the rate
  i <- 1L
  while (i <= n) {
    j <- min(i + GRID_CONC - 1L, n)
    t0 <- Sys.time()
    idx <- i:j
    chunk <- tryCatch(
      parallel::mclapply(idx, function(m)
        fetch_point(tab$lon[m], tab$lat[m], start_fun(tab$pid[m])), mc.cores = GRID_CONC),
      error = function(e) lapply(idx, function(m)
        fetch_point(tab$lon[m], tab$lat[m], start_fun(tab$pid[m]))))
    for (r in chunk) if (is.data.frame(r) && nrow(r) > 0) {
      new_rows[[length(new_rows) + 1L]] <<- r; got <- got + 1L
    }
    if (j %% 100 < GRID_CONC || j == n) cat(sprintf("  %s %d/%d (%d ok)\n", label, j, n, got))
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (pace_mult > 0 && el < min_chunk_s) Sys.sleep((min_chunk_s - el) * pace_mult)
    account_and_maybe_pause(length(idx) * cost)
    i <- j + 1L
  }
}
if (nrow(maintain) > 0) {
  last_by <- cache[, .(last = max(date)), by = pid]; setkey(last_by, pid)
  ms <- function(k) { l <- last_by[.(k), last]; if (length(l) == 0 || is.na(l)) emergence else max(emergence, l + 1) }
  do_fetch(maintain, ms, "refresh", 1)
}
if (nrow(add) > 0) do_fetch(add, function(k) emergence, "add", cost_new)

# Serial retry (clean, non-forked connection) for add points that failed under
# concurrency, so as much coverage as possible is filled this run rather than
# waiting for the next weekly run.
fetched_pids <- if (length(new_rows) > 0) unique(rbindlist(new_rows)$pid) else character(0)
failed_add <- add[!pid %in% fetched_pids]
if (nrow(failed_add) > 0) {
  cat(sprintf("  serial retry for %d failed add points\n", nrow(failed_add)))
  for (i in seq_len(nrow(failed_add))) {
    r <- fetch_point(failed_add$lon[i], failed_add$lat[i], emergence)
    if (!is.null(r) && nrow(r) > 0) new_rows[[length(new_rows) + 1L]] <- r
  }
}
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
       events = sum(dt$infect[dt$date > (end_date - win)], na.rm = TRUE))
}
pm <- cache[, { m <- model_pt(.SD); .(lon = lon[1], lat = lat[1],
              intensity = m$intensity, events = m$events) }, by = pid]
pm_epi <- pm[!is.na(intensity)]
map_spacing <- if (nrow(pm) > 0) sqrt(prod(ext[c(2,4)] - ext[c(1,3)]) / nrow(pm)) else NA
cat(sprintf("Modelled %d points; effective spacing ~%.2f deg\n", nrow(pm), map_spacing))

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
    nrow(fread(gz_file, colClasses = list(character = "pid"))) == nrow(dt)
  }, error = function(e) FALSE)
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

# Stats line for the email: points, spacing, last week's points, target spacing,
# cache format and size, so the reader sees the map sharpening and the cache size.
writeLines(sprintf("%d|%.2f|%d|%.2f|%s|%.0f",
                   length(unique(cache$pid)), map_spacing,
                   length(cached_pids), GRID_RES_FINEST, wc$fmt, wc$kb),
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

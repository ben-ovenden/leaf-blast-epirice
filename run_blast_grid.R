#!/usr/bin/env Rscript
################################################################################
# run_blast_grid.R
#
# Continental leaf blast risk HEATMAP with a persistent weather cache and
# DYNAMIC resolution.
#
# Each run:
#   - loads the weather cache (kept in the repo, pruned to the 60-day window);
#   - refreshes the newest days for points it already has (cheap);
#   - spends the remaining API budget ADDING new points (coarse-to-fine), so the
#     map sharpens over successive runs until it reaches GRID_RES_FINEST;
#   - runs EPIRICE per point, interpolates to a surface, and renders the map.
#
# Maps POTENTIAL risk (as if a crop of age CROP_AGE_DAYS grew everywhere).
################################################################################

SCRIPT_DIR <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile), winslash = "/"),
  error = function(e) normalizePath(getwd(), winslash = "/"))

source(file.path(SCRIPT_DIR, "epirice_model.R"))
source(file.path(SCRIPT_DIR, "openmeteo_wth.R"))
source(file.path(SCRIPT_DIR, "blast_config.R"))

suppressPackageStartupMessages({library(data.table); library(terra)})

# Per-point weather fetch hook (start..end for one point). Tests can inject.
if (!exists("POINT_FETCH_FN")) POINT_FETCH_FN <- get_openmeteo_wth

ext <- GRID_EXTENT
key_of <- function(lon, lat) paste0(sprintf("%.4f", lon), "_", sprintf("%.4f", lat))

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
  stop("No Australia polygon (australia_land.geojson missing?). ",
       "Add it or set LAND_ONLY <- FALSE.", call. = FALSE)
if (!is.null(land_poly)) land_poly <- terra::crop(land_poly, terra::ext(ext))

# ---- Target lattice at the finest resolution, ordered coarse-to-fine ------
build_targets <- function() {
  fin <- GRID_RES_FINEST
  g <- as.data.table(expand.grid(
    lon = seq(ext[1], ext[2], by = fin),
    lat = seq(ext[3], ext[4], by = fin)))
  if (!is.null(land_poly)) {
    pts <- terra::vect(as.matrix(g[, .(lon, lat)]), type = "points", crs = "EPSG:4326")
    inside <- !is.na(terra::extract(land_poly, pts)[, 2])
    g <- g[inside]
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

# ---- Weather window -------------------------------------------------------
end_date  <- Sys.Date() - ARCHIVE_LAG_DAYS
emergence <- end_date - CROP_AGE_DAYS
win_dates <- seq(emergence, end_date, by = "day")
cat(sprintf("Window %s to %s (%d days); target lattice %d land points at %.2f deg\n",
            emergence, end_date, length(win_dates), nrow(targets), GRID_RES_FINEST))

# ---- Load cache -----------------------------------------------------------
OUT <- file.path(SCRIPT_DIR, OUTPUT_DIR)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cache_file <- file.path(OUT, WEATHER_CACHE_FILE)
cache <- if (file.exists(cache_file))
  tryCatch(fread(cache_file, colClasses = list(character = "pid")),
           error = function(e) NULL) else NULL
if (is.null(cache) || nrow(cache) == 0)
  cache <- data.table(pid = character(), lon = numeric(), lat = numeric(),
                      date = as.Date(character()), TEMP = numeric(),
                      RHUM = numeric(), RAIN = numeric())
cache[, date := as.Date(date)]

# ---- Decide what to fetch this run ----------------------------------------
cost_new <- (CROP_AGE_DAYS / 7) * 0.75            # weighted calls for full window
cached_keys <- unique(cache$pid)
maintain <- targets[pid %in% cached_keys]         # refresh newest days
to_add   <- targets[!pid %in% cached_keys]        # brand-new points

maintain_cost <- nrow(maintain) * 1               # ~1 call (<=2 weeks) each
remaining <- max(0, TARGET_CALLS_PER_RUN - maintain_cost)
n_add <- min(nrow(to_add), floor(remaining / cost_new))
add <- if (n_add > 0) to_add[seq_len(n_add)] else to_add[0]

cat(sprintf("Cache: %d points. This run: refresh %d, add %d new (est %.0f/%d hourly calls)\n",
            length(cached_keys), nrow(maintain), nrow(add),
            maintain_cost + nrow(add) * cost_new, FREE_HOURLY_CALLS))

# ---- Fetch (per point, only the missing dates) ----------------------------
fetch_point <- function(lon, lat, start) {
  w <- tryCatch(POINT_FETCH_FN(lat = lat, lon = lon,
                               start_date = start, end_date = end_date),
                error = function(e) NULL)
  if (is.null(w) || nrow(w) == 0) return(NULL)
  data.table(pid = key_of(lon, lat), lon = lon, lat = lat,
             date = as.Date(w$YYYYMMDD), TEMP = w$TEMP, RHUM = w$RHUM, RAIN = w$RAIN)
}

new_rows <- list()
# Pace by weighted cost to stay under 600 weighted calls/min (=10/sec); use 8/sec
# for headroom. Refresh calls cost ~1, new points ~6.4. GRID_PACE (tests) scales it.
rate_per_sec <- 6
pace_mult <- if (exists("GRID_PACE")) GRID_PACE else 1
do_fetch <- function(tab, start_fun, label, cost) {
  got <- 0L
  for (i in seq_len(nrow(tab))) {
    st <- start_fun(tab$pid[i])
    r <- fetch_point(tab$lon[i], tab$lat[i], st)
    if (!is.null(r)) { new_rows[[length(new_rows) + 1L]] <<- r; got <- got + 1L }
    if (i %% 100 == 0 || i == nrow(tab))
      cat(sprintf("  %s %d/%d (%d ok)\n", label, i, nrow(tab), got))
    if (pace_mult > 0) Sys.sleep((cost / rate_per_sec) * pace_mult)
  }
}
# maintain: from the day after each point's last cached date
if (nrow(maintain) > 0) {
  last_by_key <- cache[, .(last = max(date)), by = pid]
  setkey(last_by_key, pid)
  maint_start <- function(k) {
    l <- last_by_key[.(k), last]
    if (length(l) == 0 || is.na(l)) emergence else max(emergence, l + 1)
  }
  do_fetch(maintain, maint_start, "refresh", cost = 1)
}
if (nrow(add) > 0) do_fetch(add, function(k) emergence, "add", cost = cost_new)

if (length(new_rows) > 0)
  cache <- rbindlist(list(cache, rbindlist(new_rows)), use.names = TRUE)

# ---- Prune to window, de-duplicate ---------------------------------------
cache <- cache[date >= emergence & date <= end_date]
setorder(cache, pid, date)
cache <- unique(cache, by = c("pid", "date"), fromLast = TRUE)

# ---- Run EPIRICE per point that has enough data ---------------------------
run_tag <- format(Sys.Date(), "%Y-%m-%d")
model_one <- function(dt) {
  if (nrow(dt) < 20) return(NA_real_)
  w <- data.table(YYYYMMDD = dt$date, DOY = as.integer(format(dt$date, "%j")),
                  TEMP = dt$TEMP, RHUM = dt$RHUM, RAIN = dt$RAIN,
                  LAT = dt$lat[1], LON = dt$lon[1])
  setorder(w, YYYYMMDD)
  lb <- tryCatch(predict_leaf_blast(w, emergence = emergence,
                                    duration = as.integer(min(120L, nrow(w)))),
                 error = function(e) NULL)
  if (is.null(lb) || nrow(lb) == 0) NA_real_ else lb$intensity[nrow(lb)]
}
pts_intensity <- cache[, .(lon = lon[1], lat = lat[1], intensity = model_one(.SD)),
                       by = pid]
pts_intensity <- pts_intensity[!is.na(intensity)]
cat(sprintf("Modelled %d points; effective spacing ~%.2f deg\n",
            nrow(pts_intensity),
            if (nrow(pts_intensity) > 0)
              sqrt(prod(ext[c(2, 4)] - ext[c(1, 3)]) / nrow(pts_intensity)) else NA))

# ---- Interpolate to a surface and render ----------------------------------
if (nrow(pts_intensity) < 3) {
  cat("Too few modelled points to render a surface (fetch failed and no cache?).\n",
      "Cache saved for next run; skipping map this run.\n", sep = "")
  fwrite(cache, cache_file)
  quit(save = "no", status = 0)
}
disp_res <- GRID_RES_FINEST / 2
r0 <- terra::rast(xmin = ext[1], xmax = ext[2], ymin = ext[3], ymax = ext[4],
                  resolution = disp_res, crs = "EPSG:4326")
vpts <- terra::vect(as.matrix(pts_intensity[, .(lon, lat)]), type = "points",
                    crs = "EPSG:4326")
vpts$val <- pts_intensity$intensity * 100
r <- tryCatch(
  terra::interpIDW(r0, vpts, field = "val", radius = 6, power = 2, maxPoints = 12),
  error = function(e) terra::rasterize(vpts, r0, field = "val", fun = "mean"))
names(r) <- "leaf_blast_pct"
if (!is.null(land_poly)) r <- terra::mask(r, land_poly)

if (isTRUE(WRITE_GEOTIFF))
  terra::writeRaster(r, file.path(OUT, sprintf("blast_heatmap_%s.tif", run_tag)),
                     overwrite = TRUE)

load_bundled <- function(fname) {
  f <- file.path(SCRIPT_DIR, fname)
  v <- tryCatch(if (file.exists(f)) terra::vect(f) else NULL, error = function(e) NULL)
  if (is.null(v) || nrow(v) == 0) return(NULL)
  tryCatch(terra::crop(v, terra::ext(ext)), error = function(e) v)
}
rivers <- if (isTRUE(SHOW_RIVERS)) load_bundled("australia_rivers.geojson") else NULL
roads  <- if (isTRUE(SHOW_ROADS))  load_bundled("australia_roads.geojson") else NULL

png_file <- file.path(OUT, sprintf("blast_heatmap_%s.png", run_tag))
png(png_file, width = 1000, height = 900, res = 120)
op <- par(mar = c(4, 4, 3, 5))
ramp <- grDevices::colorRampPalette(HEAT_COLOURS)(100)
rmax <- if (!is.null(HEAT_MAX)) HEAT_MAX else {
  m <- terra::global(r, "max", na.rm = TRUE)[1, 1]; if (!is.finite(m) || m <= 0) 1 else m
}
terra::plot(r, col = ramp, range = c(0, rmax), xlab = "Longitude", ylab = "Latitude",
            main = sprintf("Rice leaf blast potential risk (%%)  %s", run_tag),
            plg = list(title = "intensity %"))
if (isTRUE(SHOW_RIVERS) && !is.null(rivers)) try(terra::lines(rivers, col = COL_RIVER, lwd = 0.6), silent = TRUE)
if (isTRUE(SHOW_ROADS)  && !is.null(roads))  try(terra::lines(roads,  col = COL_ROAD,  lwd = 0.5), silent = TRUE)
if (isTRUE(SHOW_COAST)  && !is.null(land_poly)) try(terra::lines(land_poly, col = COL_COAST, lwd = 1), silent = TRUE)
if (isTRUE(SHOW_TOWNS) && exists("MONITOR_TOWNS") && nrow(MONITOR_TOWNS) > 0) {
  points(MONITOR_TOWNS$lon, MONITOR_TOWNS$lat, pch = 21, bg = "white", col = COL_TOWN, cex = 1.1, lwd = 1.4)
  text(MONITOR_TOWNS$lon, MONITOR_TOWNS$lat, MONITOR_TOWNS$name, pos = 4, offset = 0.3, cex = 0.5, col = COL_TOWN)
}
par(op); dev.off()
file.copy(png_file, file.path(OUT, "blast_heatmap_latest.png"), overwrite = TRUE)

# ---- Save cache -----------------------------------------------------------
fwrite(cache, cache_file)   # .csv.gz -> fwrite compresses by extension
cat(sprintf("Cache saved: %d points x up to %d days -> %s\n",
            length(unique(cache$pid)), length(win_dates), cache_file))
cat("Heatmap:", png_file, "\nDone.\n")

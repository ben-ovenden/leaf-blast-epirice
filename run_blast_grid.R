#!/usr/bin/env Rscript
################################################################################
# run_blast_grid.R
#
# Continental leaf blast risk HEATMAP. Runs the EPIRICE model (Savary et al.
# 2012) on a grid of points from Open-Meteo weather, assembles a raster surface,
# masks it to the coastline so all land is covered, and renders a heatmap with
# coastline, rivers, roads and towns.
#
#   Rscript run_blast_grid.R
#
# This maps POTENTIAL risk: it colours the land as if a crop of age
# CROP_AGE_DAYS were growing everywhere, given the recent weather. It is a
# regional potential signal, not a map of actual crop or measured disease.
################################################################################

SCRIPT_DIR <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile), winslash = "/"),
  error = function(e) normalizePath(getwd(), winslash = "/"))

source(file.path(SCRIPT_DIR, "epirice_model.R"))
source(file.path(SCRIPT_DIR, "openmeteo_wth.R"))
source(file.path(SCRIPT_DIR, "blast_config.R"))

suppressPackageStartupMessages({library(data.table); library(terra)})

if (!exists("GRID_FETCH_FN")) GRID_FETCH_FN <- get_openmeteo_grid

ext <- GRID_EXTENT

# Land polygon. Try the bundled GeoJSON first, read with terra alone (no sf,
# which can be broken on some runners). Fall back to package sources only if the
# bundled file is missing. Returns a SpatVector or NULL.
get_land <- function() {
  # 0) bundled file, terra-only (most reliable on CI)
  f <- file.path(SCRIPT_DIR, "australia_land.geojson")
  v <- tryCatch(if (file.exists(f)) terra::vect(f) else NULL,
                error = function(e) NULL)
  if (!is.null(v) && nrow(v) > 0) return(v)

  if (requireNamespace("sf", quietly = TRUE)) try(sf::sf_use_s2(FALSE), silent = TRUE)

  # 1) ozmaps
  v <- tryCatch({
    if (requireNamespace("ozmaps", quietly = TRUE) &&
        requireNamespace("sf", quietly = TRUE))
      terra::vect(sf::st_union(sf::st_make_valid(ozmaps::ozmap_country)))
    else NULL
  }, error = function(e) NULL)
  if (!is.null(v)) return(v)

  # 2) rnaturalearth
  v <- tryCatch({
    if (requireNamespace("rnaturalearth", quietly = TRUE))
      terra::vect(rnaturalearth::ne_countries(country = "Australia",
                                              returnclass = "sf"))
    else NULL
  }, error = function(e) NULL)
  if (!is.null(v)) return(v)

  # 3) maps
  v <- tryCatch({
    if (requireNamespace("maps", quietly = TRUE) &&
        requireNamespace("sf", quietly = TRUE)) {
      m <- maps::map("world", "Australia", plot = FALSE, fill = TRUE)
      terra::vect(sf::st_make_valid(sf::st_as_sf(m)))
    } else NULL
  }, error = function(e) NULL)
  v
}
land_poly <- if (isTRUE(LAND_ONLY)) get_land() else NULL

if (isTRUE(LAND_ONLY) && is.null(land_poly)) {
  stop("LAND_ONLY is TRUE but no Australia polygon could be built (ozmaps, ",
       "rnaturalearth and maps all unavailable). Without the mask the run would ",
       "fetch the whole ocean box and exceed the free API budget. Install one of ",
       "those packages, or set LAND_ONLY <- FALSE deliberately.", call. = FALSE)
}

# ---- Build grid: keep every cell that TOUCHES land (no coastal gaps) -------
r0 <- terra::rast(xmin = ext[1], xmax = ext[2], ymin = ext[3], ymax = ext[4],
                  resolution = GRID_RES, crs = "EPSG:4326")

if (!is.null(land_poly)) {
  land_poly <- terra::crop(land_poly, terra::ext(r0))
  cover <- terra::rasterize(land_poly, r0, cover = TRUE)
  cv <- terra::values(cover)
  keep_cells <- which(!is.na(cv) & cv > 0)
} else {
  keep_cells <- seq_len(terra::ncell(r0))
}
xy <- terra::xyFromCell(r0, keep_cells)
grid <- data.table(cell = keep_cells, lon = xy[, 1], lat = xy[, 2])

cat(sprintf("Grid: %d cells at %.2f deg over [%g,%g]x[%g,%g]%s\n",
            nrow(grid), GRID_RES, ext[1], ext[2], ext[3], ext[4],
            if (!is.null(land_poly)) " (land, incl. coast)" else ""))

# ---- Budget self-check -----------------------------------------------------
.free_calls <- if (exists("FREE_DAILY_CALLS")) FREE_DAILY_CALLS else 10000
calls_per_cell <- (CROP_AGE_DAYS / 7) * 0.75
est_calls <- nrow(grid) * calls_per_cell
cat(sprintf("Estimated API usage: %.0f call units (%.0f%% of %d/day free)\n",
            est_calls, 100 * est_calls / .free_calls, .free_calls))
if (is.null(land_poly) && isTRUE(LAND_ONLY))
  cat("  note: ocean not masked (install ozmaps for land-only, which cuts this).\n")
if (est_calls > .free_calls)
  warning(sprintf("Estimated usage %.0f exceeds free daily allowance (%d).",
                  est_calls, .free_calls), call. = FALSE)

# ---- Weather window --------------------------------------------------------
end_date  <- Sys.Date() - ARCHIVE_LAG_DAYS
emergence <- end_date - CROP_AGE_DAYS
cat(sprintf("Weather %s to %s  (crop age %d days)\n",
            emergence, end_date, CROP_AGE_DAYS))

# ---- Fetch and model -------------------------------------------------------
cat("Fetching gridded weather from Open-Meteo ...\n")
wth_list <- GRID_FETCH_FN(lats = grid$lat, lons = grid$lon,
                          start_date = emergence, end_date = end_date)

cat("Running EPIRICE leaf blast per cell ...\n")
grid[, intensity := NA_real_]
ok <- 0L
for (k in seq_len(nrow(grid))) {
  w <- wth_list[[k]]
  if (is.null(w) || nrow(w) < 20) next
  dur <- as.integer(min(120L, nrow(w)))
  lb <- tryCatch(predict_leaf_blast(w, emergence = emergence, duration = dur),
                 error = function(e) NULL)
  if (!is.null(lb) && nrow(lb) > 0) {
    grid$intensity[k] <- lb$intensity[nrow(lb)]
    ok <- ok + 1L
  }
}
cat(sprintf("Modelled %d of %d cells\n", ok, nrow(grid)))

# ---- Assemble raster, clip to coast ---------------------------------------
OUT <- file.path(SCRIPT_DIR, OUTPUT_DIR)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
run_tag <- format(Sys.Date(), "%Y-%m-%d")

r <- r0
terra::values(r) <- NA_real_
r[grid$cell] <- grid$intensity * 100
names(r) <- "leaf_blast_pct"
if (!is.null(land_poly)) r <- terra::mask(r, land_poly)

if (isTRUE(WRITE_GEOTIFF))
  terra::writeRaster(r, file.path(OUT, sprintf("blast_heatmap_%s.tif", run_tag)),
                     overwrite = TRUE)

# ---- Monitoring-town values and rolling trends CSV ------------------------
# Sample the risk surface at each monitoring town, then keep a wide table of the
# last HISTORY_RUNS runs (one column per run date) so trends are visible per row.
if (exists("MONITOR_TOWNS") && nrow(MONITOR_TOWNS) > 0) {
  mt <- as.data.table(MONITOR_TOWNS)
  tpts <- terra::vect(as.matrix(mt[, .(lon, lat)]), type = "points",
                      crs = "EPSG:4326")
  tvals <- terra::extract(r, tpts)[, 2]
  today <- data.table(town = mt$name)
  today[[run_tag]] <- round(tvals, 3)

  hist_file <- file.path(OUT, "town_trends.csv")
  hist <- if (file.exists(hist_file))
    fread(hist_file, header = TRUE, colClasses = list(character = "town")) else
    data.table(town = character())
  if (run_tag %in% names(hist)) hist[, (run_tag) := NULL]   # re-run same day

  hist <- merge(hist, today, by = "town", all = TRUE)

  # rows in MONITOR_TOWNS order (any retired towns kept at the end)
  ord <- c(mt$name, setdiff(hist$town, mt$name))
  hist <- hist[match(ord, town)]

  # keep only the most recent HISTORY_RUNS date columns, in date order
  keep_n <- if (exists("HISTORY_RUNS")) HISTORY_RUNS else 10L
  date_cols <- setdiff(names(hist), "town")
  date_cols <- date_cols[order(as.Date(date_cols))]
  if (length(date_cols) > keep_n) date_cols <- tail(date_cols, keep_n)
  hist <- hist[, c("town", date_cols), with = FALSE]

  fwrite(hist, hist_file)
  cat(sprintf("Town trends: %d towns x %d runs -> %s\n",
              nrow(hist), length(date_cols), hist_file))
}

r_disp <- r
if (SMOOTH_FACTOR > 1L) {
  r_disp <- terra::disagg(r, fact = SMOOTH_FACTOR, method = "bilinear")
  if (!is.null(land_poly)) r_disp <- terra::mask(r_disp, land_poly)
}

# ---- Bundled Natural Earth overlays, read with terra (no sf, no download) ---
load_bundled <- function(fname) {
  f <- file.path(SCRIPT_DIR, fname)
  v <- tryCatch(if (file.exists(f)) terra::vect(f) else NULL,
                error = function(e) NULL)
  if (is.null(v) || nrow(v) == 0) return(NULL)
  tryCatch(terra::crop(v, terra::ext(r0)), error = function(e) v)
}
rivers <- if (isTRUE(SHOW_RIVERS)) load_bundled("australia_rivers.geojson") else NULL
roads  <- if (isTRUE(SHOW_ROADS))  load_bundled("australia_roads.geojson") else NULL

# ---- Render ----------------------------------------------------------------
png_file <- file.path(OUT, sprintf("blast_heatmap_%s.png", run_tag))
png(png_file, width = 1000, height = 900, res = 120)
op <- par(mar = c(4, 4, 3, 5))
ramp <- grDevices::colorRampPalette(HEAT_COLOURS)(100)
rmax <- if (!is.null(HEAT_MAX)) HEAT_MAX else {
  m <- terra::global(r, "max", na.rm = TRUE)[1, 1]
  if (!is.finite(m) || m <= 0) 1 else m
}
terra::plot(r_disp, col = ramp, range = c(0, rmax),
            xlab = "Longitude", ylab = "Latitude",
            main = sprintf("Rice leaf blast potential risk (%%)  %s", run_tag),
            plg = list(title = "intensity %"))

if (isTRUE(SHOW_RIVERS) && !is.null(rivers))
  try(terra::lines(rivers, col = COL_RIVER, lwd = 0.6), silent = TRUE)
if (isTRUE(SHOW_ROADS) && !is.null(roads))
  try(terra::lines(roads, col = COL_ROAD, lwd = 0.5), silent = TRUE)
if (isTRUE(SHOW_COAST) && !is.null(land_poly))
  try(terra::lines(land_poly, col = COL_COAST, lwd = 1), silent = TRUE)

if (isTRUE(SHOW_TOWNS) && exists("MONITOR_TOWNS") && nrow(MONITOR_TOWNS) > 0) {
  points(MONITOR_TOWNS$lon, MONITOR_TOWNS$lat, pch = 21, bg = "white",
         col = COL_TOWN, cex = 1.1, lwd = 1.4)
  text(MONITOR_TOWNS$lon, MONITOR_TOWNS$lat, MONITOR_TOWNS$name, pos = 4,
       offset = 0.3, cex = 0.5, col = COL_TOWN)
}

par(op); dev.off()
file.copy(png_file, file.path(OUT, "blast_heatmap_latest.png"), overwrite = TRUE)

cat("Heatmap: ", png_file, "\n")
if (isTRUE(WRITE_GEOTIFF))
  cat("GeoTIFF: ", file.path(OUT, sprintf("blast_heatmap_%s.tif", run_tag)), "\n")
cat("\nDone.\n")

#!/usr/bin/env Rscript
################################################################################
# run_blast_grid.R
#
# Continental leaf blast risk HEATMAP. Runs the EPIRICE model (Savary et al.
# 2012) on a grid of points from Open-Meteo weather, assembles a raster surface,
# and renders a heatmap with an Australian coastline and your monitoring sites.
#
#   Rscript run_blast_grid.R
#
# This maps POTENTIAL risk: it colours the extent as if a crop of age
# CROP_AGE_DAYS were growing at every cell, given the recent weather. It is a
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

# ---- Build grid -----------------------------------------------------------
ext <- GRID_EXTENT
lon_seq <- seq(ext[1], ext[2], by = GRID_RES)
lat_seq <- seq(ext[3], ext[4], by = GRID_RES)
grid <- as.data.table(expand.grid(lon = lon_seq, lat = lat_seq))

# Optional land mask: keep cells inside an Australia polygon
land_poly <- NULL
if (isTRUE(LAND_ONLY)) {
  land_poly <- tryCatch({
    if (requireNamespace("ozmaps", quietly = TRUE)) {
      terra::vect(sf::st_union(ozmaps::ozmap_country))
    } else if (requireNamespace("rnaturalearth", quietly = TRUE)) {
      terra::vect(rnaturalearth::ne_countries(country = "Australia",
                                              returnclass = "sf"))
    } else NULL
  }, error = function(e) NULL)

  if (!is.null(land_poly)) {
    pts <- terra::vect(as.matrix(grid[, .(lon, lat)]), type = "points",
                       crs = "EPSG:4326")
    inside <- !is.na(terra::extract(land_poly, pts)[, 2])
    grid <- grid[inside]
  }
}
cat(sprintf("Grid: %d cells at %.2f deg over [%g,%g]x[%g,%g]%s\n",
            nrow(grid), GRID_RES, ext[1], ext[2], ext[3], ext[4],
            if (!is.null(land_poly)) " (land only)" else ""))

# ---- Budget self-check ----------------------------------------------------
# Open-Meteo counts long windows as fractional extra calls: about 0.75 of a call
# per week of data per cell. Warn before spending if the run exceeds the free
# daily allowance.
.free_calls <- if (exists("FREE_DAILY_CALLS")) FREE_DAILY_CALLS else 10000
calls_per_cell <- (CROP_AGE_DAYS / 7) * 0.75
est_calls <- nrow(grid) * calls_per_cell
cat(sprintf("Estimated API usage: %.0f call units (%.0f%% of %d/day free)\n",
            est_calls, 100 * est_calls / .free_calls, .free_calls))
if (is.null(land_poly) && isTRUE(LAND_ONLY)) {
  cat("  note: ocean not masked (install ozmaps for land-only, which cuts this).\n")
}
if (est_calls > .free_calls) {
  warning(sprintf(paste0("Estimated usage %.0f exceeds the free daily allowance ",
                         "(%d). Coarsen GRID_RES, shrink GRID_EXTENT, or lower ",
                         "CROP_AGE_DAYS."), est_calls, .free_calls), call. = FALSE)
}

# ---- Weather window -------------------------------------------------------
end_date   <- Sys.Date() - ARCHIVE_LAG_DAYS
emergence  <- end_date - CROP_AGE_DAYS
cat(sprintf("Weather %s to %s  (crop age %d days)\n",
            emergence, end_date, CROP_AGE_DAYS))

# ---- Fetch weather for all cells -----------------------------------------
cat("Fetching gridded weather from Open-Meteo ...\n")
wth_list <- GRID_FETCH_FN(lats = grid$lat, lons = grid$lon,
                          start_date = emergence, end_date = end_date)

# ---- Run EPIRICE per cell -------------------------------------------------
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

# ---- Assemble raster ------------------------------------------------------
OUT <- file.path(SCRIPT_DIR, OUTPUT_DIR)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
run_tag <- format(Sys.Date(), "%Y-%m-%d")

r <- terra::rast(xmin = ext[1] - GRID_RES/2, xmax = ext[2] + GRID_RES/2,
                 ymin = ext[3] - GRID_RES/2, ymax = ext[4] + GRID_RES/2,
                 resolution = GRID_RES, crs = "EPSG:4326")
r <- terra::rasterize(as.matrix(grid[, .(lon, lat)]), r,
                      values = grid$intensity * 100, fun = "mean")
names(r) <- "leaf_blast_pct"

if (isTRUE(WRITE_GEOTIFF)) {
  terra::writeRaster(r, file.path(OUT, sprintf("blast_heatmap_%s.tif", run_tag)),
                     overwrite = TRUE)
}

# Smooth only for display
r_disp <- r
if (SMOOTH_FACTOR > 1L) {
  r_disp <- terra::disagg(r, fact = SMOOTH_FACTOR, method = "bilinear")
}

# ---- Render heatmap PNG ---------------------------------------------------
sites <- as.data.table(SITES)
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
# Coastline overlay
if (!is.null(land_poly)) {
  terra::lines(land_poly, col = "#22272B", lwd = 1)
} else if (requireNamespace("maps", quietly = TRUE)) {
  maps::map("world", "Australia", add = TRUE, col = "#22272B", lwd = 1)
}
# Monitoring sites
points(sites$lon, sites$lat, pch = 21, bg = "white", col = "#22272B",
       cex = 1.2, lwd = 1.5)
text(sites$lon, sites$lat, sites$name, pos = 3, cex = 0.7, col = "#22272B")
par(op); dev.off()
file.copy(png_file, file.path(OUT, "blast_heatmap_latest.png"), overwrite = TRUE)

cat("Heatmap: ", png_file, "\n")
if (isTRUE(WRITE_GEOTIFF))
  cat("GeoTIFF: ", file.path(OUT, sprintf("blast_heatmap_%s.tif", run_tag)), "\n")
cat("\nDone.\n")

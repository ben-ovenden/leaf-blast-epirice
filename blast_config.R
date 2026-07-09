################################################################################
# blast_config.R
#
# Edit this file to configure the weekly leaf blast run, then commit to git.
################################################################################

# ---- 1. Monitoring sites ---------------------------------------------------
# name, latitude, longitude. Optionally add an `emergence` column (YYYY-MM-DD)
# to set a different transplanting/emergence date per site; otherwise the single
# EMERGENCE_DATE below is used for all sites.
SITES <- data.frame(
  name = c("Cairns", "Townsville", "Mackay", "Bundaberg"),
  lat  = c(-16.88, -19.25, -21.14, -24.87),
  lon  = c(145.75, 146.82, 149.19, 152.35),
  stringsAsFactors = FALSE
)

# ---- 2. Season -------------------------------------------------------------
# Emergence / transplanting date. The model simulates from here up to the most
# recent weather available, capped at the 120 day season length.
EMERGENCE_DATE <- "2025-12-01"

# ---- 3. Risk bands ---------------------------------------------------------
# EPIRICE reports leaf blast intensity as a proportion of diseased sites (0-1).
# These cut points are PROVISIONAL and should be calibrated against your own
# field observations before being relied on. They only affect labelling.
INTENSITY_LOW_MAX      <- 0.05   # below this = low
INTENSITY_MODERATE_MAX <- 0.15   # below this = moderate, at or above = high

# ---- 4. Data window --------------------------------------------------------
# The Open-Meteo archive (ERA5) lags real time by about 5 days. End the fetch a
# few days back so the most recent days are actually populated.
ARCHIVE_LAG_DAYS <- 5

# Minimum days of weather after emergence before a run is meaningful. Disease
# onset in the leaf blast model is day 15.
MIN_DAYS <- 16

# ---- 5. Output -------------------------------------------------------------
OUTPUT_DIR <- "blast_outputs"
MAKE_MAP   <- TRUE          # write a simple PNG risk map
MAP_WIDTH  <- 900
MAP_HEIGHT <- 800

# NSW Government brand palette (https://designsystem.nsw.gov.au)
COL_LOW      <- "#00AA45"   # NSW green (Success)
COL_MODERATE <- "#DC5800"   # NSW orange (Warning)
COL_HIGH     <- "#B81237"   # NSW red (Error)
COL_NODATA   <- "#C8CDD0"   # NSW grey

# ---- 6. Email (optional) ---------------------------------------------------
# Store credentials as environment variables / GitHub Secrets, never in git.
EMAIL_FROM     <- Sys.getenv("BLAST_EMAIL_USER", "")
EMAIL_PASSWORD <- Sys.getenv("BLAST_EMAIL_PASS", "")
EMAIL_TO       <- Sys.getenv("BLAST_EMAIL_TO", "")
EMAIL_SMTP_HOST <- "smtp.gmail.com"
EMAIL_SMTP_PORT <- 587
SEND_EMAIL <- FALSE          # set TRUE once secrets are configured

################################################################################
# 7. Heatmap grid (for run_blast_grid.R)
################################################################################
# The heatmap runs the EPIRICE model on a grid of points and renders the result
# as a continuous risk surface, in the spirit of Savary et al. (2012), who
# mapped POTENTIAL epidemics. It colours the whole extent as if rice were grown
# everywhere, so most of the continent is a potential signal, not actual crop.

# Extent: c(lon_min, lon_max, lat_min, lat_max). Default = Australian continent.
GRID_EXTENT <- c(112, 154, -44, -10)

# Resolution in degrees. 0.75 is tuned to keep a land-only continental run (about
# 1,240 cells over a 60 day window) at roughly 80% of the free Open-Meteo budget
# (10,000 calls/day), with headroom. Finer than ~0.70 over the whole continent
# exceeds the free tier; go finer only over a smaller extent.
GRID_RES <- 0.75

# Free Open-Meteo daily call allowance, used by the budget self-check.
FREE_DAILY_CALLS <- 10000

# Crop age framing: each week, assume a crop this many days old everywhere, and
# run the model over the trailing weather window. This makes the surface
# spatially comparable and updates with recent weather. Also bounds API cost.
CROP_AGE_DAYS <- 60L

# Drop ocean cells using an Australia land polygon (needs ozmaps or rnaturalearth
# at run time). If unavailable, all cells are kept.
LAND_ONLY <- TRUE

# Display smoothing factor for the rendered PNG (model still runs at GRID_RES;
# this only resamples for a smoother-looking image). 1 = no smoothing.
SMOOTH_FACTOR <- 4L

# Also write a GeoTIFF of the risk surface alongside the PNG.
WRITE_GEOTIFF <- TRUE

# Heatmap colour ramp (low -> high), NSW brand-aligned.
HEAT_COLOURS <- c("#00AA45", "#FFB223", "#DC5800", "#B81237")

# Heatmap colour scale maximum (%). NULL = auto-scale to the current week's data
# (always readable, but colours are not comparable between weeks). Set a fixed
# number once you know your typical range to make weeks comparable.
HEAT_MAX <- NULL

################################################################################
# 8. Citation (added to the summary and to the heatmap footer)
################################################################################
# The implemented disease model is EPIRICE (leaf blast), via Adam Sparks' epicrop
# package. BLASTL is cited as the related Japanese epidemic-progression model
# (conceptual lineage), not as the implemented code.
CITATION <- paste(
  "Model and data citation",
  "------------------------------------------------------------",
  "Model: EPIRICE. Savary, S., Nelson, A., Willocquet, L., Pangga, I.,",
  "  and Aunario, J. (2012). Modeling and mapping potential epidemics of",
  "  rice diseases globally. Crop Protection 34: 6-17.",
  "  doi:10.1016/j.cropro.2011.11.009",
  "Software: epicrop R package (leaf blast parameters and SEIR engine),",
  "  Adam H. Sparks and colleagues; adapted from cropsim (Hijmans et al. 2009).",
  "SEIR framework: Zadoks, J.C. (1971). Phytopathology 61: 600-610.",
  "Weather: Open-Meteo ERA5 archive (data CC BY 4.0).",
  "Related model (BLASTL, conceptual lineage): Hashimoto, A., Hirano, K.,",
  "  and Matsumoto, K. (1984). Studies on the forecasting of rice leaf blast",
  "  development by application of the computer simulation. Special Bulletin",
  "  of Fukushima Prefecture Agricultural Experiment Station 2: 1-104.",
  sep = "\n"
)

# Compact one/two-line credit for the map footer
CITATION_MAP <- paste(
  "Model: EPIRICE (Savary et al. 2012, Crop Prot. 34:6-17; doi:10.1016/j.cropro.2011.11.009).",
  "Code: epicrop (A.H. Sparks). Weather: Open-Meteo ERA5 (CC BY 4.0). Lineage: BLASTL (Hashimoto et al. 1984).",
  sep = "\n"
)

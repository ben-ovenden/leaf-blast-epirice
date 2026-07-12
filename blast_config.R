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
# Calibrated so that 1% (0.01) intensity, taken as a significant epidemic, is the
# top of the "high" band. Adjust these against your own field observations.
INTENSITY_LOW_MAX      <- 0.002  # below 0.2% = low
INTENSITY_MODERATE_MAX <- 0.01   # 0.2% to <1% = moderate; 1% and above = high

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

# Resolution. The map now refines over consecutive runs: a weather cache
# (WEATHER_CACHE_FILE, kept in the repo) means each run only fetches the newest
# days for points it already has, which frees budget to ADD new points. Coverage
# fills coarse-to-fine toward GRID_RES_FINEST, then holds there once maintaining
# the points costs about the per-run budget.
GRID_RES_FINEST <- 0.5              # finest resolution the map refines toward
GRID_RES_LEVELS <- c(2.0, 1.0, 0.5) # coarse-to-fine fill order; each must be a
                                    # whole-number multiple of GRID_RES_FINEST
TARGET_CALLS_PER_RUN <- 4000        # weighted-call budget per run (< 5000/hour)
GRID_CONC           <- 4            # grid points fetched concurrently (Linux runner)
GRID_TARGET_PER_MIN <- 400          # weighted-call rate cap while fetching (< 600/min)
WEATHER_CACHE_FILE   <- "weather_cache.csv"     # plain CSV (git-safe), committed

# Legacy fixed resolution (no longer used by the dynamic runner; kept for
# reference and any single-shot use).
GRID_RES <- 1.25

# Free Open-Meteo limits used by the budget self-check. The hourly cap is the
# binding one for a single weekly run.
FREE_HOURLY_CALLS <- 5000
FREE_DAILY_CALLS  <- 10000

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

# Heatmap colour ramp, low -> high: light blue through to dark red.
HEAT_COLOURS <- c("#BFE0F5", "#FFF6B0", "#FDB147", "#E8492B", "#8B0000")

# BLASTAM infection-days heatmap ramp, low -> high: light blue through to dark
# red (same scheme as EPIRICE; the title and legend distinguish the two maps).
BLASTAM_HEAT_COLOURS <- c("#BFE0F5", "#FFF6B0", "#FDB147", "#E8492B", "#8B0000")

# BLASTAM reporting window: the map, town table and trends count FAVOURABLE
# infection days over the most recent BLASTAM_WINDOW_DAYS, to show where infection
# pressure is building now (not a whole-season total). EPIRICE keeps its own
# CROP_AGE_DAYS window; only the cached daily flags are shared.
BLASTAM_WINDOW_DAYS <- 21

# Deepest colour at this many favourable days within the 21-day window (the max
# possible is the window length). PROVISIONAL: lower it (e.g. 14) to make
# building pressure show up more strongly. NULL auto-scales each week.
BLASTAM_HEAT_MAX <- 21

# Overlay layers on the heatmap
SHOW_COAST  <- TRUE
SHOW_TOWNS  <- TRUE
SHOW_ROADS  <- TRUE    # major roads from Natural Earth (downloaded at run time)
SHOW_RIVERS <- TRUE    # rivers from Natural Earth (downloaded at run time)

COL_COAST <- "#22272B"
COL_ROAD  <- "#8A6D3B"
COL_RIVER <- "#2E75B6"
COL_TOWN  <- "#111111"

# Monitoring towns: modelled each run, tracked in the trends CSV, and highlighted
# on the map. name, lon, lat. Coordinates for the remote roadhouses (Archer
# River, Lakeland, Timber Creek) are approximate; edit any as needed.
MONITOR_TOWNS <- data.frame(
  name = c("Malanda", "Lismore", "Gympie", "Dalby", "Warwick", "Griffith",
           "Deniliquin", "Moree", "Emerald", "Clermont", "Biloela", "Marian",
           "Proserpine", "Home Hill", "Lakeland", "Archer River", "Bamaga",
           "Croydon", "Burketown", "Borroloola", "Katherine", "Humpty Doo",
           "Jabiru", "Timber Creek", "Kununurra",
           "Roma", "Rolleston", "Goondiwindi", "Theodore",
           "Fitzroy Crossing", "Dubbo"),
  state = c("QLD", "NSW", "QLD", "QLD", "QLD", "NSW",
            "NSW", "NSW", "QLD", "QLD", "QLD", "QLD",
            "QLD", "QLD", "QLD", "QLD", "QLD",
            "QLD", "QLD", "NT", "NT", "NT",
            "NT", "NT", "WA",
            "QLD", "QLD", "QLD", "QLD",
            "WA", "NSW"),
  lon  = c(145.596, 153.277, 152.665, 151.262, 152.034, 146.040,
           144.958, 149.841, 148.159, 147.639, 150.503, 148.947,
           148.581, 147.415, 144.851, 142.940, 142.386,
           142.240, 139.546, 136.307, 132.264, 131.281,
           132.836, 130.480, 128.741,
           148.787, 148.630, 150.310, 150.076,
           125.565, 148.601),
  lat  = c(-17.354, -28.814, -26.190, -27.183, -28.219, -34.288,
           -35.532, -29.462, -23.527, -22.826, -24.403, -21.150,
           -20.401, -19.663, -15.855, -13.435, -10.892,
           -18.204, -17.744, -16.070, -14.465, -12.584,
           -12.671, -15.660, -15.772,
           -26.570, -24.462, -28.549, -24.948,
           -18.197, -32.257),
  stringsAsFactors = FALSE
)

# Number of recent runs to keep in the trends CSV (one column per run).
HISTORY_RUNS <- 10

# Town weather fetches run concurrently (they are slow, cold, independent hourly
# requests). This caps the concurrency. Linux only (the GitHub runner); Windows
# runs serial automatically.
TOWN_FETCH_CORES <- 4L

# Heatmap colour scale maximum (%). A FIXED number makes every week's colours
# directly comparable (a given colour = the same intensity each week), which is
# usually what you want for a weekly product. NULL would auto-scale each week.
#
# Set to 2 so 2% intensity and above is the deepest red (a severe epidemic);
# ~1% reads as strong orange-red. Values grade from light blue up.
HEAT_MAX <- 2

################################################################################
# 8. Citation (added to the summary and to the heatmap footer)
################################################################################
# The implemented disease model is EPIRICE (leaf blast), via Adam Sparks' epicrop
# package. BLASTL is cited as the related Japanese epidemic-progression model
# (conceptual lineage), not as the implemented code.
CITATION <- paste(
  "Models and data citation",
  "------------------------------------------------------------",
  "EPIRICE (intensity): Savary, S., Nelson, A., Willocquet, L., Pangga, I.,",
  "  and Aunario, J. (2012). Modeling and mapping potential epidemics of",
  "  rice diseases globally. Crop Protection 34: 6-17.",
  "  doi:10.1016/j.cropro.2011.11.009",
  "  Software: epicrop R package (SEIR engine), A.H. Sparks and colleagues,",
  "  adapted from cropsim (Hijmans et al. 2009). Framework: Zadoks (1971).",
  "BLASTAM (infection days): infection-warning model of Koshimizu, Y. (1988),",
  "  A forecasting method for occurrence of rice leaf blast with AMeDAS data,",
  "  Bull. Tohoku Natl. Agric. Exp. Stn. 78: 67-121; program in Hayashi, T. and",
  "  Koshimizu, Y. (1988), ibid. 78: 123-138 [in Japanese]. A day is favourable",
  "  when leaf wetness >=10 h, the mean temperature during wetness is within",
  "  bounds, and the preceding 5-day mean is within bounds. Original Japanese",
  "  bounds 15-25 C and 20-25 C; the upper bounds are RAISED here to 15-32 C and",
  "  20-30 C for warmer Australian conditions (a deliberate deviation, since",
  "  blast infects up to ~32 C given leaf wetness). Leaf wetness estimated from",
  "  hourly humidity and rainfall.",
  "Weather: Open-Meteo ERA5 archive (data CC BY 4.0).",
  sep = "\n"
)

# Compact one/two-line credit for the map footer
CITATION_MAP <- paste(
  "Model: EPIRICE (Savary et al. 2012, Crop Prot. 34:6-17; doi:10.1016/j.cropro.2011.11.009).",
  "Code: epicrop (A.H. Sparks). Weather: Open-Meteo ERA5 (CC BY 4.0). Lineage: BLASTL (Hashimoto et al. 1984).",
  sep = "\n"
)

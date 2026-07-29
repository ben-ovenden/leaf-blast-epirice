################################################################################
# blast_config.R
#
# Edit this file to configure the weekly leaf blast run, then commit to git.
#
# UNITS WARNING. Two different quantities are easy to confuse and confusing them
# was the cause of several past bugs:
#   * a FETCH is one grid point brought up to date;
#   * a WEIGHTED CALL is what Open-Meteo actually charges.
# They are not the same. One 61 day fetch of three hourly variables costs about
# 4.4 weighted calls, not 1. Every budget below is labelled with its unit.
################################################################################

# ---- 1. Season -------------------------------------------------------------
# Risk bands. EPIRICE reports leaf blast intensity as a proportion of diseased
# sites (0-1). Calibrated so 1% (0.01) intensity, taken as a significant
# epidemic, is the top of the "high" band. Adjust against field observations.
INTENSITY_LOW_MAX      <- 0.002  # below 0.2% = low
INTENSITY_MODERATE_MAX <- 0.01   # 0.2% to <1% = moderate; 1% and above = high

# Crop age framing: each run, assume a crop this many days old everywhere and
# run the model over the trailing weather window. Makes the surface spatially
# comparable and bounds the API cost.
#
# NOTE for the trends CSVs: because emergence is (end_date - CROP_AGE_DAYS), it
# moves with every run. The trends columns are therefore a ROLLING 60 day window,
# not a cumulative season total, and must not be read as season progress.
CROP_AGE_DAYS <- 60L

# Minimum days of weather after emergence before a run is meaningful. Disease
# onset in the leaf blast model is day 15. With the rolling emergence above this
# never binds, but it still guards a fixed per site emergence override.
MIN_DAYS <- 16

# ---- 2. Data window --------------------------------------------------------
# The Open-Meteo archive (ERA5) lags real time by about 5 days. End the fetch a
# few days back so the most recent days are actually populated.
ARCHIVE_LAG_DAYS <- 6

# ---- 3. Output -------------------------------------------------------------
OUTPUT_DIR <- "blast_outputs"
MAKE_MAP   <- FALSE         # off: the workflow emails the two grid heatmaps;
                            # this simple town point map is not delivered
MAP_WIDTH  <- 900
MAP_HEIGHT <- 800

# NSW Government brand palette (https://designsystem.nsw.gov.au)
COL_LOW      <- "#00AA45"   # NSW green (Success)
COL_MODERATE <- "#DC5800"   # NSW orange (Warning)
COL_HIGH     <- "#B81237"   # NSW red (Error)
COL_NODATA   <- "#C8CDD0"   # NSW grey

################################################################################
# 4. Weather source
################################################################################
# Hourly variables fed to blastam_daily_from_hourly(). KEEP THIS AT TEN OR FEWER:
# at eleven the weighted cost per location gains a variable multiplier.
OM_HOURLY_VARS <- c("temperature_2m", "relative_humidity_2m", "precipitation")

# Pin the reanalysis explicitly. Leaving models= unset lets the API choose, and
# it can serve different variables from different products, so temperature might
# come from ERA5-Land (0.1 deg) while humidity comes from ERA5 (0.25 deg). That
# makes the map's stated resolution unverifiable.
#
# "era5"      0.25 deg, all three variables available. Safe default.
# "era5_land" 0.1 deg, finer, but CONFIRM relative_humidity_2m is served before
#             switching, and re-state the map resolution if you do.
#
# The driver's native resolution is the real floor on map resolution. At "era5"
# a 0.3 deg lattice is already at that floor and finer grids buy interpolation
# rather than information.
OPENMETEO_MODEL <- "era5"

# Locations per request. Open-Meteo accepts comma separated coordinate lists.
# Weighted cost is unchanged (it scales with locations) but round trips fall by
# a factor of OM_BATCH_SIZE, which is what actually governs wall clock. 25 keeps
# each request to roughly a minute of weighted budget and the URL well short of
# any length limit.
OM_BATCH_SIZE <- 25L

# Per REQUEST timeout, not per point. A 25 location, 61 day hourly request
# returns a few MB, so allow room, but not so much that a stalled connection can
# hold the run. Failures now retry with jittered backoff instead of being lost.
OM_TIMEOUT_S <- 60L
OM_MAX_ATTEMPTS <- 3L

# Free tier ceilings, all in WEIGHTED CALLS.
FREE_DAILY_CALLS   <- 10000
FREE_HOURLY_CALLS  <- 5000
FREE_MINUTELY_CALLS <- 600

################################################################################
# 5. Heatmap grid (for run_blast_grid.R)
################################################################################
# The heatmap runs the models on a grid of points and renders the result as a
# continuous risk surface, in the spirit of Savary et al. (2012), who mapped
# POTENTIAL epidemics. It colours the whole extent as if rice were grown
# everywhere, so most of the continent is a potential signal, not actual crop.

# Extent: c(lon_min, lon_max, lat_min, lat_max). Default = Australian continent.
GRID_EXTENT <- c(112, 154, -44, -10)

# Resolution. The lattice fills coarse to fine over successive runs, so any
# partial cache is still a uniform sample of the whole continent rather than a
# geographically biased block.
#
# WHAT SETS THE CEILING. Under GRID_WINDOW_MODE = "latest" every mapped cell must
# reach the same end_date, so the whole grid must be refreshed inside one day's
# 10,000 weighted calls. A refresh costs 1 weighted call per point (see the 14
# day floor below), so the same-date ceiling is about 10,000 points, which over
# the Australian land mass is about 0.26 deg. 0.3 deg is therefore close to the
# practical limit for this mode.
#
# To go finer, switch GRID_WINDOW_MODE to "coverage" and let the refresh spread
# across the week. Every cell is still truncated to one common end date, so the
# map is not a mosaic, but that common date sits up to a week further back.
# The weekly ceiling is then about 70,000 points, or roughly 0.1 deg, which is
# also ERA5-Land's native resolution. 0.15 deg (about 31,000 points, 44% of the
# weekly quota) is the sensible target with headroom.
GRID_RES_FINEST <- 0.3              # finest lattice spacing
GRID_RES_LEVELS <- c(1.2, 0.6, 0.3) # coarse-to-fine fill order

# Point cap per run (unit: FETCHES, i.e. points).
GRID_MAX_FETCHES_PER_RUN <- 8500L
# Refresh points stale by at least this many days. 0 brings every point to
# end_date each run, which is what "latest" window mode needs.
REFRESH_MIN_STALE_DAYS   <- 0L
# Budget per run (unit: WEIGHTED CALLS). Below the 10,000/day free ceiling.
DAILY_WEIGHTED_CAP       <- 9000
TARGET_CALLS_PER_RUN     <- DAILY_WEIGHTED_CAP

# Sustained request pacing (unit: WEIGHTED CALLS PER MINUTE).
# The hourly ceiling (5,000/h = 83/min) binds long before the per minute one
# (600/min), so pace to the hourly figure with a little margin. Do NOT read this
# as fetches per minute: at cost ~4.4 per point, 80 weighted/min is about 18
# points per minute.
GRID_TARGET_PER_MIN <- 80           # 4,800/hour, under the 5,000/hour cap

# Refresh always requests a fixed tail rather than only the missing days. The API
# charges a MINIMUM of 14 days, so a 1 day top up and a 14 day one cost exactly
# the same. Taking the full 14 costs nothing, self heals gaps, and picks up the
# ERA5T revisions applied to the most recent days. It also gives every refresh
# point a single shared start date, so they all batch into one request shape.
REFRESH_TAIL_DAYS <- 14L

# Wall clock budgets (unit: MINUTES).
GRID_MAX_MINUTES     <- 200L   # total fetch budget, measured from run start
GRID_RESERVE_MINUTES <- 25L    # modelling, rendering, saving and committing
# Held back from the add phase so the retry pass is reachable. Previously the add
# phase ran to the full deadline and the retry could never start, which stranded
# 617 attempted points on the 2026-07-28 run.
GRID_RETRY_RESERVE_MINUTES <- 20L
# Share of the WEIGHTED budget held back for the retry pass. Without it the add
# phase plans right up to DAILY_WEIGHTED_CAP, the retry is handed a budget of
# zero and stops immediately, which is the quota-flavoured version of the
# wall-clock bug above.
GRID_RETRY_WEIGHT_FRAC <- 0.05
# The workflow timeout must exceed GRID_MAX_MINUTES + GRID_RESERVE_MINUTES.

# Share of the fetch budget reserved for adding new cells. Refresh runs first.
# Without a reserve the grid stalls: once refreshing every cached cell fills the
# budget there is never time to add more.
GRID_ADD_RESERVE_FRAC <- 0.2

# Window mode.
# "latest":   the map window always ends at the archive edge; cells that missed
#             the refresh are absent. Freshness fixed, cell count varies.
# "coverage": pull the window back so GRID_WINDOW_COVERAGE of cells are included;
#             nearly all cells mapped, but the window is staler. Required if you
#             want to go finer than about 0.26 deg (see the note above).
GRID_WINDOW_MODE     <- "latest"
GRID_WINDOW_COVERAGE <- 0.98

# Rendering. The IDW search radius is now derived from the achieved spacing
# rather than fixed: the old value of 6 was six DEGREES, about 660 km and roughly
# fifteen times the point spacing, which interpolated smoothly straight across
# unsampled regions so a gap was indistinguishable from data.
IDW_RADIUS_MULT <- 1.5   # search radius = this many times the mean land spacing
# ABSOLUTE CAP on that radius, in degrees. Without it the radius is self
# defeating: sparse coverage raises the mean spacing, which widens the radius,
# which lets the interpolation smear across exactly the gaps the mask exists to
# expose. 1.0 deg is about 110 km, a defensible smoothing scale for a weather
# driven index. Past that you are inventing rather than interpolating.
IDW_RADIUS_MAX  <- 1.0
IDW_POWER       <- 2
IDW_MAX_POINTS  <- 12
# Blank cells further than the search radius from any modelled point, so gaps in
# the lattice read as gaps.
MASK_UNCOVERED  <- TRUE

# ---- Weather cache ---------------------------------------------------------
WEATHER_CACHE_GZ   <- "weather_cache.csv.gz"   # primary: gzipped
WEATHER_CACHE_CSV  <- "weather_cache.csv"      # fallback: plain CSV
# Write saves the gz to a temp path, reads it back to verify, then renames, so a
# cancelled job cannot leave a truncated cache in place.
#
# KEEP_CSV is now FALSE. The plain CSV is 16.4 MB against 2.1 MB for the gz, and
# BOTH were being committed on every run, so each run added about 18.5 MB of
# permanently unprunable git history. The gz has verified cleanly for many runs
# and the reader still falls back to a CSV if one is present.
WEATHER_CACHE_KEEP_CSV <- FALSE

# Retain weather beyond the modelling window, so a past map can be reproduced and
# retrospective analysis does not need the whole grid fetched again. Weather
# already fetched costs nothing to re-acquire.
#
# SIZE. Measured on the real cache: 61 days x 4,592 points is 2.1 MB gzipped, so
# roughly 0.034 MB per day of history at the current point count, and about three
# times that once the 0.3 deg grid is full. 120 days is therefore about 4 MB per
# commit, which is liveable weekly. Do NOT raise this much further while the
# cache is committed to git: 400 days at a full grid would be roughly 40 MB per
# run, or 2 GB of unprunable history a year. For longer retention, move the cache
# to a release asset or an orphan branch and raise this then.
CACHE_KEEP_HISTORY <- TRUE
CACHE_HISTORY_DAYS <- 120L

# Cache schema version. BUMP THIS whenever a change alters the VALUES stored in
# the cache, not just its columns. On a mismatch the cache is discarded and
# rebuilt, because a partial refresh would otherwise leave the map mixing old and
# new definitions with no way to tell them apart.
#
#  1  original
#  2  BLASTAM night window moved from UTC to local solar time, daily aggregates
#     moved to local solar days, and unjudgeable nights now store NA rather than
#     0. Every infect/semi/wet_hours/temp_wet value written under version 1 is
#     wrong and cannot be reused.
CACHE_SCHEMA_VERSION <- 2L
CACHE_VERSION_FILE   <- "cache_version.txt"

# Ledger of points whose fetch failed, so repeat failures can be deprioritised
# rather than retried in the same order every run.
FAIL_LEDGER_FILE <- "fetch_failures.csv"
FAIL_LEDGER_MAX_STRIKES <- 4L   # skip a point for a while after this many

# Drop ocean cells using the bundled Australia land polygon.
LAND_ONLY <- TRUE

# Also write a GeoTIFF of the risk surface alongside the PNG.
WRITE_GEOTIFF <- TRUE

# Heatmap colour ramps, low -> high.
HEAT_COLOURS         <- c("#BFE0F5", "#FFF6B0", "#FDB147", "#E8492B", "#8B0000")
BLASTAM_HEAT_COLOURS <- c("#BFE0F5", "#FFF6B0", "#FDB147", "#E8492B", "#8B0000")

# Heatmap colour scale maximum (%). A FIXED number makes every week's colours
# directly comparable. 2 means 2% intensity and above is the deepest red.
HEAT_MAX <- 2

# BLASTAM reporting window: count FAVOURABLE infection days over the most recent
# BLASTAM_WINDOW_DAYS, to show where infection pressure is building now.
BLASTAM_WINDOW_DAYS <- 21
# End the window this many days before the data edge. A night needs the FOLLOWING
# morning's hours, and how much of that morning falls inside a UTC fetch depends
# on longitude, so without this the westernmost cells lose their final night and
# the map gains a spurious east-west gradient of about 1/21. One day makes the
# window longitude-neutral across the continent.
BLASTAM_END_LAG_DAYS <- 1L
# Deepest colour at this many favourable days within the window. PROVISIONAL:
# lower it (e.g. 14) to make building pressure show up more strongly.
BLASTAM_HEAT_MAX <- 21

# Overlay layers on the heatmap
SHOW_COAST  <- TRUE
SHOW_TOWNS  <- TRUE
SHOW_ROADS  <- TRUE
SHOW_RIVERS <- TRUE

COL_COAST <- "#22272B"
COL_ROAD  <- "#8A6D3B"
COL_RIVER <- "#2E75B6"
COL_TOWN  <- "#111111"

################################################################################
# 6. Monitoring towns
################################################################################
# Modelled each run, tracked in the trends CSVs, and highlighted on the map.
# Coordinates for the remote roadhouses (Archer River, Lakeland, Timber Creek)
# are approximate; edit any as needed.
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

# Number of recent runs to keep in the trends CSVs (one column per run).
HISTORY_RUNS <- 10

# Town fetches are batched into a single request now, so no fork pool is needed.
# Retained only for the serial fallback path.
TOWN_FETCH_CORES <- 4L

################################################################################
# 7. Citation (added to the summary and to the heatmap footer)
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
  "Weather: Open-Meteo ERA5 archive (data CC BY 4.0), non-commercial research use.",
  sep = "\n"
)

# Compact credit for the map footer
CITATION_MAP <- paste(
  "Model: EPIRICE (Savary et al. 2012, Crop Prot. 34:6-17; doi:10.1016/j.cropro.2011.11.009).",
  "Code: epicrop (A.H. Sparks). Weather: Open-Meteo ERA5 (CC BY 4.0). Lineage: BLASTL (Hashimoto et al. 1984).",
  sep = "\n"
)

################################################################################
# 8. Removed settings
################################################################################
# The following were present but unused, and are removed to stop them misleading
# a future reader:
#   SITES            superseded by MONITOR_TOWNS, which is what run_blast.R reads
#   EMERGENCE_DATE   superseded by the rolling (end_date - CROP_AGE_DAYS)
#   GRID_RES         legacy fixed resolution, unused by the dynamic runner
#   SMOOTH_FACTOR    never read
#   EMAIL_*, SEND_EMAIL, mailR path: the workflow sends mail via send_email.py
#   GRID_CONC, GRID_CONC_MIN, GRID_CONC_MAX, GRID_FAIL_BACKOFF, GRID_PROBE_CHUNKS
#                    the adaptive concurrency controller is gone. Its pacer made
#                    throughput algebraically independent of concurrency, so it
#                    was hill climbing on noise while the failure rate rose.
#   GRID_BATCHES, GRID_BATCH_CALLS, GRID_BATCH_WAIT_S
#                    the hour-long sleep machinery, unnecessary with token bucket
#                    pacing.
#   GRID_FETCH_TIMEOUT_S  replaced by OM_TIMEOUT_S (per request, not per point)

# Below this success fraction the retry pass is skipped: a wholesale failure is a
# bad day at the archive, not transient flakiness, and retrying doubles the waste.
GRID_RETRY_MIN_OK_FRAC <- 0.25

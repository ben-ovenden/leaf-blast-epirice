################################################################################
# blast_config.R
#
# Edit this file to configure the weekly leaf blast run, then commit to git.
# This file is the SINGLE SOURCE OF TRUTH for every tunable parameter. It is
# sourced BEFORE the model files, and the model files now guard their own
# defaults with if (!exists(...)), so a value set here always wins.
#
# UNITS WARNING. Two different quantities are easy to confuse and confusing them
# was the cause of several past bugs:
#   * a FETCH is one grid point brought up to date;
#   * a WEIGHTED CALL is what Open-Meteo actually charges.
# They are not the same. The 68 day fetch of three hourly variables used when a
# new point is added costs about 4.86 weighted calls, not 1. Every budget below
# is labelled with its unit.
################################################################################

################################################################################
# 0. Run identity
################################################################################
# ONE run date for the whole pipeline. A grid run takes hours, so calling
# Sys.Date() at several points let a run straddle local midnight and produce a
# map dated one day and a town table dated the next. The 2026-07-30 email showed
# exactly that: the maps said "weather to 2026-07-23" while the body said
# "weather to 24 Jul 2026".
#
# The workflow exports BLAST_RUN_DATE once and both scripts read it here, so the
# two products cannot disagree. Locally, leave it unset and today's date is used.
blast_run_date <- function() {
  v <- Sys.getenv("BLAST_RUN_DATE", "")
  d <- if (nzchar(v)) suppressWarnings(as.Date(v)) else as.Date(NA)
  if (is.na(d)) Sys.Date() else d
}

# The UTC calendar date the job started on, pinned once by the workflow as
# BLAST_UTC_DATE for the same reason as the run date: a grid run lasts long
# enough to cross 00:00 UTC, and the town script must not see a later date than
# the grid script did.
blast_utc_date <- function() {
  v <- Sys.getenv("BLAST_UTC_DATE", "")
  d <- if (nzchar(v)) suppressWarnings(as.Date(v)) else as.Date(NA)
  if (is.na(d)) as.Date(format(Sys.time(), "%Y-%m-%d", tz = "UTC")) else d
}

# The last day FETCHED. The archive's availability is a function of the UTC
# clock, not of the Sydney calendar, so the lag is counted from whichever of the
# two dates is EARLIER. See ARCHIVE_LAG_DAYS in section 2 for what went wrong
# when it was counted from the Sydney date alone.
blast_archive_edge <- function(run_date = blast_run_date()) {
  min(as.Date(run_date), blast_utc_date()) - ARCHIVE_LAG_DAYS
}

# The arithmetic above is an UPPER bound. Open-Meteo publishes each new ERA5 day
# at about 00:30 to 01:00 UTC, so early in a UTC day the arithmetic edge is not
# there yet. The workflows run probe_archive_edge.R, which asks the archive, and
# export its answer as BLAST_DATA_END; it may only pull the window BACK. Unset
# (a local run), the arithmetic edge is used, as before.
blast_data_end <- function(run_date = blast_run_date()) {
  edge <- blast_archive_edge(run_date)
  v <- Sys.getenv("BLAST_DATA_END", "")
  d <- if (nzchar(v)) tryCatch(suppressWarnings(as.Date(v)), error = function(e) as.Date(NA)) else as.Date(NA)
  if (!is.na(d) && d <= edge) d else edge
}

################################################################################
# 1. Season and risk bands
################################################################################
# EPIRICE reports leaf blast intensity as a proportion of diseased sites (0-1).
# 1% (0.01) intensity is taken as the threshold for a significant epidemic, so it
# is the BOTTOM of the "high" band: anything at or above 1% is high, and the band
# is open ended above. Calibrate against field observations.
INTENSITY_LOW_MAX      <- 0.002  # below 0.2%          = low
INTENSITY_MODERATE_MAX <- 0.01   # 0.2% to below 1%    = moderate
                                 # 1% and above        = high (no upper bound)

# Crop age framing: each run, assume a crop this many days old everywhere and
# run the model over the trailing weather window. Makes the surface spatially
# comparable and bounds the API cost.
#
# The modelled series runs emergence to end_date inclusive, which is
# CROP_AGE_DAYS + 1 rows, so the FINAL day carries a crop age of CROP_AGE_DAYS.
#
# NOTE for the trends CSVs: because emergence is (end_date - CROP_AGE_DAYS), it
# moves with every run. The trends columns are therefore a ROLLING 60 day window,
# not a cumulative season total, and must not be read as season progress.
CROP_AGE_DAYS <- 60L

# Minimum days of weather after emergence before a run is meaningful. Disease
# onset in the leaf blast model is day 15. With the rolling emergence above this
# never binds, but it still guards a fixed per site emergence override.
MIN_DAYS <- 16

################################################################################
# 2. Data window
################################################################################
# The Open-Meteo archive (ERA5) lags real time by about 5 days. End the fetch a
# few days back so the most recent days are actually populated.
#   data_end  = min(run date, UTC date) - ARCHIVE_LAG_DAYS   the last day FETCHED
#   end_date  = data_end - DAY_CUT_LAG_DAYS                  the last day MODELLED
#
# WHY THE UTC DATE IS IN THERE. The Monday job fires at about 22:30 UTC on
# SUNDAY, so the Sydney run date is a day ahead of the UTC date and "run date
# minus 6" is only five days behind the archive's clock. That day was not in the
# archive yet on any 22:30 UTC run through August and September 2026. With the
# model day cut at 10:00 local solar, a cell at longitude L needs (10 - L/15)
# hours of that missing day to complete end_date; the completeness rule allows
# two, so every cell west of 120 E landed one day short: 841 cells, 10.9% of the
# grid, hence "only 89% of cached cells reached end_date" on every scheduled run
# and a "Degraded run" banner on every Monday email. A week later those cells
# were more than a tail behind and cost a full refetch each. Counting the lag
# from the UTC date costs one day of freshness on the Monday email and removes
# the whole chain. Manual runs after 10:00 Sydney are unaffected.
ARCHIVE_LAG_DAYS <- 6L

# WHEN THE NEW DAY APPEARS. "UTC date minus 6" is published at about 00:30 to
# 01:00 UTC, not at midnight: it was still missing at 00:08 and 00:11 UTC on
# 2 October 2026, a run that began at 00:26 UTC on 25 August caught it a few
# minutes in (38 western cells a day short, 649 not), and every run that started
# at 05:09 UTC or later had it. So the workflows ask (probe_archive_edge.R). The
# probe requests one far-western point, because that is where the last day's
# hours are needed, and the weekly run waits for the update because the extra day
# is the reason it runs on Monday UTC at all.
ARCHIVE_EDGE_PROBE_LON    <- 114.0
ARCHIVE_EDGE_PROBE_LAT    <- -26.0
ARCHIVE_EDGE_WAIT_MAX_MIN <- 90    # weekly run: wait at most this long for the new day
ARCHIVE_EDGE_POLL_MIN     <- 10    # and ask again this often (1 weighted call each)

# Why the extra day. The model day now runs from BLASTAM_DAY_CUT_HOUR local solar
# (see section 4a), so the last fetched day is only partially covered, and how
# much of it is covered depends on longitude. Dropping one day makes the modelled
# window identical at every longitude by construction, rather than correcting for
# it afterwards. This replaces the old BLASTAM_END_LAG_DAYS fudge, which is now 0.
DAY_CUT_LAG_DAYS <- 1L

################################################################################
# 3. Output
################################################################################
OUTPUT_DIR <- "blast_outputs"
MAKE_MAP   <- FALSE         # off: the workflow emails the two grid heatmaps;
                            # this simple town point map is not delivered
MAP_WIDTH  <- 900
MAP_HEIGHT <- 800

# NSW Government brand palette (https://designsystem.nsw.gov.au/core/colour/)
NSW_BRAND_BLUE <- "#002664"
NSW_BRAND_DARK <- "#22272B"
NSW_GREY_01    <- "#F4F7F9"
NSW_GREY_02    <- "#DEE3E5"
NSW_GREY_04    <- "#495054"
COL_LOW        <- "#00AA45"   # NSW success green
COL_MODERATE   <- "#DC5800"   # NSW warning orange
COL_HIGH       <- "#B81237"   # NSW error red
COL_NODATA     <- "#C8CDD0"   # NSW grey

# Band styling for the HTML email. Tints with dark text rather than white text on
# orange, which does not reach 4.5:1 contrast at this font size. The full
# strength brand colour appears as a dot beside the label.
BAND_STYLE <- list(
  low        = list(bg = "#E5F6EC", fg = NSW_BRAND_DARK, dot = COL_LOW),
  moderate   = list(bg = "#FDEBD9", fg = "#8F3B00",      dot = COL_MODERATE),
  high       = list(bg = "#FBE3E9", fg = "#8B0E29",      dot = COL_HIGH),
  `no data`  = list(bg = "#EBEDEF", fg = NSW_GREY_04,    dot = COL_NODATA)
)
# Single hue blue ramp for the BLASTAM day count column, low to high.
BLASTAM_CELL_BG <- c("#F4F7F9", "#EBF1F8", "#C7DBEF", "#93B8DE", "#4B7EBD", NSW_BRAND_BLUE)

CONTACT_NAME  <- "Ben Ovenden"
CONTACT_EMAIL <- "ben.ovenden@dpird.nsw.gov.au"
SENDER_LABEL  <- "WWAI Cereal Pathology: blast models"

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
OPENMETEO_MODEL <- "era5"

# Locations per request. Weighted cost is unchanged (it scales with locations)
# but round trips fall by a factor of OM_BATCH_SIZE, which is what actually
# governs wall clock.
OM_BATCH_SIZE <- 25L

# Per REQUEST timeout, not per point.
OM_TIMEOUT_S <- 60L
OM_MAX_ATTEMPTS <- 3L

# HTTP 429 handling. Open-Meteo enforces three ceilings (per minute, per hour,
# per day) and says which one in the response body. Only the DAILY one means the
# run is over. The 2026-09-07 grid run hit the HOURLY ceiling at exactly 5,000
# weighted calls, about 60 minutes in, and the code then treated every 429 as
# "quota spent": it abandoned the remaining 42% of the grid with 140 minutes of
# wall clock still in hand, and the town fetch that followed in the same hour was
# refused too, so that email carried 31 towns of "no data". A minutely or hourly
# 429 is now waited out and the fetch resumes; the wait is bounded by the total
# below and by the run's own deadline, and a 429 is never charged.
OM_QUOTA_WAIT_S        <- 300   # seconds between probes after an hourly 429
OM_QUOTA_WAIT_MAX_MIN  <- 60    # total 429 waiting a GRID run will tolerate
TOWN_QUOTA_WAIT_MAX_MIN <- 10   # and a TOWN run, which has no fetch deadline of
                                # its own and must leave time to commit and email

# Free tier ceilings, all in WEIGHTED CALLS. Reference values, used by the shared
# spend ledger below to keep the grid run and the town run inside one day's quota.
FREE_DAILY_CALLS    <- 10000
FREE_HOURLY_CALLS   <- 5000
FREE_MINUTELY_CALLS <- 600

# Combined ceiling across every script in one UTC day. run_blast_grid.R plans to
# DAILY_WEIGHTED_CAP, then run_blast.R reads the ledger and caps its own town
# fetch at whatever is left below this figure. Without the ledger the town run
# spent its ~150 weighted calls unbudgeted on top of the grid's 9,000.
DAILY_WEIGHTED_HARD_CAP <- 9500
SPEND_LEDGER_FILE       <- "weighted_spend.csv"

################################################################################
# 4a. BLASTAM parameters (single source of truth)
################################################################################
# blastam_model.R reads these. They used to live only in the model file, which
# meant the README's advice to set them here silently did nothing: the model file
# is sourced second and overwrote them.

# Temperature bounds. Koshimizu 1988 originals in brackets. Deviation 2: the
# upper bounds are raised for tropical Australia; the lower bounds are not.
BLASTAM_TWET_MIN  <- 15    # wetness-period mean lower bound (C)  [Japan: 15]
BLASTAM_TWET_MAX  <- 32    # wetness-period mean upper bound (C)  [Japan: 25]
BLASTAM_PREV5_MIN <- 20    # preceding 5-day mean lower bound (C) [Japan: 20]
BLASTAM_PREV5_MAX <- 30    # preceding 5-day mean upper bound (C) [Japan: 25]

# Leaf wetness proxy (Deviation 4)
BLASTAM_RH_WET     <- 90    # RH (%) at or above this = leaf wet
BLASTAM_RAIN_WET   <- 0.2   # hourly rain (mm) at or above this = leaf wet
BLASTAM_RAIN_HEAVY <- 4.0   # hourly rain (mm) at or above this = spores washed
                            # off, hour excluded from the infection-wet count
                            # (Yoshino 1988). Set to Inf to disable.

# Wetness duration threshold (Deviation 1)
BLASTAM_USE_BJ_THRESHOLD <- TRUE   # TRUE: Barksdale & Jones 1965 curve
BLASTAM_WET_HOURS_FIXED  <- 10L    # used when the above is FALSE

# Night window, in LOCAL SOLAR hours (Deviation 5)
BLASTAM_NIGHT_START <- 15L
BLASTAM_NIGHT_END   <- 9L

# Where the 24 hour MODEL DAY starts, in local solar hours.
#
# WHY THIS IS NOT MIDNIGHT. Schema 2 moved both the BLASTAM night window and the
# EPIRICE daily aggregates onto local solar days cut at midnight. The night
# window fix was right; moving the daily aggregates with it was not. EPIRICE's
# rainlim gate is a daily SUM, so cutting at local midnight splits a nocturnal
# rain event across two days and halves the peak daily total. On a synthetic
# series with 6 mm falling between 22:00 and 03:00 local every third night, the
# number of days reaching the 5 mm gate went from 24 to 0 and final intensity
# from 0.0616% to 0.0000%, on identical physical rainfall. The measured symptom
# was Malanda dropping from 0.374% to 0.006% between the 2026-07-28 and
# 2026-07-29 runs, with every other tropical town falling to zero at the same
# step while the southern NSW sites did not move.
#
# Cutting at 10:00 local solar puts the whole nocturnal wet and rain period
# inside one model day, which is what the UTC day did by accident before schema 2
# and what a rain day does by convention. A model day labelled 23 July runs from
# 10:00 on 23 July to 09:59 on 24 July, local solar.
BLASTAM_DAY_CUT_HOUR <- 10L

# Lead-in days the caller must fetch and discard: 1 for the local solar shift
# (the first model day of any fetch is partial) plus 5 for the lagged preceding
# 5-day mean temperature.
BLASTAM_LEADIN_DAYS <- 6L

# Completeness. A night used to be judged on a single evening hour plus a single
# morning hour, and hours with missing humidity counted as dry, so a mostly empty
# night was scored "not favourable" rather than "not judged".
BLASTAM_MIN_EVE_HOURS  <- 7L    # of the 9 hours from 15:00 to 23:59
BLASTAM_MIN_MORN_HOURS <- 7L    # of the 9 hours from 00:00 to 08:59
BLASTAM_MAX_NA_FRAC    <- 0.10  # more missing than this and the night is NA

# BLASTAM reporting window: count FAVOURABLE infection days over the most recent
# BLASTAM_WINDOW_DAYS, to show where infection pressure is building now.
BLASTAM_WINDOW_DAYS <- 21L
BLASTAM_RECENT_DAYS <- 7L
# Now 0. DAY_CUT_LAG_DAYS in section 2 removes the partial final day for every
# longitude at once, so there is nothing left for this to correct. Kept as a knob.
BLASTAM_END_LAG_DAYS <- 0L

################################################################################
# 4b. EPIRICE parameters
################################################################################
# Optimum temperature for the RcT infection-rate curve, in C. Either 25 or 20.
#
#  25  The value published in Table 2 of Savary et al. 2012. This is the default.
#      Hashioka (1965), the data behind the Lanoiselet et al. (2002) DYMEX model
#      for the NSW rice belt, measured minimum germination and penetration times
#      of 6 h at 24 C, 8 h at 28 C and 10 h at 32 C, which puts the infection
#      rate optimum at 24 to 25 C. The broader literature reports 25 to 28 C for
#      germination, infection, lesion formation and sporulation.
#
#  20  The epicrop package value. Its source carries a comment asserting that the
#      published 25 C is a typo. That is a package author's reading of the paper,
#      not a published erratum.
#
# THE CHOICE IS LARGE. Holding RH at 92% so the wetness gate is open every day,
# final intensity over 61 days differs by roughly 6 to 20 fold at 26 to 30 C:
#
#     TEMP C        20   22     24     26     28     30
#     peak 20    4.73%  2.65%  1.34%  0.59%  0.21%  0.05%
#     peak 25    0.91%  1.91%  3.58%  3.58%  1.91%  0.91%
#
# Below about 22 C the 20 C curve gives the higher answer, so this is NOT a
# uniform scaling and the two are not interchangeable. Note also that sporulation
# peaks cooler than infection (Kato and Kozaka 1974, Phytopathology 64: 828-830:
# 399 spores per lesion per day at 20 C, 271 at 25 C, 131 at 32 C) and EPIRICE
# uses one RcT for both, so
# either choice is a compromise. 25 favours the infection step, which is the
# mechanistic basis for the SEIR transition and the reason for the default.
EPIRICE_RCT_PEAK <- 25L

################################################################################
# 5. Heatmap grid (for run_blast_grid.R)
################################################################################
# Extent: c(lon_min, lon_max, lat_min, lat_max). Default = Australian continent.
# build_targets() rounds the extent OUT to a whole number of cells, so the
# northern row is no longer dropped: seq(-44, -10, by = 0.3) stopped at -10.1.
GRID_EXTENT <- c(112, 154, -44, -10)

# Resolution. The lattice fills coarse to fine over successive runs, so any
# partial cache is still a uniform sample of the whole continent.
#
# WHAT SETS THE CEILING. Under GRID_WINDOW_MODE = "latest" every mapped cell must
# reach the same end_date, so the whole grid must be refreshed inside one day's
# 10,000 weighted calls. A refresh costs 1 weighted call per point (see the 14
# day floor below), so the same-date ceiling is about 10,000 points, which over
# the Australian land mass is about 0.26 deg. 0.3 deg is close to the practical
# limit for this mode.
#
# To go finer, switch GRID_WINDOW_MODE to "coverage". That path now works: it
# used to hand SEIR the run's global emergence date while truncating the weather
# to an earlier model_end, so the alignment check threw for every point and the
# EPIRICE map came out empty while the BLASTAM map still rendered.
GRID_RES_FINEST <- 0.3              # finest lattice spacing
GRID_RES_LEVELS <- c(1.2, 0.6, 0.3) # coarse-to-fine fill order

# Point cap per run (unit: FETCHES, i.e. points).
GRID_MAX_FETCHES_PER_RUN <- 8500L
# Refresh points stale by at least this many days. 0 brings every point to
# end_date each run, which is what "latest" window mode needs.
REFRESH_MIN_STALE_DAYS   <- 0L
# Budget per run (unit: WEIGHTED CALLS). Below the 10,000/day free ceiling.
DAILY_WEIGHTED_CAP       <- 9000
# BLAST_WEIGHTED_CAP can only LOWER it. The top-up workflow's manual dispatch
# passes it, so a change to the workflows can be proved end to end (restore the
# cache, probe, fetch a batch, save, commit) in three minutes and thirty calls
# instead of two hours and a day's quota.
.cap_env <- suppressWarnings(as.numeric(Sys.getenv("BLAST_WEIGHTED_CAP", "")))
if (isTRUE(is.finite(.cap_env)) && .cap_env >= 0)
  DAILY_WEIGHTED_CAP <- min(DAILY_WEIGHTED_CAP, .cap_env)

# Sustained request pacing (unit: WEIGHTED CALLS PER MINUTE).
# The hourly ceiling (5,000/h = 83/min) binds long before the per minute one
# (600/min). Do NOT read this as fetches per minute: at ~4.86 weighted per new
# point, 70 weighted/min is about 14 points per minute.
#
# This was 80, which is 96% of the hourly ceiling with nothing in hand for the
# pacer's opening burst, for retried requests (which were not paced at all) or
# for the ~150 weighted the town run spends in the grid run's final hour. The
# 2026-09-07 run crossed the line at exactly 5,000. At 70 an hour of grid
# fetching plus the town run comes to about 4,400.
GRID_TARGET_PER_MIN <- 70           # 4,200/hour, under the 5,000/hour cap

# Every cached cell is refreshed from the day after its newest cached row, plus
# BLASTAM_LEADIN_DAYS of lead-in and the day cut lag, so it pays for exactly the
# days it is missing. The API charges a MINIMUM of 14 days, so a cell refreshed on
# the weekly cadence (7 days missing) costs exactly 1.00 weighted:
#
# ** REFRESH_TAIL_DAYS + BLASTAM_LEADIN_DAYS + DAY_CUT_LAG_DAYS MUST COME TO 14. **
# 7 + 6 + 1 = 14. At 15 days the weight is 15/14 = 1.07, and a 7% surcharge on
# every weekly refresh costs about 540 weighted calls once the grid is full. The
# runner checks this sum and warns. REFRESH_TAIL_DAYS is the weekly cadence that
# the check describes; the planner itself no longer uses a fixed tail.
#
# Cells further behind are grouped into cohorts by how far behind they are and
# each cohort is fetched over its own window: 8 days missing costs 1.07, 15 days
# 1.57, 22 days 2.07, and never more than a brand new cell (4.86). Previously any
# cell more than 7 days behind was refetched over the FULL crop window at 4.86,
# which is why one interrupted run in September 2026 took the grid three weeks
# to recover from instead of one.
REFRESH_TAIL_DAYS <- 7L

# Wall clock budgets (unit: MINUTES). The workflow's timeout-minutes must exceed
# GRID_MAX_MINUTES + GRID_RESERVE_MINUTES + TOWN_QUOTA_WAIT_MAX_MIN plus a margin
# for the town fetch, commit and email.
GRID_MAX_MINUTES     <- 200L   # total fetch budget, measured from run start
GRID_RESERVE_MINUTES <- 25L    # modelling, rendering, saving and committing
GRID_RETRY_RESERVE_MINUTES <- 20L
# Share of the WEIGHTED budget held back for the retry pass. This is now enforced
# on the fetch itself, not only in planning: the 2026-07-30 run reported 8,667
# weighted spent against a planned 8,550, because charged retries ate the slice.
GRID_RETRY_WEIGHT_FRAC <- 0.05

# Share of the fetch budget reserved for adding new cells. Refresh runs first.
GRID_ADD_RESERVE_FRAC <- 0.2

# Window mode. "latest": the map window always ends at the archive edge; cells
# that missed the refresh are absent. "coverage": pull the window back so
# GRID_WINDOW_COVERAGE of cells are included, at the cost of a staler window.
GRID_WINDOW_MODE     <- "latest"
GRID_WINDOW_COVERAGE <- 0.98

# GRACEFUL DEGRADATION FOR "latest".
#
# "latest" pins the window to end_date and drops any cell that has not reached
# it. That is right when the run refreshed normally, and catastrophic when it
# could not: if nothing was refreshed, NO cell reaches end_date, nothing is
# modelled, no map renders and the email step fails on a missing attachment. The
# 2026-07-31 run hit exactly this after the weighted ledger correctly capped it
# at 4 calls following an earlier run the same UTC day. The cache held 1,874
# perfectly good points and the run produced nothing from them.
#
# So: use end_date only if at least this fraction of cached cells reach it.
# Otherwise step back to the newest date that this fraction DO reach. The
# same-date guarantee is preserved either way; only the date moves.
GRID_WINDOW_MIN_COVERAGE <- 0.90
# Warn loudly (do not cap) if the fallback lands more than this far back. A very
# stale map that says so is more useful than no map, but it should be obvious.
GRID_WINDOW_WARN_FALLBACK_DAYS <- 10L
# ...but the fallback is no longer unconditional. When fewer than
# GRID_WINDOW_MIN_COVERAGE reach end_date, the window may step back at most this
# many days to bring every cell in: a complete map a day or two old beats a
# current one with holes. Further back than that, and provided at least
# GRID_WINDOW_MIN_DRAW_COVERAGE of cells are current, the map is drawn at
# end_date and the stale cells are drawn GREY, with their count and age in the
# footer and the email. Below that fraction the grid is mostly stale and the old
# complete map is the more useful thing. See grid_window.R.
#
# In September 2026 the unlimited fallback drew every cell at a 16 day old window
# three weeks running while 76% of the grid was current.
GRID_WINDOW_MAX_FALLBACK_DAYS <- 3L
GRID_WINDOW_MIN_DRAW_COVERAGE <- 0.50
COL_STALE <- COL_NODATA        # cells not refreshed to the window are drawn in this
# Under "coverage" the common window ends earlier than end_date, so the EPIRICE
# window STARTS earlier too, and the cache must already hold weather before the
# current run's emergence date. New points are therefore fetched with this many
# extra days of lookback in coverage mode, at a cost of about 0.5 weighted calls
# each. Without it the first coverage run reports NA intensity for every point,
# for a reason unrelated to the alignment bug this replaced.
GRID_WINDOW_MAX_LAG_DAYS <- 7L

# Grid top-ups (.github/workflows/midweek_topup.yml): two extra fetch days a
# week, Thursday and Sunday UTC, run with BLAST_MIDWEEK=1, which fetches, saves
# the cache and stops. Thursday recovers from a bad Monday and speeds a cold
# fill; Sunday is Monday's insurance, leaving unreached cells one day behind
# instead of three. run_blast.R reports the last top-up in the Monday email and
# complains when it is older than this, which on a Monday means Sunday's did not
# run.
MIDWEEK_MAX_AGE_DAYS <- 2L

# Rendering. The IDW search radius is derived from the achieved spacing.
IDW_RADIUS_MULT <- 1.5   # search radius = this many times the mean land spacing
IDW_RADIUS_MAX  <- 1.0   # absolute cap, degrees (~110 km)
IDW_POWER       <- 2
IDW_MAX_POINTS  <- 12
MASK_UNCOVERED  <- TRUE

# Coastal cells. ERA5 cells on the coastal fringe are partly marine, so their
# humidity is not representative of any paddock, yet they carry most of the
# BLASTAM signal on the delivered maps and pull the eye to country where rice is
# not grown. Set a distance in km to blank cells that close to the coastline.
# 0 disables. Applied as a raster mask at render time only; the cells are still
# fetched, cached and written to the GeoTIFF.
COAST_MASK_KM <- 0

# Overlay layers. australia_roads.geojson contains at least one feature with a
# 3.25 deg jump from Victoria straight to Tasmania, which drew as a line across
# Bass Strait on every map. Line parts are split at jumps longer than this.
OVERLAY_MAX_SEGMENT_DEG <- 1.0

SHOW_COAST  <- TRUE
SHOW_TOWNS  <- TRUE
SHOW_ROADS  <- TRUE
SHOW_RIVERS <- TRUE
# Town labels. Each one tries twelve positions around its marker and takes the
# first that stays on the plot, touches no other label and covers no town marker
# (map_labels.R). The earlier scheme fixed every label to one side and dropped any
# whose MARKER was within 0.9 degrees of another: "Humpty Doo" still ran through
# Jabiru and "Kununurra" through Timber Creek, 1.6 and 1.7 degrees away on the
# same latitude, while five towns that had room on another side went unlabelled.
LABEL_CEX          <- 0.5
# White halo behind each name, in inches (about one pixel at 120 dpi), so a label
# stays readable where it crosses the coast, a road or a river. 0 turns it off.
LABEL_HALO_IN      <- 0.010
# Extra longitude drawn east of GRID_EXTENT, all of it sea, so that labels for
# the east coast towns have somewhere to go. Display only: nothing is fetched or
# modelled there.
MAP_EAST_PAD_DEG   <- 2

COL_COAST <- NSW_BRAND_DARK
COL_ROAD  <- "#8A6D3B"
COL_RIVER <- "#2E75B6"
COL_TOWN  <- "#111111"

# ---- Weather cache ---------------------------------------------------------
WEATHER_CACHE_GZ   <- "weather_cache.csv.gz"
WEATHER_CACHE_CSV  <- "weather_cache.csv"
WEATHER_CACHE_KEEP_CSV <- FALSE
# How far back the cache reaches. The models need CROP_AGE_DAYS + 1 = 61 days (68
# in coverage mode); 90 leaves room for a fallback window three weeks back and
# nothing reads anything older. It was 120. The outputs that matter over a season
# are in the trends CSVs, not here.
#
# This bounds the FILE. What bounds the REPOSITORY is that the cache is no longer
# committed to main at all (see "Restore weather cache" in the workflows): git
# keeps every version of a committed file, however short each one is.
CACHE_KEEP_HISTORY <- TRUE
CACHE_HISTORY_DAYS <- 90L

# Cache schema version. BUMP THIS whenever a change alters the VALUES stored in
# the cache, not just its columns. On a mismatch the cache is discarded and
# rebuilt, because a partial refresh would otherwise leave the map mixing old and
# new definitions with no way to tell them apart.
#
#  1  original
#  2  BLASTAM night window moved from UTC to local solar time, daily aggregates
#     moved to local solar midnight days, unjudgeable nights store NA not 0.
#  3  Three value changes, all of which invalidate version 2 rows:
#       * the model day is cut at BLASTAM_DAY_CUT_HOUR (10:00 local solar), not
#         local midnight, so every TEMP, RHUM and RAIN value changes;
#       * the preceding 5-day mean is now genuinely preceding (lagged one day),
#         so infect and semi change;
#       * night completeness now requires a minimum number of hours and treats
#         missing humidity as unknown rather than dry.
CACHE_SCHEMA_VERSION <- 3L
CACHE_VERSION_FILE   <- "cache_version.txt"

# Ledger of points whose fetch failed.
FAIL_LEDGER_FILE <- "fetch_failures.csv"
FAIL_LEDGER_MAX_STRIKES <- 4L

LAND_ONLY <- TRUE
WRITE_GEOTIFF <- TRUE

# Heatmap colour ramps, low -> high: near white for nothing, pale green for a
# little, then yellow, NSW warning orange and NSW error red. Green is what the
# email's risk table already uses for "low".
#
# The bottom of the ramp used to be a pale grey running into light blue, so the
# first step up from "nothing" was a COOLER colour, and a patch of low risk read
# as less than the grey around it. The yellow, orange and red are unchanged. Two
# alternatives were rendered on the 2026-09-28 data and rejected: white straight
# to yellow (low values all but vanish) and blue as the zero colour (low values
# become pale holes in a blue continent).
HEAT_COLOURS         <- c("#F6F6F4", "#CFE8C8", "#FFF6B0", COL_MODERATE, COL_HIGH)
BLASTAM_HEAT_COLOURS <- HEAT_COLOURS

# Heatmap colour scale maxima. FIXED numbers make every week's colours directly
# comparable. HEAT_MAX is EPIRICE intensity in %; BLASTAM_HEAT_MAX is days.
HEAT_MAX         <- 2
BLASTAM_HEAT_MAX <- 21

# Colour stretch exponent. 1 = linear. Below 1 expands the low end.
#
# THIS IS NOW WIRED UP. It was documented in the README as a working control but
# no script read it, so both delivered maps rendered as one flat pale blue:
# EPIRICE peaked near 0.03% against a 2% scale, using about 1.5% of the ramp.
# The anchor stays fixed (HEAT_MAX is still the deepest red every week) and the
# legend is labelled with TRUE values, so nothing is misrepresented; only the
# spacing of the colours changes. Set to 1 for linear.
HEAT_STRETCH    <- 0.4
BLASTAM_STRETCH <- 0.6
# Print the observed maximum in the map footer, so a reader can tell a genuinely
# flat map from a broken one.
SHOW_OBSERVED_MAX <- TRUE

################################################################################
# 6. Monitoring towns
################################################################################
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

# Number of recent runs to keep in the trends CSVs (one column per DATA end
# date, not per run date, so three test runs on consecutive days no longer
# consume three columns of history describing the same weather).
HISTORY_RUNS <- 12
# Sidecar recording one row per run: run date, data window, schema version and
# model options. Without it a methodological break in the trends series looks
# like a real epidemiological collapse.
RUN_LOG_FILE <- "run_log.csv"

################################################################################
# 6a. Run health: when is a run degraded? (run_health.R)
################################################################################
# One verdict per run, written to RUN_STATUS_FILE by run_blast.R. It prefixes the
# email subject with "[DEGRADED]" and turns the workflow run red in its final
# step, after everything has been committed, emailed and uploaded. Until this
# existed, a run that delivered 31 towns of "no data" (2026-09-07) or a map three
# weeks stale (09-14, 09-21) exited 0 and showed a green tick.
#
# These fail the run:
HEALTH_MIN_TOWN_FRAC       <- 0.9   # fewer towns modelled than this share
HEALTH_MAX_MAP_BEHIND_DAYS <- 1L    # map window more than this behind the towns
HEALTH_MAX_GREY_FRAC       <- 0.02  # more than this share of cells drawn grey
HEALTH_MAX_MAP_DROP_FRAC   <- 0.10  # mapped cells fell by more than this vs last run
# Anything the grid fetch recorded as its reason (a 429, the deadline, a budget
# shortfall) is carried as a WARNING, not a failure: it is the product that is
# judged, not the weather API's day.
RUN_STATUS_FILE <- "run_status.txt"

################################################################################
# 7. Citation (added to the summary and to the heatmap footer)
################################################################################
CITATION <- paste(
  "Models and data citation",
  "------------------------------------------------------------",
  "EPIRICE (intensity): Savary, S., Nelson, A., Willocquet, L., Pangga, I.,",
  "  and Aunario, J. (2012). Modeling and mapping potential epidemics of",
  "  rice diseases globally. Crop Protection 34: 6-17.",
  "  doi:10.1016/j.cropro.2011.11.009",
  "  Software: epicrop R package (SEIR engine), A.H. Sparks and colleagues,",
  "  adapted from cropsim (Hijmans et al. 2009). Framework: Zadoks (1971).",
  sprintf("  RcT infection optimum set to %d C (see blast_config.R section 4b).",
          EPIRICE_RCT_PEAK),
  "  Temperature response evidence: Kato, H. and Kozaka, T. (1974). Effect of",
  "  temperature on lesion enlargement and sporulation of Pyricularia oryzae in",
  "  rice leaves. Phytopathology 64: 828-830. doi:10.1094/Phyto-64-828.",
  "  Kato, H. (1974). Epidemiology of rice blast disease. Review of Plant",
  "  Protection Research 7: 1-20.",
  "BLASTAM (infection days): infection-warning model of Koshimizu, Y. (1988),",
  "  A forecasting method for occurrence of rice leaf blast with AMeDAS data,",
  "  Bull. Tohoku Natl. Agric. Exp. Stn. 78: 67-121; Hayashi & Koshimizu (1988)",
  "  ibid. 78: 123-138 [in Japanese]. A night is favourable when (1) leaf wetness",
  "  >= a temperature-dependent minimum (Barksdale & Jones 1965 curve, about",
  "  12 h at 16 C falling to about 8 h at 27 C; Koshimizu used a fixed 10 h),",
  sprintf("  (2) mean temp during wetness %d-%d C, and (3) the preceding 5-day mean",
          BLASTAM_TWET_MIN, BLASTAM_TWET_MAX),
  sprintf("  temp %d-%d C.", BLASTAM_PREV5_MIN, BLASTAM_PREV5_MAX),
  "  Deviations from Koshimizu: temperature-dependent wetness threshold (Barksdale",
  "  & Jones 1965); upper bounds raised from 25 C for tropical Australia; heavy",
  "  rain (>=4 mm/h) excluded from the wet count (Yoshino 1988); wetness from",
  "  hourly ERA5 RH (>=90%) rather than AMeDAS energy balance; night window in",
  "  local solar time.",
  "Prior Australian modelling: Lanoiselet, V., Cother, E.J. and Ash, G.J. (2002).",
  "  CLIMEX and DYMEX simulations of the potential occurrence of rice blast",
  "  disease in south-eastern Australia. Australasian Plant Pathology 31: 1-7.",
  "Weather: Open-Meteo ERA5 archive (data CC BY 4.0), non-commercial research use.",
  sep = "\n"
)

CITATION_MAP <- paste(
  "Model: EPIRICE (Savary et al. 2012, Crop Prot. 34:6-17). Code: epicrop (A.H. Sparks).",
  "BLASTAM: Koshimizu 1988. Weather: Open-Meteo ERA5 (CC BY 4.0).",
  sep = "\n"
)

# Standing caveat carried on every product. Lanoiselet et al. (2002) measured
# in-canopy RH at Yanco at least 20 percentage points above ambient, so both
# models read as conservative lower bounds in irrigated paddocks.
CAVEAT_CANOPY <- paste0(
  "Both models are driven by ERA5 ambient humidity, roughly a Stevenson screen ",
  "at field-bank height. In-canopy RH in irrigated rice was measured at Yanco ",
  "at least 20 percentage points above ambient (Lanoiselet et al. 2002), so ",
  "these outputs are better read as conservative lower bounds than as estimates."
)

################################################################################
# 8. Removed settings
################################################################################
#   SITES, EMERGENCE_DATE, GRID_RES, SMOOTH_FACTOR, EMAIL_*, SEND_EMAIL,
#   GRID_CONC*, GRID_FAIL_BACKOFF, GRID_PROBE_CHUNKS, GRID_BATCHES,
#   GRID_BATCH_CALLS, GRID_BATCH_WAIT_S, GRID_FETCH_TIMEOUT_S
#     as before.
#   TARGET_CALLS_PER_RUN   duplicated DAILY_WEIGHTED_CAP and was never read.
#   TOWN_FETCH_CORES       the fork pool is gone; towns are one batched request.

# Below this success fraction the retry pass is skipped: a wholesale failure is a
# bad day at the archive, not transient flakiness, and retrying doubles the waste.
GRID_RETRY_MIN_OK_FRAC <- 0.25

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
#
# CHANGES IN THIS REVISION
#   * ONE run date, from blast_run_date(). run_tag used to come from a second
#     Sys.Date() called after the fetch, so a run that started before local
#     midnight produced a map dated 30 July over weather to 23 July while the
#     town table, started after midnight, said 24 July.
#   * model_pt() derives emergence from model_end, not from the run's global
#     emergence. Under GRID_WINDOW_MODE = "coverage" the two differ, SEIR's
#     alignment check threw for every point, and the EPIRICE map came out empty
#     while the BLASTAM map still rendered.
#   * HEAT_STRETCH and BLASTAM_STRETCH are applied. They were documented as
#     working controls but no script read them, so both delivered maps rendered
#     as one flat pale blue against a fixed 2% and 21 day scale.
#   * The observed maximum is printed in the footer and passed to the email.
#   * The refresh phase is capped by the weighted budget, not only by the fetch
#     count, and the retry reserve is enforced on the fetch rather than only in
#     planning.
#   * Overlay line parts are split at long jumps. australia_roads.geojson holds a
#     feature with a 3.25 deg step from Victoria to Tasmania, which drew as a line
#     across Bass Strait on every map.
#   * Town labels are placed by a declutter pass and pushed inward near the east
#     coast, so "Gympie" is no longer clipped to "Gym".
#   * Optional COAST_MASK_KM blanks the partly marine coastal fringe, which was
#     carrying most of the BLASTAM signal in country where rice is not grown.
#   * The target lattice extent is rounded out to a whole number of cells, so the
#     northernmost row is no longer dropped by seq(-44, -10, by = 0.3).
#
# The test hook is GRID_ON_POINT(pid, lon, lat, hourly_dt) -> cache rows.
################################################################################

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

RUN_DATE <- blast_run_date()
run_tag  <- format(RUN_DATE, "%Y-%m-%d")
blastam_check_fetch_arithmetic()

key_of <- function(lon, lat) paste0(sprintf("%.4f", lon), "_", sprintf("%.4f", lat))
CACHE_COLS <- c("pid","lon","lat","date","TEMP","RHUM","RAIN","wet_hours","temp_wet","infect","semi")

# Round the requested extent OUT to a whole number of cells, so the lattice
# covers it. seq(-44, -10, by = 0.3) stops at -10.1, silently dropping the
# northernmost row of cells; the extra cells added here are ocean and are masked.
ext <- {
  fin <- GRID_RES_FINEST
  c(GRID_EXTENT[1],
    GRID_EXTENT[1] + ceiling((GRID_EXTENT[2] - GRID_EXTENT[1]) / fin) * fin,
    GRID_EXTENT[3],
    GRID_EXTENT[3] + ceiling((GRID_EXTENT[4] - GRID_EXTENT[3]) / fin) * fin)
}

# Hourly JSON for one point -> the daily cache rows both models read.
#
#  * lon is passed so the model day is cut at BLASTAM_DAY_CUT_HOUR local solar
#    and the BLASTAM night window sits inside it.
#  * keep_from drops the lead-in days, whose infect is NA because the lagged
#    preceding 5-day mean is not yet available. Writing them would overwrite good
#    cached values with NA on every short refresh.
#  * keep_to caps the rows at end_date. The fetch runs one day further, to
#    data_end, purely so the final model day is complete at every longitude.
make_on_point <- function(keep_from, keep_to) function(pid, lon, lat, hw) {
  w <- tryCatch(blastam_daily_from_hourly(hw, lon = lon), error = function(e) NULL)
  if (is.null(w) || nrow(w) == 0) return(NULL)
  w <- w[date >= keep_from & date <= keep_to]
  if (nrow(w) == 0) return(NULL)
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

  # PROGRESSIVE WITHIN-LEVEL ORDERING. A bit-reversed Morton (Z-order) index
  # makes ANY prefix a spatially uniform sample of the level, so a level that
  # runs out of budget partway thins the map evenly instead of leaving the
  # tropics on the coarse lattice while the south is fine.
  g[, `:=`(ci = as.integer(round((lon - ext[1]) / fin)),
           ri = as.integer(round((lat - ext[3]) / fin)))]
  nbits <- max(1L, as.integer(ceiling(log2(max(g$ci, g$ri) + 1))))
  bitrev <- function(x, bits) {
    r <- integer(length(x)); x <- as.integer(x)
    for (b in seq_len(bits)) {
      r <- bitwOr(bitwShiftL(r, 1L), bitwAnd(x, 1L)); x <- bitwShiftR(x, 1L)
    }
    r
  }
  morton <- function(i, j, bits) {
    k <- numeric(length(i))
    for (b in 0:(bits - 1L))
      k <- k + bitwAnd(bitwShiftR(i, b), 1L) * 2^(2 * b + 1) +
               bitwAnd(bitwShiftR(j, b), 1L) * 2^(2 * b)
    k
  }
  g[, ord := morton(bitrev(ci, nbits), bitrev(ri, nbits), nbits)]
  setorder(g, lvl, ord)
  g[, c("ci", "ri", "ord") := NULL]

  g[, pid := key_of(lon, lat)]
  g[]
}
targets <- build_targets()
LEVEL_RES <- c(GRID_RES_LEVELS, GRID_RES_FINEST)[seq_len(max(targets$lvl))]

# ---- Window ----------------------------------------------------------------
# data_end is the last day FETCHED; end_date is the last day MODELLED. They
# differ by DAY_CUT_LAG_DAYS because the model day is cut at
# BLASTAM_DAY_CUT_HOUR local solar, so the final fetched day is only partly
# covered and by an amount that depends on longitude. Dropping it makes the
# window identical at every longitude by construction.
data_end  <- RUN_DATE - ARCHIVE_LAG_DAYS
end_date  <- data_end - DAY_CUT_LAG_DAYS
emergence <- end_date - CROP_AGE_DAYS
win_dates <- seq(emergence, end_date, by = "day")
cat(sprintf("Run %s. Fetch to %s, model %s to %s (%d days); target %d land points at %.2f deg\n",
            run_tag, data_end, emergence, end_date, length(win_dates),
            nrow(targets), GRID_RES_FINEST))

# ---- Load cache ------------------------------------------------------------
OUT <- file.path(SCRIPT_DIR, OUTPUT_DIR)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
gz_file  <- file.path(OUT, WEATHER_CACHE_GZ)
csv_file <- file.path(OUT, WEATHER_CACHE_CSV)

# fread() cannot read .gz unless R.utils is installed (it is not in the runner
# container), so gz is decompressed with a base R gzfile connection first.
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
ver_file <- file.path(OUT, if (exists("CACHE_VERSION_FILE")) CACHE_VERSION_FILE else "cache_version.txt")
want_ver <- if (exists("CACHE_SCHEMA_VERSION")) CACHE_SCHEMA_VERSION else 1L
have_ver <- suppressWarnings(as.integer(tryCatch(readLines(ver_file, warn = FALSE)[1],
                                                 error = function(e) NA)))
if (is.na(have_ver)) have_ver <- 1L

cache_in <- read_cache()
cache <- cache_in$data
read_fmt <- if (is.na(cache_in$fmt)) "none" else cache_in$fmt
if (!is.null(cache) && have_ver != want_ver) {
  cat(sprintf("Cache schema is version %d, this code writes version %d; discarding %d cached point(s) and rebuilding.\n",
              have_ver, want_ver, length(unique(cache$pid))))
  cat("  (version 3 moved the model day cut to BLASTAM_DAY_CUT_HOUR, lagged the\n")
  cat("   preceding 5-day mean and tightened night completeness, so every TEMP,\n")
  cat("   RHUM, RAIN, infect and semi value from version 2 is superseded.)\n")
  cache <- NULL
  read_fmt <- "none"
}
if (is.null(cache) || !all(c("infect","semi","wet_hours") %in% names(cache))) {
  if (!is.null(cache)) cat("Cache lacks BLASTAM columns; rebuilding from scratch.\n")
  cache <- data.table(pid = character(), lon = numeric(), lat = numeric(),
                      date = as.Date(character()), TEMP = numeric(), RHUM = numeric(),
                      RAIN = numeric(), wet_hours = numeric(), temp_wet = numeric(),
                      infect = integer(), semi = integer())
}
cache[, date := as.Date(date)]

# Drop cached points that are not on the current target grid.
if (nrow(cache) > 0) {
  orphans <- setdiff(unique(cache$pid), targets$pid)
  if (length(orphans) > 0) {
    cat(sprintf("Dropping %d cached point(s) not on the current %.2f deg grid.\n",
                length(orphans), GRID_RES_FINEST))
    cache <- cache[!pid %in% orphans]
  }
}

# ---- Failure ledger --------------------------------------------------------
ledger_file <- file.path(OUT, if (exists("FAIL_LEDGER_FILE")) FAIL_LEDGER_FILE else "fetch_failures.csv")
max_strikes <- if (exists("FAIL_LEDGER_MAX_STRIKES")) FAIL_LEDGER_MAX_STRIKES else 4L
fails <- if (file.exists(ledger_file))
  tryCatch(fread(ledger_file, colClasses = list(character = "pid")), error = function(e) NULL) else NULL
if (is.null(fails) || !all(c("pid","strikes","last_try","last_status") %in% names(fails)))
  fails <- data.table(pid = character(), strikes = integer(),
                      last_try = as.Date(character()), last_status = character())
fails[, last_try := as.Date(last_try)]
benched <- fails[strikes >= max_strikes & last_try > (RUN_DATE - 2^pmin(strikes, 8L)), pid]
if (length(benched) > 0)
  cat(sprintf("%d point(s) benched after repeated fetch failures; they will be retried later.\n",
              length(benched)))

# ---- Decide what to fetch --------------------------------------------------
# The API charges max(1, nvars/10) * max(1, ndays/14) per location, with a 14 day
# FLOOR. Every fetch carries BLASTAM_LEADIN_DAYS of extra history that is
# computed and then discarded (1 day for the local solar shift, 5 for the lagged
# preceding 5-day mean) plus DAY_CUT_LAG_DAYS at the leading edge.
lead <- BLASTAM_LEADIN_DAYS
# In coverage mode the modelled window can end up to GRID_WINDOW_MAX_LAG_DAYS
# earlier than end_date, so it also STARTS that much earlier and the cache has to
# hold the extra days.
win_back <- if (identical(GRID_WINDOW_MODE, "coverage") &&
                exists("GRID_WINDOW_MAX_LAG_DAYS")) as.integer(GRID_WINDOW_MAX_LAG_DAYS) else 0L
add_keep_from   <- emergence - win_back
add_fetch_from  <- add_keep_from - lead
n_days_add <- as.integer(data_end - add_fetch_from) + 1L
cost_new   <- om_weight_per_location(n_days_add, length(OM_HOURLY_VARS))

tail_start      <- max(emergence, end_date - (REFRESH_TAIL_DAYS - 1L))
ref_keep_from   <- tail_start
ref_fetch_from  <- tail_start - lead
cost_ref   <- om_weight_per_location(as.integer(data_end - ref_fetch_from) + 1L,
                                     length(OM_HOURLY_VARS))

cached_pids <- unique(cache$pid)
last_by <- if (nrow(cache) > 0) cache[, .(last = max(date)), by = pid] else
  data.table(pid = character(), last = as.Date(character()))

max_fetch  <- GRID_MAX_FETCHES_PER_RUN
wt_cap     <- DAILY_WEIGHTED_CAP
stale_days <- REFRESH_MIN_STALE_DAYS

# Respect anything already spent today by another script or an earlier run.
spend_ledger <- file.path(OUT, SPEND_LEDGER_FILE)
already <- om_spend_read(spend_ledger)
wt_cap <- max(0, min(wt_cap, DAILY_WEIGHTED_HARD_CAP - already))
if (already > 0)
  cat(sprintf("Weighted ledger: %.0f already spent today, so this run is capped at %.0f.\n",
              already, wt_cap))

eligible <- last_by[last < end_date & last <= (end_date - stale_days)][order(last)]

# Hold back a slice of the weighted budget so the retry pass has something to
# spend. plan_cap is used for planning AND is now handed to the fetch, so a
# charged retry cannot silently eat the reserve: the 2026-07-30 run reported
# 8,667 weighted spent against a planned 8,550 for exactly that reason.
retry_wfrac <- if (exists("GRID_RETRY_WEIGHT_FRAC")) GRID_RETRY_WEIGHT_FRAC else 0.05
plan_cap <- wt_cap * (1 - retry_wfrac)

# A short tail can only close a gap it actually spans, so a point whose last row
# predates ref_keep_from - 1 is refetched over the FULL window instead.
tail_ok   <- eligible[last >= (ref_keep_from - 1L)]
tail_late <- eligible[last <  (ref_keep_from - 1L)]

# THE REFRESH PHASE IS NOW BUDGETED. It used to be min(nrow(tail_ok), max_fetch)
# with no weighted check, unlike the add and refetch phases. At cost_ref = 1.00
# and max_fetch = 8500 that fitted inside plan_cap by 50 calls; raise the tail so
# cost_ref becomes 1.07 and the planner would schedule a refresh it cannot pay
# for, the fetch would stop on "budget" partway through, and nothing would be
# added for the rest of the run.
n_refresh <- max(0L, min(nrow(tail_ok), max_fetch, floor(plan_cap / cost_ref)))
if (n_refresh < nrow(tail_ok))
  cat(sprintf("Refresh limited to %d of %d eligible points by the %s.\n",
              n_refresh, nrow(tail_ok),
              if (max_fetch <= floor(plan_cap / cost_ref)) "fetch count cap"
              else "weighted budget"))
maintain  <- targets[pid %in% tail_ok$pid[seq_len(n_refresh)]]
maintain_cost <- nrow(maintain) * cost_ref

n_restale <- max(0L, min(nrow(tail_late),
                         max_fetch - nrow(maintain),
                         floor(max(0, plan_cap - maintain_cost) / cost_new)))
restale <- if (n_restale > 0) targets[pid %in% tail_late$pid[seq_len(n_restale)]] else targets[0]
restale_cost <- nrow(restale) * cost_new

to_add <- targets[!pid %in% cached_pids & !pid %in% benched]
n_add <- max(0L, min(nrow(to_add),
                     max_fetch - nrow(maintain) - nrow(restale),
                     floor(max(0, plan_cap - maintain_cost - restale_cost) / cost_new)))
add <- if (n_add > 0) to_add[seq_len(n_add)] else to_add[0]

skipped_fresh <- length(cached_pids) - nrow(eligible)
cat(sprintf("Cache: %d points (%d fresh, skipped). Refresh %d @ %.2f, refetch %d stale @ %.2f, add %d new @ %.2f (%d/%d fetches, ~%.0f of %.0f weighted)\n",
            length(cached_pids), skipped_fresh, nrow(maintain), cost_ref,
            nrow(restale), cost_new, nrow(add), cost_new,
            nrow(maintain) + nrow(restale) + nrow(add), max_fetch,
            maintain_cost + restale_cost + nrow(add) * cost_new, wt_cap))
if (nrow(tail_late) > nrow(restale))
  cat(sprintf("  %d point(s) fell behind the %d day tail and are queued for a full refetch on a later run.\n",
              nrow(tail_late) - nrow(restale), REFRESH_TAIL_DAYS))

# ---- Deadlines -------------------------------------------------------------
reserve_min   <- GRID_RESERVE_MINUTES
retry_reserve <- if (exists("GRID_RETRY_RESERVE_MINUTES")) GRID_RETRY_RESERVE_MINUTES else 20L
if (!is.na(GRID_MAX_MINUTES) && GRID_MAX_MINUTES > 0)
  retry_reserve <- min(as.numeric(retry_reserve), as.numeric(GRID_MAX_MINUTES) * 0.25)

fetch_deadline <- if (!is.na(GRID_MAX_MINUTES) && GRID_MAX_MINUTES > 0)
  RUN_T0 + GRID_MAX_MINUTES * 60 else as.POSIXct(NA)
add_deadline <- if (is.na(fetch_deadline)) fetch_deadline else fetch_deadline - retry_reserve * 60
add_frac <- GRID_ADD_RESERVE_FRAC
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

run_phase <- function(tab, fetch_from, keep_from, deadline, label, cap) {
  if (nrow(tab) == 0L || isTRUE(quota_hit)) return(invisible())
  r <- fetch_points_batched(tab[, .(pid, lon, lat)], fetch_from, data_end,
                            on_point = make_on_point(keep_from, end_date),
                            deadline = deadline,
                            budget = max(0, cap - spent),
                            label = label)
  if (length(r$rows) > 0) new_rows <<- c(new_rows, r$rows)
  if (nrow(r$ledger) > 0) all_ledger[[length(all_ledger) + 1L]] <<- r$ledger
  spent <<- spent + r$spent
  if (identical(r$stopped, "quota")) quota_hit <<- TRUE
  invisible()
}

# Refresh first: under "latest" window mode a stale cell drops off the map.
run_phase(maintain, ref_fetch_from, ref_keep_from, refresh_deadline, "refresh", plan_cap)
run_phase(restale,  add_fetch_from, add_keep_from, refresh_deadline, "refetch", plan_cap)
run_phase(add,      add_fetch_from, add_keep_from, add_deadline,     "add",     plan_cap)

# ---- Retry -----------------------------------------------------------------
# Only points that were ATTEMPTED and failed with a transport or HTTP error.
# "empty" means the API answered and had nothing. The retry alone may draw on the
# reserved slice, which is why its cap is wt_cap rather than plan_cap.
led <- if (length(all_ledger)) rbindlist(all_ledger) else
  data.table(pid = character(), status = character(), code = integer())
retryable <- led[status %in% c("http", "transport"), unique(pid)]
attempted <- nrow(led)
ok_frac <- if (attempted > 0) nrow(led[status == "ok"]) / attempted else 1
retry_min_ok <- if (exists("GRID_RETRY_MIN_OK_FRAC")) GRID_RETRY_MIN_OK_FRAC else 0.25
if (length(retryable) > 0 && !quota_hit && ok_frac < retry_min_ok) {
  cat(sprintf("Skipping the retry pass: only %.0f%% of %d attempted points succeeded, which looks like a systemic failure rather than transient flakiness.\n",
              100 * ok_frac, attempted))
} else if (length(retryable) > 0 && !quota_hit) {
  rt_ref <- maintain[pid %in% retryable]
  rt_add <- rbind(restale[pid %in% retryable], add[pid %in% retryable])
  cat(sprintf("Retrying %d refresh and %d add/refetch point(s) that failed with a transport or HTTP error.\n",
              nrow(rt_ref), nrow(rt_add)))
  run_phase(rt_ref, ref_fetch_from, ref_keep_from, fetch_deadline, "refresh-retry", wt_cap)
  run_phase(rt_add, add_fetch_from, add_keep_from, fetch_deadline, "add-retry", wt_cap)
}

om_spend_add(spend_ledger, spent, "grid")

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
    upd[, `:=`(strikes = strikes + 1L, last_try = RUN_DATE, last_status = status)]
    fails <- rbind(fails[!pid %in% upd$pid],
                   upd[, .(pid, strikes, last_try, last_status)], fill = TRUE)
  }
  fwrite(fails, ledger_file)
  cat(sprintf("Failure ledger: %d point(s) carrying at least one strike.\n", nrow(fails)))
}

if (length(new_rows) > 0)
  cache <- rbindlist(list(cache, rbindlist(new_rows)), use.names = TRUE)

# Retain history rather than pruning to the modelling window. Weather already
# fetched costs nothing to keep.
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
writeLines(run_tag, file.path(OUT, "run_date.txt"))

# COMMON WINDOW, REAL WEATHER ONLY. Cells differ in how current they are, so the
# map's window ends at a date essentially every cell has reached, and each cell's
# series is TRUNCATED to it.
pt_end <- if (nrow(cache) > 0) cache[, .(mx = max(date)), by = pid] else
  data.table(pid = character(), mx = as.Date(character()))
n_cache_pts <- nrow(pt_end)
wmode <- GRID_WINDOW_MODE

# The newest date that at least `cover` of the cached cells have reached, never
# later than cap_date.
pick_window_end <- function(pt_end, cap_date, cover) {
  if (nrow(pt_end) == 0L) return(cap_date)
  mx <- sort(pt_end$mx, decreasing = TRUE)
  need <- max(1L, ceiling(cover * length(mx)))
  min(mx[need], cap_date)
}

min_cov  <- if (exists("GRID_WINDOW_MIN_COVERAGE")) GRID_WINDOW_MIN_COVERAGE else 0.90
reach_now <- if (n_cache_pts > 0L) sum(pt_end$mx >= end_date) / n_cache_pts else 1
window_note <- ""

if (identical(wmode, "latest") && reach_now >= min_cov) {
  model_end <- end_date
} else {
  # FALLBACK. Under "latest" this fires only when the run could not refresh
  # enough of the grid, which in practice means the weighted quota was already
  # spent. The map is still ONE window across every cell, just an older one.
  cover <- if (identical(wmode, "latest")) min_cov else GRID_WINDOW_COVERAGE
  model_end <- pick_window_end(pt_end, end_date, cover)
  if (identical(wmode, "latest") && model_end < end_date) {
    behind <- as.integer(end_date - model_end)
    window_note <- sprintf("window fell back %d %s to %s (only %.0f%% of cached cells reached %s)",
                           behind, if (behind == 1L) "day" else "days",
                           format(model_end), 100 * reach_now, format(end_date))
    cat(sprintf("WINDOW FALLBACK: %s.\n", window_note))
    cat("  This is the degraded-but-useful path: the run could not refresh the grid,\n")
    cat("  usually because the daily weighted quota was already spent, so the map is\n")
    cat("  built from cache at the newest date the cells actually share.\n")
    warn_days <- if (exists("GRID_WINDOW_WARN_FALLBACK_DAYS")) GRID_WINDOW_WARN_FALLBACK_DAYS else 10L
    if (behind > warn_days)
      cat(sprintf("  WARNING: that is more than %d days behind the archive edge. Check the quota ledger and the failure ledger.\n",
                  warn_days))
  }
}
# EVERY MAPPED CELL SHARES THIS EMERGENCE DATE. It used to be the run's global
# `emergence`, which equals model_end - CROP_AGE_DAYS only under "latest". Under
# "coverage" the weather was truncated to an earlier model_end while SEIR was
# still handed the later emergence, so its alignment check threw for every point
# and the EPIRICE map rendered empty while BLASTAM rendered normally.
model_start <- model_end - CROP_AGE_DAYS

current_pids <- pt_end[mx >= model_end, pid]
held_out <- n_cache_pts - length(current_pids)
# NOTE the >= . With `>` this yields CROP_AGE_DAYS rows starting one day after
# model_start, and SEIR's alignment check then throws for every point.
model_cache <- cache[pid %in% current_pids & date >= model_start & date <= model_end]
cat(sprintf("Window [%s] ends %s: %d of %d cells reach it, %d absent (%d days behind the archive edge).\n",
            wmode, format(model_end), length(current_pids), n_cache_pts, held_out,
            as.integer(end_date - model_end)))

# Distinguish "no EPIRICE because the alignment is wrong" from "no EPIRICE because
# the cache does not reach back to model_start yet". The second is expected on the
# first few coverage-mode runs and resolves itself as CACHE_HISTORY_DAYS fills.
if (nrow(model_cache) > 0) {
  reach <- model_cache[, .(first = min(date)), by = pid]
  n_reach <- sum(reach$first <= model_start)
  if (n_reach == 0L)
    cat(sprintf(paste0("No cell holds weather back to %s, so EPIRICE will be NA everywhere. ",
                       "The earliest cached day is %s. In coverage mode this is expected until ",
                       "the cache has accumulated GRID_WINDOW_MAX_LAG_DAYS (%s) of extra ",
                       "history; it is not the SEIR alignment fault.\n"),
                format(model_start), format(min(reach$first)),
                if (exists("GRID_WINDOW_MAX_LAG_DAYS")) GRID_WINDOW_MAX_LAG_DAYS else 0L))
}

# SEIR indexes the weather vector by POSITION, so a missing day would shift every
# later day by one and be modelled silently. Drop any point whose window is not a
# continuous run of dates. (SEIR now also refuses a gappy series itself.)
shape <- if (nrow(model_cache) > 0L)
  model_cache[, .(n = .N, span = as.integer(max(date) - min(date)) + 1L), by = pid] else
  data.table(pid = character(), n = integer(), span = integer())
gappy <- shape[n != span, pid]
if (length(gappy) > 0) {
  cat(sprintf("Dropping %d point(s) with gaps in the modelling window; they will be refetched.\n",
              length(gappy)))
  model_cache <- model_cache[!pid %in% gappy]
}
short <- if (nrow(shape) > 0L) shape[n == span & n < (CROP_AGE_DAYS + 1L), .N] else 0L
if (short > 0)
  cat(sprintf("%d point(s) have a short but continuous window; EPIRICE will report NA for them.\n", short))

# EPIRICE is reported only for a cell holding the FULL common window, so every
# coloured cell on the map carries the same crop age and the same number of days.
model_pt <- function(dt) {
  epi <- NA_real_
  if (nrow(dt) >= (CROP_AGE_DAYS + 1L) && min(dt$date) <= model_start) {
    w <- data.table(YYYYMMDD = dt$date, DOY = as.integer(format(dt$date, "%j")),
                    TEMP = dt$TEMP, RHUM = dt$RHUM, RAIN = dt$RAIN,
                    LAT = dt$lat[1], LON = dt$lon[1])
    setorder(w, YYYYMMDD)
    lb <- tryCatch(predict_leaf_blast(w, emergence = model_start,
                                      duration = as.integer(min(120L, nrow(w)))),
                   error = function(e) NULL)
    if (!is.null(lb) && nrow(lb) > 0) epi <- lb$intensity[nrow(lb)]
  }
  win  <- BLASTAM_WINDOW_DAYS
  lagd <- if (exists("BLASTAM_END_LAG_DAYS")) BLASTAM_END_LAG_DAYS else 0L
  bend <- model_end - lagd
  bs <- blastam_score(dt$infect, dt$semi, dt$date, bend,
                      window = win, recent = BLASTAM_RECENT_DAYS)
  list(intensity = epi, events = bs$events, unjudged = bs$unjudged)
}

pm <- model_cache[, { m <- model_pt(.SD); .(lon = lon[1], lat = lat[1],
              intensity = m$intensity, events = m$events, unjudged = m$unjudged) }, by = pid]
pm_epi <- pm[!is.na(intensity)]

# ---- Coverage reporting ----------------------------------------------------
land_area_deg2 <- nrow(targets) * GRID_RES_FINEST^2
map_spacing <- if (nrow(pm) > 0) sqrt(land_area_deg2 / nrow(pm)) else NA_real_

lvl_target <- targets[, .(N_target = .N), by = lvl]
lvl_done   <- targets[pid %in% pm$pid, .(N_done = .N), by = lvl]
lvl_tab <- merge(lvl_target, lvl_done, by = "lvl", all.x = TRUE)[order(lvl)]
lvl_tab[is.na(N_done), N_done := 0L]
lvl_tab[, frac := N_done / N_target]
complete_lvl <- suppressWarnings(max(lvl_tab[frac >= 0.995, lvl]))
res_complete <- if (is.finite(complete_lvl)) LEVEL_RES[complete_lvl] else NA_real_
obs_max_epi <- if (nrow(pm_epi) > 0) max(pm_epi$intensity, na.rm = TRUE) * 100 else NA_real_
obs_max_bl  <- if (nrow(pm) > 0) max(pm$events, na.rm = TRUE) else NA_real_
cat(sprintf("Level coverage: %s\n",
            paste(sprintf("L%d(%.2f deg) %d/%d", lvl_tab$lvl, LEVEL_RES[lvl_tab$lvl],
                          lvl_tab$N_done, lvl_tab$N_target), collapse = "; ")))
cat(sprintf("Modelled %d points to %s; %s; mean land spacing ~%s deg\n",
            nrow(pm), format(model_end),
            if (is.na(res_complete)) "no lattice level complete yet"
            else sprintf("complete to %.2f deg", res_complete),
            if (is.na(map_spacing)) "NA" else sprintf("%.2f", map_spacing)))
cat(sprintf("Observed maxima: EPIRICE %s, BLASTAM %s. Colour ceilings: %.2f%% and %d days.\n",
            if (is.na(obs_max_epi)) "NA" else sprintf("%.4f%%", obs_max_epi),
            if (is.na(obs_max_bl)) "NA" else sprintf("%.0f days", obs_max_bl),
            HEAT_MAX, as.integer(BLASTAM_HEAT_MAX)))

# ---- Render ----------------------------------------------------------------
# Overlay loader. australia_roads.geojson contains a feature whose last vertex
# jumps 3.25 deg from Victoria (144.67, -38.38) straight to Tasmania
# (146.33, -41.17), which drew as a line across Bass Strait on every map. Rather
# than editing the bundled data, any line part is split at a jump longer than
# OVERLAY_MAX_SEGMENT_DEG, which also catches two smaller artefacts.
split_long_segments <- function(v, maxd) {
  g <- tryCatch(as.data.frame(terra::geom(v)), error = function(e) NULL)
  if (is.null(g) || !nrow(g)) return(v)
  g$grp <- paste(g$geom, g$part, sep = "_")
  parts <- lapply(split(g, g$grp), function(sub) {
    if (nrow(sub) < 2L) return(NULL)
    step <- c(0, sqrt(diff(sub$x)^2 + diff(sub$y)^2))
    sub$seg <- cumsum(step > maxd)
    sub
  })
  g2 <- do.call(rbind, parts[!vapply(parts, is.null, logical(1))])
  if (is.null(g2) || !nrow(g2)) return(v)
  g2$obj <- as.integer(factor(paste(g2$grp, g2$seg)))
  keep <- g2$obj %in% as.integer(names(which(table(g2$obj) >= 2L)))
  g2 <- g2[keep, , drop = FALSE]
  if (!nrow(g2)) return(v)
  m <- cbind(object = as.integer(factor(g2$obj)), part = 1L, x = g2$x, y = g2$y, hole = 0L)
  out <- tryCatch(terra::vect(m, type = "lines", crs = terra::crs(v)),
                  error = function(e) NULL)
  if (is.null(out)) v else out
}

load_bundled <- function(fname) {
  f <- file.path(SCRIPT_DIR, fname)
  v <- tryCatch(if (file.exists(f)) terra::vect(f) else NULL, error = function(e) NULL)
  if (is.null(v) || nrow(v) == 0) return(NULL)
  v <- tryCatch(terra::crop(v, terra::ext(ext)), error = function(e) v)
  if (is.null(v) || nrow(v) == 0) return(NULL)
  maxd <- if (exists("OVERLAY_MAX_SEGMENT_DEG")) OVERLAY_MAX_SEGMENT_DEG else Inf
  if (is.finite(maxd)) v <- split_long_segments(v, maxd)
  v
}
rivers <- if (isTRUE(SHOW_RIVERS)) load_bundled("australia_rivers.geojson") else NULL
roads  <- if (isTRUE(SHOW_ROADS))  load_bundled("australia_roads.geojson") else NULL

# Greedy label declutter, south to north so the order is deterministic. Labels
# that would sit within LABEL_MIN_SEP_DEG of one already placed are dropped, and
# labels near the eastern edge are placed to the left so they are not clipped.
declutter_labels <- function(lon, lat, minsep) {
  keep <- logical(length(lon))
  px <- numeric(0); py <- numeric(0)
  for (i in order(lat)) {
    if (length(px) == 0L ||
        all(sqrt((lon[i] - px)^2 + (lat[i] - py)^2) >= minsep)) {
      keep[i] <- TRUE; px <- c(px, lon[i]); py <- c(py, lat[i])
    }
  }
  keep
}

# Apply the colour stretch. The ANCHOR IS UNCHANGED: hmax is still the deepest
# colour every week. Only the spacing changes, and the legend is labelled with
# true values, so nothing is misrepresented.
stretch_raster <- function(r, hmax, s) {
  if (!is.finite(s) || abs(s - 1) < 1e-9) return(r)
  z <- terra::clamp(r, 0, hmax, values = TRUE) / hmax
  z^s * hmax
}
legend_ticks <- function(hmax, s, n = 5L) {
  true_at <- pretty(c(0, hmax), n = n)
  true_at <- true_at[true_at >= 0 & true_at <= hmax]
  disp_at <- hmax * (true_at / hmax)^s
  list(at = disp_at, labels = format(true_at, trim = TRUE))
}

render_map <- function(pts, valcol, title, colours, hmax, stretch, legend, base,
                       obs_max = NA_real_, obs_fmt = "%.3f") {
  if (nrow(pts) < 3) {
    cat(sprintf("Too few points to render %s (%d with a value). No PNG this run.\n",
                base, nrow(pts)))
    return(FALSE)
  }
  r0 <- terra::rast(xmin = ext[1], xmax = ext[2], ymin = ext[3], ymax = ext[4],
                    resolution = GRID_RES_FINEST / 2, crs = "EPSG:4326")
  v <- terra::vect(as.matrix(pts[, .(lon, lat)]), type = "points", crs = "EPSG:4326")
  v$val <- pts[[valcol]]

  idw_rad <- min(if (exists("IDW_RADIUS_MAX")) IDW_RADIUS_MAX else 1.0,
                 max(GRID_RES_FINEST, IDW_RADIUS_MULT * map_spacing))
  r <- tryCatch(terra::interpIDW(r0, v, field = "val", radius = idw_rad,
                                 power = IDW_POWER, maxPoints = IDW_MAX_POINTS),
                error = function(e) terra::rasterize(v, r0, field = "val", fun = "mean"))
  names(r) <- valcol
  if (!is.null(land_poly)) r <- terra::mask(r, land_poly)
  if (isTRUE(MASK_UNCOVERED)) {
    dmask <- tryCatch(terra::distance(r0, v), error = function(e) NULL)
    if (!is.null(dmask)) r[dmask > idw_rad * 111320] <- NA
  }
  # Optional coastal fringe mask. ERA5 cells on the coast are partly marine, so
  # their humidity is not representative of any paddock, yet they carried most of
  # the BLASTAM signal on the delivered maps. Off by default.
  cmk <- if (exists("COAST_MASK_KM")) COAST_MASK_KM else 0
  if (is.finite(cmk) && cmk > 0 && !is.null(land_poly)) {
    cl <- tryCatch(terra::as.lines(land_poly), error = function(e) NULL)
    dcoast <- if (is.null(cl)) NULL else
      tryCatch(terra::distance(r0, cl), error = function(e) NULL)
    if (!is.null(dcoast)) r[dcoast < cmk * 1000] <- NA
  }

  # The GeoTIFF carries TRUE values; only the PNG is stretched.
  if (isTRUE(WRITE_GEOTIFF))
    terra::writeRaster(r, file.path(OUT, sprintf("%s_%s.tif", base, run_tag)),
                       overwrite = TRUE)

  rmax <- if (!is.null(hmax)) hmax else {
    m <- terra::global(r, "max", na.rm = TRUE)[1,1]; if (!is.finite(m) || m <= 0) 1 else m }
  rs <- stretch_raster(r, rmax, stretch)
  tk <- legend_ticks(rmax, stretch)

  # Filenames stay on the RUN date so send_email.py keeps finding them; the title
  # and footer carry the DATA date.
  png_file <- file.path(OUT, sprintf("%s_%s.png", base, run_tag))
  png(png_file, width = 1000, height = 900, res = 120)
  op <- par(mar = c(4.4, 4, 3, 5))
  ramp <- grDevices::colorRampPalette(colours)(200)
  plotted <- tryCatch({
    terra::plot(rs, col = ramp, range = c(0, rmax), xlab = "Longitude", ylab = "Latitude",
                main = sprintf("%s  weather to %s", title, format(model_end)),
                plg = list(title = legend, at = tk$at, labels = tk$labels))
    TRUE
  }, error = function(e) FALSE)
  if (!plotted)
    terra::plot(rs, col = ramp, range = c(0, rmax), xlab = "Longitude", ylab = "Latitude",
                main = sprintf("%s  weather to %s", title, format(model_end)),
                plg = list(title = legend))

  foot <- sprintf("run %s | %d cells | complete to %s deg | driver %s | colour scale 0 to %s%s",
                  run_tag, nrow(pts),
                  if (is.na(res_complete)) "no level" else sprintf("%.2f", res_complete),
                  OPENMETEO_MODEL, format(rmax, trim = TRUE),
                  if (abs(stretch - 1) > 1e-9) sprintf(", stretch %.2f", stretch) else "")
  if (isTRUE(SHOW_OBSERVED_MAX) && is.finite(obs_max))
    foot <- paste0(foot, " | observed max ", sprintf(obs_fmt, obs_max))
  mtext(foot, side = 1, line = 3.2, cex = 0.62, col = NSW_GREY_04)

  if (isTRUE(SHOW_RIVERS) && !is.null(rivers)) try(terra::lines(rivers, col = COL_RIVER, lwd = 0.6), silent = TRUE)
  if (isTRUE(SHOW_ROADS)  && !is.null(roads))  try(terra::lines(roads,  col = COL_ROAD,  lwd = 0.5), silent = TRUE)
  # Draw the coastline as polygon borders rather than terra::lines(), which is
  # more tolerant of multipart geometry.
  if (isTRUE(SHOW_COAST) && !is.null(land_poly))
    try(terra::plot(land_poly, add = TRUE, col = NA, border = COL_COAST, lwd = 1), silent = TRUE)

  if (isTRUE(SHOW_TOWNS) && exists("MONITOR_TOWNS") && nrow(MONITOR_TOWNS) > 0) {
    tw <- as.data.frame(MONITOR_TOWNS)
    points(tw$lon, tw$lat, pch = 21, bg = "white", col = COL_TOWN, cex = 1.0, lwd = 1.3)
    minsep <- if (exists("LABEL_MIN_SEP_DEG")) LABEL_MIN_SEP_DEG else 0.9
    keep <- declutter_labels(tw$lon, tw$lat, minsep)
    lab <- tw[keep, , drop = FALSE]
    # Push labels left near the eastern edge so they are not clipped.
    pos <- ifelse(lab$lon > (ext[2] - 6), 2, 4)
    text(lab$lon, lab$lat, lab$name, pos = pos, offset = 0.35,
         cex = if (exists("LABEL_CEX")) LABEL_CEX else 0.5, col = COL_TOWN)
    if (sum(!keep) > 0)
      mtext(sprintf("%d town label(s) suppressed to avoid overprinting; all %d towns are plotted.",
                    sum(!keep), nrow(tw)),
            side = 1, line = 3.9, cex = 0.55, col = NSW_GREY_04)
  }
  par(op); dev.off()
  file.copy(png_file, file.path(OUT, sprintf("%s_latest.png", base)), overwrite = TRUE)
  cat("Heatmap:", png_file, "\n")
  TRUE
}

pm_epi[, intensity_pct := intensity * 100]
rendered_epi <- render_map(pm_epi, "intensity_pct", "EPIRICE potential risk (%)", HEAT_COLOURS,
           HEAT_MAX, HEAT_STRETCH, "intensity %", "epirice_heatmap",
           obs_max = obs_max_epi, obs_fmt = "%.4f%%")
rendered_bl <- render_map(pm, "events", sprintf("BLASTAM infection days (last %d)", BLASTAM_WINDOW_DAYS),
           BLASTAM_HEAT_COLOURS, BLASTAM_HEAT_MAX, BLASTAM_STRETCH, "days",
           "blastam_heatmap", obs_max = obs_max_bl, obs_fmt = "%.0f days")

# ---- Save cache ------------------------------------------------------------
csvdt <- copy(cache)
csvdt[, `:=`(TEMP = round(TEMP, 1), RHUM = round(RHUM, 0), RAIN = round(RAIN, 1),
             temp_wet = round(temp_wet, 1), wet_hours = round(wet_hours, 0),
             lon = round(lon, 4), lat = round(lat, 4))]

# fwrite() decides whether to compress from the FILE EXTENSION, so a temp file
# named ".tmp" is written as plain text and then renamed to ".gz". gzfile()
# reads it back happily, so the verify passes and an 8x larger file is committed
# under a .gz name. Hence both the ".tmp.gz" naming and this explicit check.
is_gzip <- function(f) {
  con <- file(f, "rb"); on.exit(close(con), add = TRUE)
  identical(as.integer(readBin(con, "raw", 2L)), c(31L, 139L))
}

write_cache <- function(dt) {
  tmp_gz <- paste0(gz_file, ".tmp.gz")     # extension must survive, see is_gzip()
  ok_gz <- tryCatch({
    fwrite(dt, tmp_gz, na = "NA")
    if (!is_gzip(tmp_gz)) {
      cat("gz cache was NOT compressed (data.table built without zlib?); using plain CSV.\n")
      FALSE
    } else nrow(read_gz_dt(tmp_gz)) == nrow(dt)
  }, error = function(e) { cat("gz verify error:", conditionMessage(e), "\n"); FALSE })
  if (ok_gz) {
    file.rename(tmp_gz, gz_file)
    if (WEATHER_CACHE_KEEP_CSV) {
      tmp_csv <- paste0(csv_file, ".tmp.csv"); fwrite(dt, tmp_csv, na = "NA")
      file.rename(tmp_csv, csv_file); fmt <- "gz+csv"
    } else {
      if (file.exists(csv_file)) file.remove(csv_file); fmt <- "gz"
    }
    return(list(fmt = fmt, kb = file.info(gz_file)$size / 1024))
  }
  unlink(tmp_gz)
  cat("gz cache write/verify failed; falling back to plain CSV.\n")
  tmp_csv <- paste0(csv_file, ".tmp.csv"); fwrite(dt, tmp_csv, na = "NA")
  file.rename(tmp_csv, csv_file)
  if (file.exists(gz_file)) file.remove(gz_file)
  list(fmt = "csv", kb = file.info(csv_file)$size / 1024)
}
wc <- write_cache(csvdt)
writeLines(as.character(want_ver), ver_file)
cat(sprintf("Cache saved: %d points, %d rows as %s (%.0f KB)\n",
            length(unique(cache$pid)), nrow(cache), wc$fmt, wc$kb))

# ---- Stats line for the email ----------------------------------------------
# Field 3 is the number MAPPED on the previous run, read back before this run
# overwrites the file. Fields 11 and 12 are the observed maxima, so the email can
# say whether a flat map is flat weather or a broken scale.
stats_file <- file.path(OUT, "map_stats.txt")
prev_mapped <- if (!file.exists(stats_file)) NA_integer_ else tryCatch({
  s <- strsplit(readLines(stats_file, warn = FALSE)[1], "\\|")[[1]]
  as.integer(s[1])
}, error = function(e) NA_integer_, warning = function(w) NA_integer_)

# Fields 13 and 14 tell run_blast.R whether the window fell back and which maps
# actually rendered, so a degraded run explains itself in the email instead of
# silently arriving without attachments.
rendered <- paste(c(if (isTRUE(rendered_epi)) "epirice", if (isTRUE(rendered_bl)) "blastam"),
                  collapse = "+")
writeLines(sprintf("%d|%.2f|%d|%.2f|%s|%.0f|%s|%s|%s|%.0f|%s|%s|%s|%s",
                   nrow(pm), map_spacing,
                   if (is.na(prev_mapped)) 0L else prev_mapped,
                   GRID_RES_FINEST, wc$fmt, wc$kb, read_fmt,
                   format(model_end),
                   if (is.na(res_complete)) "" else sprintf("%.2f", res_complete),
                   spent,
                   if (is.na(obs_max_epi)) "" else sprintf("%.4f", obs_max_epi),
                   if (is.na(obs_max_bl))  "" else sprintf("%.0f", obs_max_bl),
                   gsub("[|\n]", " ", window_note),
                   rendered),
           stats_file)

if (Sys.getenv("BLAST_MIDWEEK") == "1") {
  added <- length(unique(cache$pid)) - length(cached_pids)
  writeLines(sprintf("%s|%d|%d|%d|%s|%.0f",
                     format(RUN_DATE), length(unique(cache$pid)),
                     length(cached_pids), added, wc$fmt, wc$kb),
             file.path(OUT, "midweek_status.txt"))
  cat(sprintf("Midweek fetch-only run: +%d points, cache now %d.\n",
              added, length(unique(cache$pid))))
}
cat(sprintf("Done. Spent ~%.0f of %.0f weighted calls this run.\n", spent, wt_cap))

#!/usr/bin/env Rscript
################################################################################
# map_window_edge.R   (IFS trial only; the product has no such script)
#
# WHICH DAY THE WEEKLY MAP IS DRAWN AT. The trial's weekly run, on Sunday UTC,
# fetches the 31 towns to the newest day there is and does NOT fetch the grid:
# a full map refresh costs about 7,700 weighted calls and Sunday's free allowance
# is the product's (see section 6b of blast_config.R). The map is drawn from the
# cache as the last top-up left it, normally Friday's.
#
# Left to itself the grid runner would aim at the towns' window, find that no
# cell reaches it, step back, and report the step as a fault every week. This
# script names the window the cache actually supports: the newest model day that
# GRID_WINDOW_COVERAGE of the cached cells have reached. The workflow hands it to
# the grid step as BLAST_DATA_END, so the runner draws a complete map at that day
# and says nothing untrue. If Friday's top-up did not run the day is older, and
# the health check (HEALTH_MAX_MAP_BEHIND_DAYS) says so.
#
# Prints EDGE=YYYY-MM-DD, or EDGE= when there is no cache to read.
################################################################################
suppressPackageStartupMessages(library(data.table))

# The data_end that makes the grid runner draw at the newest model day `cover`
# of the cached cells have reached: that day plus the day-cut lag, and never
# later than the run's own window. pt_end is data.table(pid, mx).
map_window_data_end <- function(pt_end, cap_end_date, cover, cut_lag = 1L) {
  if (is.null(pt_end) || nrow(pt_end) == 0L) return(as.Date(NA))
  pick_window_end(pt_end, cap_end_date, cover) + as.integer(cut_lag)
}

# pid and date of every cached row; NULL when the file is missing or unreadable.
# fread() cannot read .gz without R.utils, so it is decompressed first.
read_cache_days <- function(f) {
  if (!file.exists(f)) return(NULL)
  tryCatch({
    if (grepl("\\.gz$", f)) {
      tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
      inp <- gzfile(f, "rb"); outp <- file(tmp, "wb")
      repeat {
        chunk <- readBin(inp, "raw", 1048576L)
        if (length(chunk) == 0L) break
        writeBin(chunk, outp)
      }
      close(inp); close(outp); f <- tmp
    }
    d <- fread(f, select = c("pid", "date"), colClasses = list(character = "pid"))
    if (nrow(d) == 0L) NULL else d
  }, error = function(e) NULL)
}

if (sys.nframe() == 0L) {
  SCRIPT_DIR <- tryCatch(
    normalizePath(dirname(sys.frame(1)$ofile), winslash = "/"),
    error = function(e) normalizePath(getwd(), winslash = "/"))
  source(file.path(SCRIPT_DIR, "blast_config.R"))
  # THE IFS TRIAL'S COPY (see README.md in this directory). It must be run from
  # ifs/, where blast_config.R is the trial's; refuse anywhere else.
  if (!isTRUE(get0("IFS_TRIAL")))
    stop("this is the IFS trial's copy: run it from the ifs/ directory", call. = FALSE)
  source(file.path(SCRIPT_DIR, "grid_window.R"))
  say <- function(...) cat(..., "\n", sep = "", file = stderr())

  out <- file.path(SCRIPT_DIR, OUTPUT_DIR)
  d <- read_cache_days(file.path(out, WEATHER_CACHE_GZ))
  if (is.null(d)) d <- read_cache_days(file.path(out, WEATHER_CACHE_CSV))
  if (is.null(d)) {
    say("No readable cache; the map step keeps its default window.")
    cat("EDGE=\n")
  } else {
    pt_end <- d[, .(mx = max(as.Date(date))), by = pid]
    cap <- blast_data_end() - DAY_CUT_LAG_DAYS
    de <- map_window_data_end(pt_end, cap, GRID_WINDOW_COVERAGE, DAY_CUT_LAG_DAYS)
    say(sprintf("Cache: %d cells; %.0f%% of them reach %s, so the map is drawn at that day (the towns' window ends %s).",
                nrow(pt_end), 100 * GRID_WINDOW_COVERAGE, format(de - DAY_CUT_LAG_DAYS), format(cap)))
    cat("EDGE=", format(de), "\n", sep = "")
  }
}

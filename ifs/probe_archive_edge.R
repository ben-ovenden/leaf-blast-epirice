#!/usr/bin/env Rscript
################################################################################
# probe_archive_edge.R
#
# Ask the archive which day it can actually serve, and print it as
#
#   EDGE=YYYY-MM-DD
#
# on stdout (everything else goes to stderr). The workflows export the answer as
# BLAST_DATA_END, and blast_data_end() in blast_config.R then uses it for both
# the grid and the town run.
#
# WHY. The last day fetched used to be pure arithmetic: the earlier of the run
# date and the UTC date, minus ARCHIVE_LAG_DAYS. Open-Meteo publishes one more
# ERA5 day at about 00:30 to 01:00 UTC, so for the first part of every UTC day
# that arithmetic names a day that is not there yet. A cell at longitude L needs
# (10 - L/15) hours of the last day to complete its final model day, and the
# completeness rule allows two to be missing, so every cell west of 120 E then
# lands a day short while the east looks fine. That was the "89% of cached cells
# reached end_date" on every scheduled run of August and September 2026. Rather
# than guess when the update lands or when GitHub will fire the job, ask.
#
# The weekly run WAITS for the new day, up to ARCHIVE_EDGE_WAIT_MAX_MIN, because
# that day is the point of running on Monday UTC. A top-up sets
# BLAST_EDGE_WAIT_MIN=0 and takes whatever is published.
#
# If the probe cannot get an answer at all, the arithmetic edge is printed, which
# is the previous behaviour.
################################################################################
SCRIPT_DIR <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile), winslash = "/"),
  error = function(e) normalizePath(getwd(), winslash = "/"))

suppressPackageStartupMessages(library(data.table))
source(file.path(SCRIPT_DIR, "blast_config.R"))
# THE IFS TRIAL'S COPY of this script (see README.md in this directory). It must
# be run from ifs/, where blast_config.R is the trial's. From anywhere else it
# would load another config and write into another directory; refuse instead.
if (!isTRUE(get0("IFS_TRIAL")))
  stop("this is the IFS trial's copy: run it from the ifs/ directory", call. = FALSE)
source(file.path(SCRIPT_DIR, "openmeteo_batch.R"))

say <- function(...) cat(..., "\n", sep = "", file = stderr())

edge <- blast_archive_edge()        # the arithmetic edge; never the pinned value
wait_max <- suppressWarnings(as.numeric(Sys.getenv("BLAST_EDGE_WAIT_MIN", "")))
if (!isTRUE(is.finite(wait_max))) wait_max <- ARCHIVE_EDGE_WAIT_MAX_MIN
poll <- ARCHIVE_EDGE_POLL_MIN
# Hours of the last fetched day that the westernmost cell needs for its final
# model day: the day cut in local solar time, less the longitude in hours.
need <- as.integer(ceiling(BLASTAM_DAY_CUT_HOUR - GRID_EXTENT[1] / 15))

# IFS TRIAL: WHAT HAS HAPPENED, NOT WHAT IS FORECAST. The IFS archive holds
# forecast hours for days ahead, so the hourly probe below finds "data" for
# today whatever the time is. What makes today usable is that a forecast run
# started today has landed (see section 6b of blast_config.R); the archive's
# metadata names the newest run it holds. One small static file, not an API call.
# Returns that run's start in seconds since 1970, or NA.
ifs_newest_run <- function() {
  tryCatch({
    r <- curl::curl_fetch_memory(IFS_META_URL, handle = curl::new_handle(timeout = 30, connecttimeout = 20))
    if (!identical(as.integer(r$status_code), 200L)) return(NA_real_)
    as.numeric(jsonlite::fromJSON(rawToChar(r$content))$last_run_initialisation_time)
  }, error = function(e) NA_real_)
}

t0 <- Sys.time(); n_probe <- 0L; found <- as.Date(NA)
repeat {
  n_probe <- n_probe + 1L
  p <- om_probe_archive_edge(edge, ARCHIVE_EDGE_PROBE_LON, ARCHIVE_EDGE_PROBE_LAT, need)
  waited <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (identical(p$status, "ok") && !is.na(p$date)) {
    run_day <- ifs_settled_day(ifs_newest_run())
    # No word on the newest run: yesterday, whose first hours are always settled.
    settled <- if (is.na(run_day)) edge - IFS_UNPROBED_LAG_DAYS else run_day
    found <- min(p$date, settled)
    say(sprintf("Archive probe %d at %s UTC: data to %s, newest forecast run in the archive started %s, so the newest usable day is %s (wanted %s).",
                n_probe, format(Sys.time(), "%H:%M", tz = "UTC"), format(p$date),
                if (is.na(run_day)) "unknown" else format(run_day), format(found), format(edge)))
    if (found >= edge) break
  } else {
    say(sprintf("Archive probe %d at %s UTC got no answer (%s %s).",
                n_probe, format(Sys.time(), "%H:%M", tz = "UTC"), p$status, p$msg))
    if (identical(p$status, "quota")) break     # waiting on a spent quota helps nobody
  }
  if (waited + poll > wait_max) break
  say(sprintf("  waiting %d min for today's first forecast run to reach the archive (%.0f of at most %.0f min used).",
              as.integer(poll), waited, wait_max))
  Sys.sleep(poll * 60)
}

out <- if (is.na(found)) edge - IFS_UNPROBED_LAG_DAYS else min(found, edge)
if (is.na(found)) {
  say(sprintf("No usable answer from the archive; falling back to %s, yesterday, whose first hours are always settled.", format(out)))
} else if (out < edge) {
  say(sprintf("No forecast run started today has reached the archive yet; fetching to %s, so that no hour modelled is a forecast.",
              format(out)))
}
# A probe is one weighted call (seven days, three variables). Keep the ledger honest.
tryCatch(om_spend_add(file.path(SCRIPT_DIR, OUTPUT_DIR, SPEND_LEDGER_FILE), n_probe, "probe"),
         error = function(e) say("Spend ledger not updated: ", conditionMessage(e)))
cat("EDGE=", format(out), "\n", sep = "")

################################################################################
# run_health.R
#
# Is this run good enough to be called a success? One verdict, computed once at
# the end of run_blast.R from the town table and the grid's map_stats.txt, and
# written to blast_outputs/run_status.txt as key=value lines. Three readers:
#
#   * send_email.py prefixes the subject with "[DEGRADED]" and carries the town
#     window (and the map window when it differs) in the subject;
#   * the workflow's final "Verdict" step exits non-zero on degraded=1, AFTER the
#     cache is committed, the email sent and the artifact uploaded, so the run
#     shows red and GitHub notifies the workflow's owner, without losing anything;
#   * the email body, which names the town shortfall in the banner.
#
# WHY. The 2026-09-07 run delivered 31 towns of "no data" and the 09-14 and 09-21
# runs delivered a map three weeks stale, and every one of them exited 0 and
# showed a green tick. Nothing in the pipeline distinguished "ran" from "worked".
#
# The R scripts themselves still exit 0 on a degraded run, deliberately: a
# non-zero exit from run_blast_grid.R or run_blast.R would stop the workflow
# before the commit and the email, and the run that most needs explaining is the
# one that must be reported (see the 2026-07-31 note in blast_config.R). The
# loud failure comes last.
################################################################################
suppressPackageStartupMessages(library(data.table))

.rh_cfg <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default

# Returns list(degraded, reasons, warnings). `reasons` fail the run; `warnings`
# are noted but do not.
#
#   towns_modelled, towns_total   towns with an EPIRICE value / towns configured
#   map_behind_days               town end_date minus the map's window end
#   map_cells, map_cells_prev     cells modelled this run and last run
#   map_grey                      cells drawn grey (not refreshed to the window)
#   fetch_reason                  the grid fetch's own account, from map_stats
health_verdict <- function(towns_modelled, towns_total,
                           map_behind_days = 0L, map_cells = NA_integer_,
                           map_cells_prev = NA_integer_, map_grey = 0L,
                           fetch_reason = "",
                           min_town_frac = .rh_cfg("HEALTH_MIN_TOWN_FRAC", 0.9),
                           max_behind    = .rh_cfg("HEALTH_MAX_MAP_BEHIND_DAYS", 1L),
                           max_grey_frac = .rh_cfg("HEALTH_MAX_GREY_FRAC", 0.02),
                           max_drop_frac = .rh_cfg("HEALTH_MAX_MAP_DROP_FRAC", 0.10)) {
  reasons <- character(0); warnings <- character(0)
  towns_modelled <- as.integer(towns_modelled); towns_total <- as.integer(towns_total)
  if (isTRUE(towns_total > 0L)) {
    need <- as.integer(ceiling(min_town_frac * towns_total))
    if (towns_modelled < need)
      reasons <- c(reasons, sprintf("only %d of %d towns modelled", towns_modelled, towns_total))
  }
  mb <- suppressWarnings(as.integer(map_behind_days))
  if (isTRUE(mb > as.integer(max_behind)))
    reasons <- c(reasons, sprintf("map window %d day%s behind the town table", mb, if (mb == 1L) "" else "s"))
  mg <- suppressWarnings(as.integer(map_grey)); mc <- suppressWarnings(as.integer(map_cells))
  if (isTRUE(mg > 0L)) {
    tot <- if (isTRUE(is.finite(mc))) mc + mg else NA_integer_
    frac <- if (isTRUE(tot > 0L)) mg / tot else NA_real_
    if (isTRUE(frac > max_grey_frac))
      reasons <- c(reasons, sprintf("%d map cells (%.0f%%) not refreshed, drawn grey", mg, 100 * frac))
    else
      warnings <- c(warnings, sprintf("%d map cell(s) drawn grey", mg))
  }
  mp <- suppressWarnings(as.integer(map_cells_prev))
  if (isTRUE(is.finite(mc)) && isTRUE(mp > 0L) && mc < (1 - max_drop_frac) * mp)
    reasons <- c(reasons, sprintf("map cells fell from %d to %d", mp, mc))
  if (nzchar(fetch_reason %||% ""))
    warnings <- c(warnings, paste0("grid fetch: ", fetch_reason))
  list(degraded = length(reasons) > 0L, reasons = reasons, warnings = warnings)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || is.na(a[1])) b else a

# key=value lines. Values are flattened to one line; "|" and newlines removed so
# the file stays trivially parseable from bash and Python.
run_status_lines <- function(fields) {
  vapply(names(fields), function(k) {
    v <- fields[[k]]
    v <- if (is.null(v) || length(v) == 0L) "" else paste(as.character(v), collapse = "; ")
    sprintf("%s=%s", k, gsub("[\r\n|]", " ", v))
  }, character(1), USE.NAMES = FALSE)
}

read_run_status <- function(path) {
  if (!file.exists(path)) return(list())
  ln <- readLines(path, warn = FALSE)
  ln <- ln[grepl("=", ln, fixed = TRUE)]
  k <- sub("=.*$", "", ln); v <- sub("^[^=]*=", "", ln)
  setNames(as.list(v), k)
}

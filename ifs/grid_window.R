################################################################################
# grid_window.R
#
# Which date the continental map is drawn at, which cached cells are drawn GREY
# because they have not reached it, and why the run could not bring them
# current. Pure functions, so the policy is testable without a cache or a fetch.
#
# THE POLICY, for GRID_WINDOW_MODE = "latest"
#   1. If at least GRID_WINDOW_MIN_COVERAGE of cached cells reach end_date, draw
#      at end_date. The few that do not are drawn grey, not interpolated over.
#   2. Otherwise find the newest date that fraction does reach. If it is within
#      GRID_WINDOW_MAX_FALLBACK_DAYS, step back to it: a complete map a day or
#      two old beats a current one with holes.
#   3. If it is further back than that and at least GRID_WINDOW_MIN_DRAW_COVERAGE
#      of cells are current, draw at end_date and grey the stale cells, saying
#      how many there are and how old they are.
#   4. Below that the grid is mostly stale; step back to the older date, as
#      before, with the loud warning. An old complete map beats a mostly grey one.
#
# WHY. In September 2026 the fallback alone (rule 2 with no limit) fired with 76%
# of the grid current and drew EVERY cell at a 16 day old window three weeks
# running, while the email blamed a spent quota that had not been spent. The
# percentile rule discarded 5,890 fresh cells to honour 1,831 stale ones. It had
# been written for the 2026-07-31 case, where NO cell reached end_date and the
# choice was an old map or no map; that case is rule 4 and still works.
################################################################################
suppressPackageStartupMessages(library(data.table))

.gw_cfg <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default

# The newest date that at least `cover` of the cached cells have reached, never
# later than cap_date.
pick_window_end <- function(pt_end, cap_date, cover) {
  cap_date <- as.Date(cap_date)
  if (nrow(pt_end) == 0L) return(cap_date)
  mx <- sort(as.Date(pt_end$mx), decreasing = TRUE)
  need <- max(1L, ceiling(cover * length(mx)))
  min(mx[need], cap_date)
}

# pt_end: data.table(pid, mx), each cached cell's newest cached day.
# Returns model_end, the fraction reaching end_date, how far back the window
# stepped, which rule decided it ("current", "short-fallback", "draw-stale",
# "mostly-stale", "coverage") and the fallback sentence for the email, if any.
grid_choose_window <- function(pt_end, end_date, mode = "latest",
                               min_cov      = .gw_cfg("GRID_WINDOW_MIN_COVERAGE", 0.90),
                               draw_cov     = .gw_cfg("GRID_WINDOW_MIN_DRAW_COVERAGE", 0.50),
                               max_fallback = .gw_cfg("GRID_WINDOW_MAX_FALLBACK_DAYS", 3L),
                               coverage     = .gw_cfg("GRID_WINDOW_COVERAGE", 0.98)) {
  end_date <- as.Date(end_date)
  pt_end <- as.data.table(pt_end)
  n <- nrow(pt_end)
  reach_now <- if (n > 0L) sum(as.Date(pt_end$mx) >= end_date) / n else 1
  res <- list(model_end = end_date, reach_now = reach_now, behind = 0L,
              rule = "current", fallback_note = "")
  if (n == 0L) return(res)
  if (!identical(mode, "latest")) {
    res$model_end <- pick_window_end(pt_end, end_date, coverage)
    res$behind <- as.integer(end_date - res$model_end)
    res$rule <- "coverage"
    return(res)
  }
  if (reach_now >= min_cov) return(res)
  fb <- pick_window_end(pt_end, end_date, min_cov)
  behind <- as.integer(end_date - fb)
  if (behind <= as.integer(max_fallback) || reach_now < draw_cov) {
    res$model_end <- fb
    res$behind <- behind
    res$rule <- if (behind <= as.integer(max_fallback)) "short-fallback" else "mostly-stale"
    res$fallback_note <- sprintf("window fell back %d %s to %s (only %.0f%% of cached cells reached %s)",
                                 behind, if (behind == 1L) "day" else "days", format(fb),
                                 100 * reach_now, format(end_date))
    return(res)
  }
  res$rule <- "draw-stale"
  res
}

# Cells not reaching model_end, grouped by their newest cached day, newest first:
# "286 last updated 06 Sep, 555 30 Aug, 990 29 Aug".
grid_stale_summary <- function(pt_end, model_end, max_groups = 4L) {
  st <- as.data.table(pt_end)[as.Date(mx) < as.Date(model_end)]
  if (nrow(st) == 0L) return("")
  g <- st[, .N, by = .(mx = as.Date(mx))][order(-as.numeric(mx))]
  shown <- head(g, max_groups)
  txt <- paste(sprintf("%d %s%s", shown$N,
                       ifelse(seq_len(nrow(shown)) == 1L, "last updated ", ""),
                       format(shown$mx, "%d %b")), collapse = ", ")
  if (nrow(g) > max_groups)
    txt <- sprintf("%s, and %d older", txt, sum(g$N[-seq_len(max_groups)]))
  txt
}

# Why the run could not bring every cell current, from what the fetch actually
# recorded. The email used to assert "usually because the daily weather-API
# quota was already spent" whatever had happened; in September 2026 the quota
# had not been spent, the backlog was simply unaffordable at the old prices.
grid_fetch_reason <- function(stops = character(0), n_left = 0L, left_cost = 0, plan_cap = NA,
                              already = 0, wt_cap = NA, spent = 0, waited_s = 0,
                              n_failed = 0L, max_minutes = NA, held_out = 0L) {
  bits <- character(0)
  if (isTRUE(already > 0))
    bits <- c(bits, sprintf("an earlier run today had already spent %.0f weighted calls, capping this one at %.0f",
                            already, wt_cap))
  if (isTRUE(n_left > 0))
    bits <- c(bits, sprintf("%d cached cell(s) did not fit the %.0f weighted budget (about %.0f more needed)",
                            n_left, plan_cap, left_cost))
  if ("quota" %in% stops)
    bits <- c(bits, sprintf("the weather API refused further requests (HTTP 429) after ~%.0f weighted%s", spent,
                            if (isTRUE(waited_s > 0)) sprintf(", %.0f min of it spent waiting", waited_s / 60) else ""))
  if ("deadline" %in% stops)
    bits <- c(bits, sprintf("the fetch reached its %s minute wall-clock limit", format(max_minutes)))
  if ("budget" %in% stops)
    bits <- c(bits, "a fetch phase ran out of weighted budget part way")
  if (isTRUE(n_failed > 0))
    bits <- c(bits, sprintf("%d point(s) returned an HTTP or transport error, or no data", n_failed))
  if (length(bits) == 0L && isTRUE(held_out > 0))
    bits <- "every stale cell was requested this run, so those still behind returned incomplete or empty data"
  paste(bits, collapse = "; ")
}

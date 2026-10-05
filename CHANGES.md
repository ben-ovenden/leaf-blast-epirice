# Changes

## September 2026: three weeks of the same map

The 2026-09-07, 09-14 and 09-21 emails all carried the map "weather to
2026-08-29" (identical maxima, 0.058% and 10 days) beside a town table that was
current; the 09-07 email had 31 towns of "no data". Diagnosed from
`map_stats.txt`, `weighted_spend.csv`, `run_log.csv` and the committed cache. The
offline suite went from 62 tests to 212; all pass.

| # | Item | Fix | Guarded by |
| --- | --- | --- | --- |
| A | The Monday job fires at ~22:30 UTC on Sunday, so `run date − 6` was only five days behind the archive's clock and that day was not there yet. A cell at longitude L needs `10 − L/15` hours of it to complete `end_date`; two may be missing, so **every cell west of 120 E (841) landed one day short on every scheduled run**: "89% reached end_date" and a Degraded run banner every Monday, and a full-price refetch for each of them a week later. | `blast_data_end()` counts `ARCHIVE_LAG_DAYS` from the **earlier** of the run date and the UTC date. The workflow pins `BLAST_UTC_DATE` beside `BLAST_RUN_DATE` for the same reason the run date is pinned (the grid run crosses 00:00 UTC). Costs the Monday email one day of freshness. | test 16, including the 115 E versus 150 E geometry |
| B | Any cell more than `REFRESH_TAIL_DAYS` behind was refetched over the **full crop window at 4.86**. After 09-07 left 3,243 cells in that class, the weekly budget recovered ~800 per run and the 90% coverage rule pinned the map to the stale cohort's date: 65%, 74%, 76%. | `om_plan_refresh()`: cohorts by days behind, each fetched from the day after its newest row (8 days 1.07, 15 days 1.57, 22 days 2.07, never more than a new cell), cheapest first, one token bucket shared across the calls. The same backlog would have cost ~8,200 and fitted in one run. | test 17 |
| C | The 09-07 grid fetch hit the **hourly** 5,000 ceiling at exactly 5,000.0 weighted about an hour in; every 429 was treated as "quota spent, stop", so 42% of the grid was abandoned with 140 minutes of deadline unused, and the town fetch in the same hour was refused too. The 80/min pacer was 96% of the ceiling, retries were not paced, and the town run always lands in the grid's final hour. | Minutely and hourly 429s are waited out (`OM_QUOTA_WAIT_S`, bounded by `OM_QUOTA_WAIT_MAX_MIN` / `TOWN_QUOTA_WAIT_MAX_MIN` and the fetch deadline) and the same batch is resent; only a daily 429 stops a run. `GRID_TARGET_PER_MIN` 80 to 70. Retries go through the pacer. The town run's serial fallback, which is unpaced and unbudgeted, is skipped after a 429 that could not be waited out. Workflow timeout 240 to 255. | test 18 |
| D | With 76% of the grid current, the 90th-percentile fallback drew **every** cell at a 16 day old window, three weeks running, and IDW would otherwise have interpolated neighbours' values across the stale cells without a word. The banner ended "usually because the daily weather-API quota was already spent", which was untrue each time. | `grid_window.R`: the window steps back at most `GRID_WINDOW_MAX_FALLBACK_DAYS` (3) to bring every cell in; beyond that, if at least `GRID_WINDOW_MIN_DRAW_COVERAGE` (50%) are current, the map is drawn at `end_date` and the stale cells are **drawn grey**, masked out of the interpolation and the GeoTIFF, with count, share and ages in the footer, the legend and the email. Below 50% the old complete map still wins. `map_stats.txt` gains the grey count, the cohort summary and the fetch's own reason (budget, 429, deadline, failures, ledger), and the banner is built from those. | test 19, and test 13 now runs against the real `pick_window_end()` |
| F | Every one of those runs exited 0 and showed a green tick: 31 towns of "no data" on 09-07, a map three weeks stale on 09-14 and 09-21. The subject carried the **map** window alone, so it read "weather to 2026-08-29" above a town table modelled to a later date, and said nothing about the run being degraded. | `run_health.R`: one verdict per run from the town table and `map_stats.txt` (towns modelled below `HEALTH_MIN_TOWN_FRAC`, map more than `HEALTH_MAX_MAP_BEHIND_DAYS` behind the towns, more than `HEALTH_MAX_GREY_FRAC` grey, mapped cells down more than `HEALTH_MAX_MAP_DROP_FRAC`), written to `run_status.txt`. `send_email.py` puts the **town** window in the subject, adds the map window only when it differs, and prefixes `[DEGRADED]`. The banner names a town shortfall. The workflow gains a final `if: always()` **Verdict** step that exits non-zero on a missing, stale or degraded status, after the commit, email and upload, so the run turns red and GitHub notifies the owner without anything being lost. The R scripts still exit 0 on a degraded run, deliberately. | test 20 |
| E | The 2026-09-21 evening rerun, correctly capped by the spend ledger at 42 weighted for the towns, fell through to a **serial fallback** that was unpaced, uncharged and invisible to the ledger, fetched 29 towns off the books, and then **overwrote** the scheduled run's 31-town trends column and run log row with its worse result. The ledger itself had booked the Sunday run's whole 8,550 to Monday because the run crossed 00:00 UTC. | The second pass is a `fetch_points_batched()` call with `batch_size = 1`: same budget, pacer, ledger and 429 handling; `get_openmeteo_hourly()` is no longer called by any runner. `trends_merge_rerun()` fills a rerun's blanks from the earlier run's column and `runlog_keep_better()` keeps the row that modelled the most towns, so a rerun can improve the record and cannot worsen it. The ledger day is the pinned `BLAST_UTC_DATE`. | test 21 |
| G | The grid needs about twelve fetch days to fill from cold and several to recover from an interrupted run, and the only scheduled workflow was weekly. A `BLAST_MIDWEEK` branch existed but ran the whole pipeline, nothing invoked it, and the email line assumed a *daily* job. | `.github/workflows/midweek_topup.yml`: Thursday 00:30 UTC, same concurrency group as the weekly run, a different UTC quota day. `run_blast_grid.R` with `BLAST_MIDWEEK=1` fetches, merges, saves the cache, writes `midweek_status.txt` and quits before the window, models, maps or `map_stats.txt`. The Monday email reports the last top-up on a weekly cadence (`MIDWEEK_MAX_AGE_DAYS`). | test 22 |
| H | **Introduced by F, caught on its first live run.** The 2026-09-28 email was sent, correctly subjected `[DEGRADED] ... (weather to 2026-09-20)`, and the step then failed: `main()`'s `window_end` had been renamed `map_end` and the final "Email sent" print still used the old name, so `send_email.py` died with a `NameError` one line after `server.send_message()`. The run would have shown red every week, healthy or not. Nothing had ever executed the script before a live run; the suite only grepped it. | The print uses the town window. `BLAST_EMAIL_DRY_RUN=1` runs every line except the SMTP conversation and `BLAST_OUT_DIR` redirects the inputs; the offline suite runs the sender that way against fixture files, and the weekly workflow installs `python3` before the tests so the check is not skipped in CI. Verified with Python 3.12: the old file reproduces the failure with SMTP stubbed, the new one exits 0 through the same stub. | test 23 |

| I | The "7am" email of 28 September arrived at **11:16**. GitHub fired the 20:30 UTC cron 2 h 40 late (14 minutes to 2 h 40 late across August and September, later every week from 30 August) and the run now takes about 2 h 05. | Weekly crons moved from 06:30 to **01:47 local** (`47 15 * * 0` AEST, `47 14 * * 0` AEDT): on time lands about 04:00, three hours late about 07:00; off the hour and half hour; still Sunday in UTC so the window is unchanged. The gate job's literal schedule strings moved with them. | test 24: the crons and the gate name the same strings, the AEDT entry is one UTC hour earlier, and a run fired three hours late still lands by about 7am |
| J | **Repository growth.** Every run committed a new 6.5 MB cache to `main`, and git keeps every version of a committed file: after 38 runs the old caches were 104 MB of a 110 MB repository, and three runs a week would add about a gigabyte a year. Shortening the cache does not help; each version is still kept. | The cache and its schema marker live on a `cache-data` branch as **one parentless commit that every run replaces** (`git commit-tree`, force-push to that branch only). Each workflow restores it first and refuses to run without it, saves it straight after the fetch (so a later failure no longer costs the weather just fetched), and will not replace it with one under half its size on the same schema. `CACHE_HISTORY_DAYS` 120 to 90: the models need 61, 68 in coverage mode. | test 25 |
| K | **The archive's newest day is not there at midnight UTC.** The last day fetched was arithmetic, "UTC date minus 6", but Open-Meteo publishes that day at about 00:30 to 01:00 UTC (still missing at 00:08, 00:11, 00:16 and 00:20 UTC on 2 October; the 25 August run that began at 00:26 UTC caught it part way, 38 western cells a day short and 649 not). Early in a UTC day the run asks for a day that does not exist and every cell west of 120 E lands a day short. And the old Sunday-UTC schedule capped the map at eight days old. | `probe_archive_edge.R` asks one far-western point for the newest day whose first three hours are present and the workflows pin the answer as `BLAST_DATA_END`; `blast_data_end()` uses it when it is behind the arithmetic. The weekly run moved to **Monday UTC** (`47 0 * * 1`) and waits up to 90 minutes for the new day, so the map is modelled to the previous Monday, **seven days before the email instead of eight**. The email arrives on Monday afternoon instead of Monday morning. One cron in UTC; the two seasonal crons and the gate job are gone. | tests 24 and 26; test 16 for the geometry |
| L | One top-up a week, on Thursday, three days before the weekly run: if Monday's grid fetch was cut short the unreached cells were three days behind and the run was flagged. | Two top-ups, **Thursday and Sunday** at 03:17 UTC (after the archive's update). Thursday recovers from a bad Monday; Sunday, the latest day that does not share Monday's quota, leaves unreached cells one day behind, so the map steps back one day, complete and unflagged. `MIDWEEK_MAX_AGE_DAYS` 7 to 2, so the Monday email says when Sunday's did not run. The dispatch takes a `weighted_cap` for a three minute end to end test. | test 22 |
| M | The first weekly run on the new schedule was healthy and still carried the note "an earlier run today had already spent 1 weighted calls, capping this one at 9000": the archive probe's own call, read back from the ledger, against the configured cap. | The ledger is a reason only when it actually lowered the run's cap. | test 19 |
| N | **Town labels overprinted.** Each label sat on one fixed side of its marker and the only collision test was marker to marker distance, so "Humpty Doo" ran through Jabiru's marker and name, "Kununurra" through Timber Creek and "Goondiwindi" through Warwick, while five towns with room on another side had no label, and the footnote saying so was cut off by the bottom of the image. | `map_labels.R`: twelve candidate positions per marker, tested on what the text covers; most constrained label first, with one repair step, and checked at 85% to 115% of the measured text size because fonts differ between machines. The frame carries 2 degrees of sea on the east (`MAP_EAST_PAD_DEG`) so coastal labels have room. Each name has a one pixel white halo (`LABEL_HALO_IN`) so it stays readable across the coast, a road or a river. All 31 towns are labelled with nothing overprinted; the label note, if ever needed, is part of the single footer line. | test 27 |
| O | **The colour scale started by getting cooler**: pale grey for nothing, light blue for a little, then yellow, orange, red, so a patch of low risk read as less than the grey around it. | Near white, pale green, then the same yellow, orange and red. Chosen on renders of the 28 September data against two alternatives: white straight to yellow lost the low values, and blue as the zero colour turned them into pale holes. | test 28 |

**First live runs of the above (28 Sep and 1 Oct 2026).** Monday: map back on
the current window (20 Sep) with 7,272 cells modelled and the last 449 of the 29
August cohort drawn grey, as the planner's arithmetic says (9,479 needed against
8,550); 31 of 31 towns; no cell west of 120 E a day short. Verdict degraded, for
the 6% grey. Thursday top-up, first scheduled run: all 7,721 cells fetched for
8,330 weighted (7,272 at 1.00 plus 449 at 2.36), every cell at 24 September.

Verified end to end against a stubbed Open-Meteo in a scratch copy: a seeded
cache in the four cohorts of the live one (6, 14, 21 and 22 days behind) was
brought to `end_date` in a single run at exactly the planned 1.00 / 1.50 / 2.00 /
2.07 per cell with no calendar gaps, the leftover budget added new cells, an
injected hourly 429 was waited out and the run continued, and the town table
landed on the same window. A second run with a budget too small to clear the
backlog left 60 of 300 cells stale: the map was drawn at `end_date`, the 60 were
grey on both PNGs and NA in the GeoTIFF, and the email read "60 cells (20%) had
not reached 20 Sep and are drawn grey rather than interpolated over: 30 last
updated 30 Aug, 30 29 Aug. Why: 60 cached cell(s) did not fit the 266 weighted
budget (about 122 more needed)." That run's `run_status.txt` read `degraded=1`
with the grey-cell reason and the budget shortfall as a note; the full-recovery
run's read `degraded=0`. A third harness ran a Thursday top-up in a child process
(cache saved and grown, `midweek_status.txt` written, no maps, no `map_stats.txt`,
no `run_date.txt`, spend on the ledger) followed by the Monday run, whose email
read "Midweek top-up ran Thu 01 Oct: 375 cell(s) fetched, 75 new cell(s), cache
now 375". Before the ledger was keyed on the pinned date, that harness also
showed the town second pass fetching towns one at a time until the budget ran out
and the verdict calling the result degraded, which is what the live 09-21 rerun
should have said.

**Not done:** an "email only, no fetch" dispatch mode. The HTML body and the PNGs
live in the run artifact, not the repository, so resending a past email means
downloading that artifact; a mode that rebuilt them from the committed cache would
still need the towns' 150 weighted calls.

---

## Earlier: the code review and the email review

Every item from the code review and the email review, with where it was fixed and
how it is now guarded. The offline suite went from 20 tests to 54; all pass.

**Schema bump.** `CACHE_SCHEMA_VERSION` is now 3. The next grid run discards the
existing cache and rebuilds. Three stored quantities change: the model day cut,
the lagged preceding 5 day mean, and night completeness.

---

## Changed the numbers

| # | Item | Fix | Guarded by |
| --- | --- | --- | --- |
| 1 | EPIRICE collapse between the 2026-07-28 and 2026-07-29 trends columns. Schema 2 moved the daily aggregates onto local midnight days; `rainlim` is a daily **sum**, so nocturnal rain was split across two days and the 5 mm gate stopped opening. | `BLASTAM_DAY_CUT_HOUR` (10:00 local solar) in `blastam_model.R`. A day labelled 23 July runs 10:00 on 23 July to 09:59 on 24 July, so the night's wet period and rainfall stay in one day. | test 4 |
| 2 | Two weather dates in one email: maps said 23 Jul, body said 24 Jul. Three separate `Sys.Date()` calls, one of them after a two hour fetch. | `blast_run_date()` in `blast_config.R`, pinned once by the workflow as `BLAST_RUN_DATE` and read by both R scripts and `send_email.py`. | test 12 |
| 3 | RcT peaked at 20 C while the README argued for the published 25 C. Roughly a factor of two at northern Australian temperatures. | `EPIRICE_RCT_PEAK` (default 25) selects between two named curves in `epirice_model.R`. Recorded per run in `run_log.csv`. | test 6 |
| 4 | Town BLASTAM window was 22 days and its "7d" figure 8 days, disagreeing with the map. `blastam_score()` had no upper bound. | Window bounded at both ends inside `blastam_score()`; both products call it. | test 5 |
| 5 | The "preceding 5 day mean" included the night's own day. | Lagged one day, and computed on a complete date sequence so a gap yields NA rather than a mean spanning six or more days. Five leading NAs is now the signature. | test 2 |
| 6 | `HEAT_STRETCH` and `BLASTAM_STRETCH` were documented as working and read by nothing, so both maps rendered flat pale blue. | Applied in `render_map()` with a fixed anchor and true value legend ticks. Observed maximum printed in the footer and the email. | test 12 |
| 7 | `GRID_WINDOW_MODE = "coverage"` would render an empty EPIRICE map: SEIR got the run's global emergence while the weather was truncated to an earlier `model_end`. | `model_start <- model_end - CROP_AGE_DAYS`, plus `GRID_WINDOW_MAX_LAG_DAYS` of extra lookback so the earlier window start is actually cached. | test 7 |
| 8 | BLASTAM parameters could not be overridden from `blast_config.R`; the model file is sourced second and overwrote them unconditionally. | All parameters live in `blast_config.R`; the model file guards its fallbacks with `if (!exists(...))`. | |
| 9 | **New, found while testing:** `terra::shift` masks `data.table::shift`. The bare call failed, the caller's `tryCatch` turned it into every grid point coming back "empty", and the run produced a blank map with no error in the log. | `data.table::shift` fully qualified. | test 14, which attaches terra on purpose |
| 10 | **New, found while testing:** the run log failed to append on the second run, because dates written as character were read back as Date and `rbind` refused. | Class alignment before binding, wrapped so a metadata failure cannot kill a run whose outputs are on disk. | test 10 |

## Robustness

| # | Item | Fix |
| --- | --- | --- |
| 11 | A night was judged on one evening hour plus one morning hour, and an hour with missing humidity counted as dry. | `BLASTAM_MIN_EVE_HOURS`, `BLASTAM_MIN_MORN_HOURS`, `BLASTAM_MAX_NA_FRAC`; an hour is usable only when temperature, humidity and rain are all present. |
| 12 | A location with a null humidity column was accepted and mapped as dry weather. | `.om_hourly_dt()` rejects it, as it already did for temperature. |
| 13 | Only the grid runner screened for calendar gaps; SEIR indexes by position. | `run_blast.R` trims to the longest continuous run and reports what it dropped; `SEIR()` now refuses a gappy series itself. |
| 14 | The refresh phase was capped by fetch count only, unlike the add phases. | Capped by the weighted budget too, with an accurate message about which limit bound. |
| 15 | `GRID_RETRY_WEIGHT_FRAC` existed only in planning, so charged retries ate the reserve (8,667 spent against 8,550 planned). | `plan_cap` is handed to the fetch for the three main phases; only the retry pass may draw on `wt_cap`. |
| 16 | `DAILY_WEIGHTED_CAP` was per run, and the town fetch ran with an unlimited budget on top of it. | Shared `weighted_spend.csv` ledger keyed on the UTC quota day, combined ceiling `DAILY_WEIGHTED_HARD_CAP`. |
| 17 | `REFRESH_TAIL_DAYS` 8 + lead-in 6 = 14 no longer holds now there is a day cut lag. | Tail is 7; 7 + 6 + 1 = 14, so a refresh still costs exactly 1.00. `blastam_check_fetch_arithmetic()` warns if it drifts. |
| 18 | `seq(-44, -10, by = 0.3)` stopped at -10.1, dropping the northern row. | The lattice extent is rounded out to whole cells. At the current extent the extra cells are all ocean, so the land count is unchanged at 7,721; the lattice no longer under covers a hand edited extent. |
| 19 | `%||%` defined after its first use; `vector("list", ceiling(n / 1L))`; retry could overshoot the budget by one batch. | All three tidied in `openmeteo_batch.R`. |

## Maps

| # | Item | Fix |
| --- | --- | --- |
| 20 | A line ran across Bass Strait on every map: `australia_roads.geojson` has a feature whose last vertex jumps 3.25 deg from (144.67, -38.38) to (146.33, -41.17). | Overlay line parts are split at jumps longer than `OVERLAY_MAX_SEGMENT_DEG`, rather than editing the bundled data. Two smaller artefacts are caught too. |
| 21 | Town labels collided; Gympie was clipped to "Gym". | Greedy declutter at `LABEL_MIN_SEP_DEG`, labels pushed left near the eastern edge, and a footnote saying how many were suppressed. All towns are still plotted. |
| 22 | The BLASTAM signal sat on the partly marine coastal fringe. | `COAST_MASK_KM` blanks it at render time; cells are still fetched, cached and written to the GeoTIFF. **Off by default**, since it is a judgement call. |
| 23 | Footer illegible; no way to tell a flat map from a broken scale. | Larger, darker footer carrying the colour ceiling, the stretch and the observed maximum. |
| 24 | The coastline was drawn with `terra::lines()`. | Drawn as polygon borders, which is more tolerant of multipart geometry. |

## Email and CSVs

| # | Item | Fix |
| --- | --- | --- |
| 25 | The footnote claimed a fixed 10 h wetness threshold while the Barksdale and Jones curve was in use, and the same email also described the curve correctly. | Both prose blocks are generated from `BLASTAM_USE_BJ_THRESHOLD`. A test greps for a hard coded threshold. |
| 26 | Every town printed "0.00%" beside a "low" band; one decimal cannot resolve the 0.2% and 1% edges. | Three decimals throughout, and a note that intensity is `(diseased − removed) / (total sites − removed)`. |
| 27 | The 7 day column read "+0.00" for all 31 towns. | A dash for no change; header renamed "7 day change (pts)". |
| 28 | Config defined the NSW palette and the HTML used different hard coded colours. | Bands use NSW tints with a full strength brand dot, header on NSW Brand Blue, and the heat ramps end on NSW warning orange and error red. Contrast checked: no white text on orange. |
| 29 | `blast_unjudged` was computed and never surfaced, so "0 days" and "could not be judged" looked identical. | Asterisk per town, a count line, and a paragraph saying unjudged nights are not scored as unfavourable. |
| 30 | No BLASTAM legend to explain the purple column. | Six step legend matching the cell shading. |
| 31 | "era5" lowercase; hyphen in the sender name. | "ERA5"; sender is "WWAI Cereal Pathology: blast models". |
| 32 | The in canopy versus ambient caveat lived only in the README. | `CAVEAT_CANOPY` appears on the email and the text summary. |
| 33 | Trends columns keyed on the run date, so 28, 29 and 30 July each took a column describing near identical weather. | Keyed on the **data end date**; a re-run over the same window replaces its column. `HISTORY_RUNS` raised to 12. |
| 34 | Blank cells indistinguishable from zeros in a spreadsheet. | `na = "NA"` on every write. |
| 35 | Nothing recorded that the method changed between two columns. | `run_log.csv`: one row per run with the data window, schema version, RcT peak, day cut hour and BLASTAM bounds. Attached to the email. |
| 36 | Subject line gave no data window. | "Blast risk summary 2026-07-30 (weather to 2026-07-23)". |

## Citations

| # | Item | Fix |
| --- | --- | --- |
| 37 | The prose cited "Kato and Kozaka 1974" for the sporulation figures while the reference list held "Kato, H. (1974), *Review of Plant Protection Research* 7: 1 to 20". Two different papers; the list entry was the wrong one. | Confirmed by search: the sporulation source is **Kato, H. and Kozaka, T. (1974). Effect of temperature on lesion enlargement and sporulation of *Pyricularia oryzae* in rice leaves. *Phytopathology* 64: 828 to 830. doi:10.1094/Phyto-64-828.** Added to the README's EPIRICE references, the `CITATION` block carried on every email, `epirice_model.R` and `blastam_model.R`. The prose attribution was correct and is unchanged. |
| 38 | The Kato 1974 title was truncated to "Epidemiology of blast" in the README and given without a title in `blastam_model.R`. | Corrected to **Kato, H. (1974). Epidemiology of rice blast disease. *Review of Plant Protection Research* 7: 1 to 20.** in both places. |

Note that the second paper is in *Phytopathology*, not in a Japanese society
journal as its co-authorship might suggest. Both are now reachable from the
emailed citation block, so a reader does not have to go to the repository.

## Documentation

The README is rewritten so every claim matches the code: the intensity and
infection rate formulas, the 61 row / 60 day crop age arithmetic, the corrected
grid fill trajectory (about eleven **weekly** runs from cold, not four to nine
daily ones, since no midweek workflow exists in this repository), the three date
definitions, the weighted cost table, and a new note that the `rhlim` gate is a
24 hour mean of 90% and therefore almost never opens, so this configuration is
effectively rain driven.

`blast_config.R` section 8 no longer claims to have removed settings that were
still present: `TARGET_CALLS_PER_RUN` and `TOWN_FETCH_CORES` are gone, and the
`FREE_*` constants are now genuinely used by the spend ledger.

---

## Not done, and why

- **`SEIR()` internals.** `removed[d]` is read on the line above its assignment,
  `sum(infectious)` sums the whole pre-allocated vector, and `removed_today`
  reads `infday` carried from the previous iteration. These are reproduced from
  epicrop and are load bearing. They are now flagged in a comment. Diff against
  upstream before touching any of them.
- **`COAST_MASK_KM` left at 0.** Masking the coastal fringe is a scientific
  judgement, not a bug fix, and at 60 km it removes a great deal of the map.
- **A midweek top up workflow.** The grid needs about eleven weekly runs to fill
  from cold. `run_blast_grid.R` already has the `BLAST_MIDWEEK` branch and
  `run_blast.R` already reports on `midweek_status.txt`, so this is a workflow
  file away, but it is a new feature rather than a fix.

## Verification

R 4.3.3 with data.table 1.14.10, terra 1.7.65, jsonlite and curl. The offline
suite passes 54 of 54. Both runners were executed end to end against a stubbed
Open-Meteo returning synthetic weather in the real JSON shape, over three
consecutive weekly run dates: cold fill, then two refresh runs at exactly 1.00
weighted per point. Coverage window mode and the coastal mask were exercised
separately. Both PNGs render, the GeoTIFFs carry true values, trends accumulate
one column per data window, and the run log accumulates one row per run.

The synthetic weather is deliberately favourable, so the absolute numbers in
those test runs mean nothing. What they demonstrate is that the paths execute.

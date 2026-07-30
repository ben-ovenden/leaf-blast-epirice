# Changes

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

# Leaf blast risk: EPIRICE + BLASTAM, Open-Meteo

> **Cache schema version 3.** `run_blast_grid.R` will discard any version 1 or 2
> cache and rebuild it. Version 3 changes three stored quantities: the model day
> is now cut at 10:00 local solar rather than local midnight, so every `TEMP`,
> `RHUM` and `RAIN` value changes; the preceding 5 day mean is lagged so it
> genuinely precedes, so `infect` and `semi` change; and night completeness now
> requires a minimum number of observed hours. Expect about twelve runs
> before the 0.3 degree grid is full again, so on the weekly schedule the map
> will be coarse for a couple of months unless a midweek top up job is added.

A self contained pipeline that runs two complementary leaf blast models each week
from GitHub Actions, using free Open-Meteo ERA5 weather. Two risk maps, a 31 town
table and an HTML summary are emailed every Monday morning.

---

## What it does

One hourly weather fetch per grid point feeds both models at no extra API cost.

**EPIRICE** (Savary *et al.* 2012) is a mechanistic SEIR epidemic simulation.
From a rolling crop emergence date it steps day by day, tracking healthy, latent,
infectious and removed leaf sites, and reports **intensity**: the proportion of
leaf tissue diseased. It answers: *how much disease has the season built up?* It
is cumulative and slow to respond.

**BLASTAM** (Koshimizu 1988) is a Japanese infection warning model. For each
night it asks whether leaf wetness, night temperature and antecedent temperature
together favoured a new *Magnaporthe oryzae* infection event, and counts
**favourable infection nights in the last 21 days**. It answers: *where is
infection pressure building right now?*

The two are complementary. BLASTAM flags when infection windows open, useful for
timing fungicide decisions; EPIRICE estimates the disease that may follow. Both
are weather driven potentials, not field measurements.

### Relationship to previous Australian work

The most directly comparable prior study is Lanoiselet, Cother and Ash (2002),
which used CLIMEX and a custom DYMEX population model to ask whether *M. grisea*
could establish in the NSW rice belt. That study ran over four BOM station
locations (Finley, Griffith, Hay, Yanco) for 1988 to 1999 and found conditions
favourable in 10 of 11 seasons at one or more locations, with Yanco most at risk
and Griffith least, driven primarily by relative humidity.

This system extends that work. Coverage expands from four points in the southern
rice belt to about 7,700 land cells at 0.3 degrees across the whole continent,
including the tropical north and the wild *Oryza* country most relevant to
incursion pathways from Southeast Asia. The weather input is hourly ERA5
reanalysis rather than daily BOM min/max data, giving a real diurnal cycle. And it
runs operationally every week rather than as a one time risk assessment.

Both systems share the same fundamental limitation described in Lanoiselet's
paper: weather is measured at ambient height, but in canopy relative humidity in
irrigated paddocks averages at least 20 percentage points higher than ambient
(measured with data loggers at Yanco during the 2000 to 01 season). Both systems
therefore likely under count infection conducive hours. Quantifying the offset
requires in canopy loggers deployed alongside ERA5 driven model runs.

---

## Files

| File | Role |
| --- | --- |
| `blast_config.R` | **Every tunable parameter.** Sites, thresholds, grid, API, output, colours. Sourced first, and the model files guard their own defaults with `if (!exists(...))`, so a value set here always wins. |
| `epirice_model.R` | Vendored SEIR engine and leaf blast parameters from epicrop, with a selectable RcT optimum |
| `blastam_model.R` | BLASTAM infection warning model and the shared daily aggregator that feeds both models |
| `openmeteo_wth.R` | Single point Open-Meteo adapter (used by the `run_blast.R` fallback) |
| `openmeteo_batch.R` | Batched Open-Meteo fetcher, weighted cost model, pacer and the shared spend ledger |
| `run_blast.R` | Town table runner: fetch, model, write CSV, HTML and text summary |
| `run_blast_grid.R` | Continental heatmap runner: fill the cache, model, render maps |
| `send_email.py` | Python stdlib email sender |
| `test_offline.R` | Offline regression tests: 54 tests, no network, runs in seconds, in CI |
| `australia_land.geojson` | Land polygon for masking ocean and clipping the map |
| `australia_rivers.geojson` | River overlay |
| `australia_roads.geojson` | Road overlay |
| `.github/workflows/weekly_blast.yml` | Monday workflow: pin the run date, test, fetch, model, commit, email |

---

## Setup

### Secrets (GitHub, Settings, Secrets, Actions)

| Secret | Value |
| --- | --- |
| `MAIL_USERNAME` | Gmail address used to send |
| `MAIL_PASSWORD` | Gmail App Password (not the account password) |
| `MAIL_TO` | Comma separated recipient addresses |

The email is sent over SMTP SSL on port 465. App Passwords require 2FA on the
sending account.

### Configuration

Edit `blast_config.R` and commit. The main things to change:

- `MONITOR_TOWNS`, the 31 sites in the town table and trends CSVs.
- `CROP_AGE_DAYS`, the rolling crop age assumed everywhere (default 60).
- `INTENSITY_LOW_MAX`, `INTENSITY_MODERATE_MAX`, the band edges. **1% is the
  bottom of the high band**, and the band is open ended above. Calibrate against
  field observations; current values are provisional.
- `EPIRICE_RCT_PEAK`, the infection optimum temperature, 25 or 20. See below.
- `BLASTAM_DAY_CUT_HOUR`, where the 24 hour model day starts in local solar time.
- `HEAT_STRETCH`, `BLASTAM_STRETCH`, the colour ramp stretch.
- `COAST_MASK_KM`, optional blanking of the partly marine coastal fringe.
- `GRID_EXTENT`, the mapped area.

---

## Running locally

```
Rscript test_offline.R       # run this first
Rscript run_blast_grid.R     # continental heatmap
Rscript run_blast.R          # town table
```

Set `BLAST_RUN_DATE=YYYY-MM-DD` to pin the run date; otherwise today is used.
Required packages: `data.table`, `jsonlite`, `curl`, `terra`. The workflow uses
the `rocker/geospatial` container.

All 54 offline tests must pass before a run is meaningful. Each test guards a bug
that was actually shipped.

---

## Dates: three of them, and why

| Name | Definition | Meaning |
| --- | --- | --- |
| run date | `blast_run_date()`, pinned by the workflow | names the output files |
| `data_end` | run date minus `ARCHIVE_LAG_DAYS` (6) | the last day **fetched** |
| `end_date` | `data_end` minus `DAY_CUT_LAG_DAYS` (1) | the last day **modelled** |

`end_date` sits a day behind `data_end` because the model day is cut at 10:00
local solar, so the final fetched day is only partly covered, and by an amount
that depends on longitude. Dropping it makes the modelled window identical at
every longitude by construction rather than correcting for it afterwards.

**All three come from one pinned date.** The workflow resolves it once and
exports `BLAST_RUN_DATE`; both R scripts and `send_email.py` read it. Previously
each script called `Sys.Date()` separately, and the grid script called it twice,
once before a two hour fetch and once after, so a run straddling local midnight
produced maps titled "weather to 2026-07-23" beside body text saying "weather to
24 Jul 2026".

---

## Outputs

### Town table (`run_blast.R`)

- `blast_outputs/blast_results_latest.csv`, one row per town: EPIRICE intensity,
  BLASTAM favourable and semi favourable days, unjudged nights and a note field.
- `blast_outputs/blast_summary_latest.txt` and `.html`, the email bodies.
- `blast_outputs/town_trends.csv`, wide table, EPIRICE intensity as a percentage.
- `blast_outputs/blastam_trends.csv`, same layout, BLASTAM favourable days.
- `blast_outputs/run_log.csv`, one row per run recording the data window, cache
  schema, RcT peak, day cut hour and BLASTAM bounds.

**Trends columns are keyed on the DATA end date, not the run date**, so a re-run
over the same window replaces its column instead of adding one. Three test runs
on 28, 29 and 30 July previously took three columns describing almost the same
weather, and with a short history that evicts genuinely older columns. Blanks are
written as `NA`, because an empty cell is indistinguishable from a zero in a
spreadsheet.

Emergence is `end_date − CROP_AGE_DAYS` and therefore **moves with each run**.
The trends CSVs are a rolling 60 day window through time, not a season total.

`run_log.csv` exists because a change of method should not read as a change in
the weather. Between the 2026-07-28 and 2026-07-29 columns, Malanda fell from
0.374% to 0.006% and every other tropical town fell to zero. That was the schema
2 aggregation change, not an epidemiological collapse, and nothing in the file
said so.

### Heatmap (`run_blast_grid.R`)

- `blast_outputs/epirice_heatmap_YYYY-MM-DD.png` and `_latest.png`
- `blast_outputs/blastam_heatmap_YYYY-MM-DD.png` and `_latest.png`
- matching GeoTIFFs, which carry **true values**; only the PNG is stretched

The date in the filename is the run date; the title and footer carry the data
date. The heatmap colours every land cell as if rice were grown there.

### Internal state (committed, not emailed)

`weather_cache.csv.gz` (the cache), `cache_version.txt` (its schema version),
`fetch_failures.csv` (the failure ledger), `weighted_spend.csv` (the shared quota
ledger), `map_stats.txt` (the grid summary the email reads back) and
`run_date.txt` (the pinned run date). All are committed so the next run picks up
where this one stopped. Dated copies of the town CSV and text summary are written
alongside the `_latest` ones and are not committed.

**Colour scale.** `HEAT_MAX` (2%) and `BLASTAM_HEAT_MAX` (21 days) are fixed
ceilings so colours are comparable week to week. `HEAT_STRETCH` (0.4) and
`BLASTAM_STRETCH` (0.6) expand the low end: the anchor is unchanged, the legend
is labelled with true values, only the spacing of the colours changes.

These stretches are now applied. They were previously documented here as working
controls while no script read them, which is why the delivered maps rendered as
one flat pale blue. The highest town value in the 2026-07-30 email was Lismore
at 0.022% against a 2% ceiling, so roughly the first 1% of the ramp carried the
whole signal. The map's own maximum was not reported anywhere at the time, which
is exactly the gap the footer now fills: the observed maximum is printed on the
map and in the email, so a genuinely flat map is distinguishable from a broken
scale.

**Coastal cells.** ERA5 cells on the coastal fringe are partly marine, so their
humidity is not representative of any paddock, yet they carried most of the
BLASTAM signal on the delivered maps and drew the eye to country where rice is
not grown. Set `COAST_MASK_KM` to blank them at render time; the cells are still
fetched, cached and written to the GeoTIFF. Off by default.

**Overlays.** `australia_roads.geojson` contains a feature whose last vertex
jumps 3.25 degrees from Victoria to Tasmania, which drew as a line across Bass
Strait on every map. Line parts are split at jumps longer than
`OVERLAY_MAX_SEGMENT_DEG` rather than editing the bundled data.

---

## Weather cache and API

### Open-Meteo weighted calls

Open-Meteo charges in **weighted calls**, not requests:

```
weight per location = max(1, n_variables / 10) x max(1, n_days / 14)
```

The 14 day floor means a 1 day top up and a 14 day fetch cost identically.

| Fetch | Days | Weighted |
| --- | --- | --- |
| refresh an existing point | `REFRESH_TAIL_DAYS` 7 + lead-in 6 + day cut lag 1 = 14 | 1.00 |
| add a new point | crop window 61 + lead-in 6 + day cut lag 1 = 68 | 4.86 |

**`REFRESH_TAIL_DAYS + BLASTAM_LEADIN_DAYS + DAY_CUT_LAG_DAYS` must come to 14.**
At 15 the weight is 1.07, and a 7% surcharge on every refresh costs about 540
weighted calls once the grid is full, which is most of the headroom for adding
new cells. `blastam_check_fetch_arithmetic()` warns if the sum drifts and a test
asserts it.

Free tier limits are 10,000 weighted calls per day, 5,000 per hour and 600 per
minute. The grid run plans to `DAILY_WEIGHTED_CAP` (9,000) and paces at
`GRID_TARGET_PER_MIN` (80 weighted per minute, so 4,800 per hour). At about 4.86
weighted per new point, 80 per minute is roughly 16 points per minute.

**Shared spend ledger.** `DAILY_WEIGHTED_CAP` is per run. `run_blast.R` used to
fetch its towns with an unlimited budget on top of whatever the grid run had
already spent, and charged retries could push the grid past its own plan: the
2026-07-30 run reported 8,667 weighted against a planned 8,550. Both scripts now
append to `blast_outputs/weighted_spend.csv`, keyed on the UTC day the quota
resets, and read it back before setting their own budget, with a combined ceiling
of `DAILY_WEIGHTED_HARD_CAP` (9,500).

### Weather cache

`blast_outputs/weather_cache.csv.gz` stores daily EPIRICE inputs and BLASTAM
night judgements for every cached grid point, and is committed so each run only
fetches the latest days.

- Cached points need only the `REFRESH_TAIL_DAYS` tail, costing exactly 1.00.
- New points need the full window, costing about 4.86.
- The spare budget after refreshing is spent adding new points.
- The cache is written to a `.tmp.gz` temp path, read back to verify the row count
  and the gzip magic bytes, then renamed atomically.
- `CACHE_KEEP_HISTORY` retains `CACHE_HISTORY_DAYS` (120) beyond the modelling
  window. Do not raise this much while the cache is committed to git.
- `CACHE_SCHEMA_VERSION` (3): bump whenever a change alters the **values** stored,
  not just the columns. On a mismatch the cache is discarded entirely.

### Grid fill and resolution

The 7,721 land cell target at 0.3 degrees takes about twelve runs from cold,
because each run must refresh everything already cached before it can add
anything. Working the recursion through, adds per run are
`(plan_cap − n_cached) / 4.86`:

| Run | Points | Notes |
| --- | --- | --- |
| 1 (cold) | ~1,760 | 1.2 deg complete, 0.6 deg partial |
| 2 | ~3,160 | 0.6 deg complete |
| 4 | ~5,150 | 0.3 deg partial |
| 8 | ~7,200 | |
| 11 | ~7,710 | |
| 12 | 7,721 | 0.3 deg complete |

On the weekly schedule that is about three months. Earlier versions of this table
claimed four to nine "daily runs", which was wrong twice over: the arithmetic was
optimistic and the only scheduled workflow is weekly. `run_blast_grid.R` has a
`BLAST_MIDWEEK` branch and `run_blast.R` reports on `midweek_status.txt`, but a
top up workflow is not in this repository. Adding one is the way to make the fill
rate match the table.

Points are added in a bit reversed Morton (Z order) sequence within each
resolution level, so any partial run is a spatially uniform sample of the
continent rather than a south to north front.

`GRID_WINDOW_MODE = "latest"` (default) ends every cell at the archive edge on
the same date. `"coverage"` pulls the window back so nearly all cells are
included; it now works, and new points are fetched with
`GRID_WINDOW_MAX_LAG_DAYS` of extra lookback so the earlier window start is
actually in the cache. Previously that mode handed SEIR the run's global
emergence date while truncating the weather to an earlier `model_end`, so the
alignment check threw for every point and the EPIRICE map rendered empty while
the BLASTAM map rendered normally.

The lattice extent is rounded out to a whole number of cells before the grid is
built. `seq(-44, -10, by = 0.3)` stops at -10.1, silently dropping the northern
row. At the current extent the extra cells are all ocean, so the land cell count
is unchanged, but the lattice no longer under covers a hand edited extent.

### Failure ledger

`blast_outputs/fetch_failures.csv` tracks points that failed repeatedly. After
`FAIL_LEDGER_MAX_STRIKES` (4) consecutive failures a point is benched with
exponential cooloff.

---

## EPIRICE: parameters and deviations

The SEIR engine and leaf blast parameterisation are vendored from `epicrop`
(Sparks and colleagues) into `epirice_model.R`, so no package install is needed
and the version is fixed.

### How it works

Each day the model moves leaf sites through healthy (H), latent (E), infectious
(I) and removed (R). New infections per day are

```
infection = now_infectious x Rc x cofr^a
```

where `Rc = RcOpt x RcA x RcT x RcW` and `cofr` is the free site correction
`1 − diseased / (sites + diseased)`. Daily output is

```
intensity = (diseased − removed) / (total_sites − removed)
```

reported as a percentage. Earlier versions of this file described intensity as
`I / Sx` and the infection rate as `Rc x H x (I / Sx)^a`; neither matched the
code.

### Published parameters (Savary et al. 2012, Table 2)

| Parameter | Value | Description |
| --- | --- | --- |
| `onset` | 15 days | days after emergence before the epidemic can start |
| `duration` | 120 days | full season, run here over `CROP_AGE_DAYS + 1` = 61 rows |
| `rhlim` | 90% | daily **mean** RH at or above which leaves count as wet |
| `rainlim` | 5 mm | daily rainfall **sum** at or above which leaves count as wet |
| `H0` | 600 | initial healthy sites |
| `I0` | 1 | initial infective site at onset |
| `RcOpt` | 1.14 | optimum basic infection rate corrected for removals |
| `p` | 5 days | latent period |
| `i` | 20 days | infectious period |
| `a` | 1 | aggregation exponent on the free sites fraction |
| `Sx` | 30,000 | maximum leaf sites per unit area |
| `RRS` | 0.01 | relative senescence rate |
| `RRG` | 0.1 | relative leaf growth rate |

The simulation runs 61 rows so the **final** day carries a crop age of
`CROP_AGE_DAYS` (60).

### The wetness gate is effectively rain driven

`rhlim` is a 24 hour **mean** of 90%. That needs an all day saturated air mass: a
series with 95% nights and 58% afternoons averages about 75% and never opens the
humidity branch. In practice this configuration is driven almost entirely by the
5 mm rain branch, which is why `BLASTAM_DAY_CUT_HOUR` matters so much (below).
This is worth stating in any write up, because "EPIRICE intensity" reads as a
humidity driven quantity and here it is not.

### RcT optimum: 25 C by default, 20 C available

`EPIRICE_RCT_PEAK` in `blast_config.R` selects the curve. The default is the
published 25 C.

| Temp (C) | 10 | 15 | 20 | 25 | 30 | 35 | 40 | 45 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| peak 25 (default) | 0 | 0.5 | 0.6 | **1.0** | 0.6 | 0.2 | 0.05 | 0 |
| peak 20 (epicrop) | 0 | 0.5 | **1.0** | 0.6 | 0.2 | 0.05 | 0.01 | 0 |

The `epicrop` source carries a comment asserting that Table 2 of Savary *et al.*
2012 contains a typo and the optimum should be 20 C. This implementation uses the
published 25 C, because:

1. Empirical infection rate data support a 24 to 25 C optimum. Hashioka (1965),
   used in the Lanoiselet *et al.* (2002) DYMEX model for the Australian rice
   belt, measured the minimum time for conidial germination and penetration at
   6 h at 24 C, 8 h at 28 C and 10 h at 32 C. As a rate, infection is fastest at
   24 C and declines on both sides.
2. The broader literature reports 25 to 28 C for germination, infection, lesion
   formation and sporulation.
3. The consequence for northern Australia is large. At 28 C the 20 C curve gives
   `RcT` = 0.36 and the 25 C curve 0.76.
4. The epicrop note is a package author's reading of the paper, not a published
   erratum.

**The choice is not a uniform scaling.** Holding RH at 92% so the wetness gate is
open every day, final intensity over 61 days is:

| Daily mean temp (C) | 20 | 22 | 24 | 26 | 28 | 30 |
| --- | --- | --- | --- | --- | --- | --- |
| peak 20 | 4.73% | 2.65% | 1.34% | 0.59% | 0.21% | 0.05% |
| peak 25 | 0.91% | 1.91% | 3.58% | 3.58% | 1.91% | 0.91% |

Below about 22 C the 20 C curve gives the higher answer. Note also that
sporulation peaks cooler than infection (Kato and Kozaka 1974: 399 spores per
lesion per day at 20 C, 271 at 25 C, 131 at 32 C) and EPIRICE uses one `RcT` for
both, so either choice is a compromise. 25 favours the infection step, which is
the mechanistic basis for the SEIR transition.

**Until this revision the README argued for 25 C while the code shipped the 20 C
curve**, so the documented model and the running model differed by roughly a
factor of two at northern Australian temperatures. `run_log.csv` now records
which curve produced each trends column, and a test asserts the curve matches the
configured optimum.

### Other differences from a stock epicrop run

1. `duration` is `CROP_AGE_DAYS + 1` (61) rather than 120.
2. Weather is Open-Meteo ERA5 rather than NASA POWER.
3. `.calculate_audpc()` uses the standard trapezoidal definition. Affects only
   the AUDPC column, not the dynamics.
4. `SEIR()` now refuses a weather series with a calendar gap. It indexes the
   weather by position, so a missing day silently shifts every later day against
   crop age. Both runners screen for this; the engine now does too.

---

## BLASTAM: parameters and deviations

`blastam_model.R` implements the infection warning model of Koshimizu (1988),
operated across Japan on the AMeDAS network. It also produces the daily
aggregates EPIRICE reads, so one fetch serves both.

### How it works

For each night the model asks three questions:

1. Was the leaf wetness duration long enough?
2. Was the temperature during the wetness period within the favourable range?
3. Has the temperature of the **preceding** five days been warm enough?

If all three say yes the night is **favourable**. If the wetness bar is met but
exactly one temperature criterion fails it is **semi favourable**. Results are
stored as `infect` (0/1/NA) and `semi` (0/1/NA); NA means the night could not be
judged. The reported score counts favourable nights over the last
`BLASTAM_WINDOW_DAYS` (21) days.

### Original criteria (Koshimizu 1988)

| Criterion | Original value |
| --- | --- |
| Leaf wetness duration | at least 10 hours (fixed) |
| Temperature during wetness | 15 to 25 C |
| Preceding 5 day mean temperature | 20 to 25 C |
| Night window | 15:00 to 09:00, local Japan Standard Time |
| Leaf wetness estimation | energy balance from AMeDAS sunshine, wind, rainfall |
| Heavy rain | hours at or above 4 mm/h excluded (Yoshino 1988) |

### Deviations, and why

**1. Temperature dependent wetness threshold**, replacing the fixed 10 h.
Barksdale and Jones (1965), lower 95% CI: 12.2 h at 15.6 C, 10.9 at 18.3, 9.7 at
21.1, 8.6 at 23.9, 7.7 at 26.7. Linear interpolation with flat extrapolation, as
`blastam_bj_min_hours()`. The fixed 10 h was calibrated for Tohoku, where wetness
period temperatures are 15 to 20 C; in tropical Australia it over counts and in
the cool dry season it under counts. Restore with
`BLASTAM_USE_BJ_THRESHOLD <- FALSE`.

**2. Upper temperature bounds raised for tropical Australia.**

| Bound | Koshimizu 1988 | Here |
| --- | --- | --- |
| Wetness period mean, upper | 25 C | **32 C** |
| Preceding 5 day mean, upper | 25 C | **30 C** |
| Lower bounds | unchanged | unchanged |

The Japanese caps exclude the entire tropical wet season range where blast is
most active. Restore with `BLASTAM_TWET_MAX <- 25` and `BLASTAM_PREV5_MAX <- 25`
**in `blast_config.R`**, which now works: these parameters used to be set
unconditionally in `blastam_model.R`, which is sourced second, so anything set in
the config was silently overwritten.

**3. Heavy rain exclusion stated explicitly.** Hours at or above
`BLASTAM_RAIN_HEAVY` (4 mm/h) are excluded from the infection conducive wet
count, because heavy rain washes conidia off the leaf. Set to `Inf` to disable.

**4. Leaf wetness estimated from hourly RH.** An hour counts as wet at RH at or
above 90%, or rain at or above 0.2 mm/h, excluding heavy rain hours. ERA5 has no
leaf wetness variable. Absolute wet hour counts are provisional until calibrated
against field data.

**5. Local solar time rather than political timezone.** The 15:00 to 09:00 window
is applied in local solar time (longitude / 15 h), because political zone
boundaries would put a step discontinuity in wet hour counts across state lines.

### The model day is cut at 10:00 local solar

`BLASTAM_DAY_CUT_HOUR` (10) sets where the 24 hour model day starts. A day
labelled 23 July runs from 10:00 on 23 July to 09:59 on 24 July, local solar, so
the night window sits inside it by construction.

This matters far more than it looks. Schema 2 correctly moved the BLASTAM night
window onto local solar time but moved the **daily aggregates** with it, onto
midnight days. EPIRICE's `rainlim` gate is a daily **sum**, so cutting at local
midnight splits a nocturnal rain event across two days and halves the peak daily
total. On a synthetic series with 6 mm falling between 22:00 and 03:00 local
every third night, the number of days reaching the 5 mm gate went from 24 to 0
and final intensity from 0.0616% to 0.0000%, on identical rainfall. The measured
symptom was Malanda falling from 0.374% to 0.006% between the 2026-07-28 and
2026-07-29 runs, with every other tropical town going to zero at the same step
while the southern NSW sites did not move.

### The preceding 5 day mean genuinely precedes

`prev5` was `frollmean(TEMP, 5, align = "right")`, which covers days *i−4* to *i*
and therefore includes the night's own day. It is now lagged by one day, and
computed on a complete date sequence so a missing day yields NA rather than a
mean silently spanning six or more calendar days. Five leading NAs, not four, is
the signature of a correctly lagged window, and a test asserts it.

### Completeness and lead-in

A night is judged only when it has at least `BLASTAM_MIN_EVE_HOURS` (7) evening
and `BLASTAM_MIN_MORN_HOURS` (7) morning hours of usable data, no more than
`BLASTAM_MAX_NA_FRAC` (10%) unusable hours, and a preceding 5 day mean. An hour
is usable only when temperature, humidity and rainfall are all present.

Previously a night was judged on a single evening hour plus a single morning
hour, and an hour with missing humidity counted as **dry**, so a mostly empty
night scored "not favourable" rather than "not judged". A location whose humidity
column came back null is now rejected outright rather than mapped as dry weather.

Callers fetch `BLASTAM_LEADIN_DAYS` (6) extra days and discard them: one for the
local solar shift and five for the lagged preceding mean.

### Scores in the email

The email reports two numbers per town: favourable nights in the last 21 days,
and in parentheses the count over the last 7. Semi favourable nights and unjudged
nights are in `blast_results_latest.csv`, and unjudged nights are flagged in the
email with an asterisk.

Both the town table and the map now compute this through `blastam_score()`, whose
window is bounded at **both** ends. It used to test only `dates > (end_date −
window)`, so the town table reported a 22 day count and an 8 day "7 day" count
and disagreed with the map, which used the correct form.

---

## Known limitations

### Ambient versus in canopy humidity

Both models are driven by ERA5 ambient humidity, roughly a Stevenson screen at
field bank height. Lanoiselet *et al.* (2002) placed loggers in the rice canopy
at Yanco and found in canopy RH at least 20 percentage points above ambient. At
an ambient hourly RH of 72% the in canopy value would exceed the BLASTAM 90%
threshold. Both outputs should be read as conservative lower bounds. This caveat
now appears on the email itself, not only here.

### No validation against field outbreaks

Australia remains free of blast in cultivated rice, so direct field validation is
not possible here. Lanoiselet *et al.* (2002) validated their DYMEX model against
the California 1996 to 1999 outbreak with reasonable skill. A similar validation
of BLASTAM and EPIRICE against documented outbreaks in climatically similar
regions would strengthen confidence. The 25 C `RcT` choice is best tested against
tropical incidence data from the Philippines or northern Thailand.

### Rolling emergence

Emergence moves with every run, so neither output can be read as a season to date
progression. For season total analysis a fixed emergence date would be needed.

### Coastal ERA5 cells

Cells on the coastal fringe are partly marine. `COAST_MASK_KM` exists but is off
by default, so the delivered maps still include them.

---

## Workflow

The Monday workflow runs:

1. **Resolve run date**, pinned once and exported as `BLAST_RUN_DATE`.
2. **Offline tests**, `Rscript test_offline.R`. 54 tests, no network. terra is
   attached inside the suite on purpose, because `terra::shift` masks
   `data.table::shift` and that masking once turned every grid point into a
   silent "empty" and produced a blank map with no error in the log.
3. **Continental heatmaps**, `Rscript run_blast_grid.R`.
4. **Town table**, `Rscript run_blast.R`.
5. **Commit** the cache, trends, run log, spend ledger, failure ledger and map
   stats. The cache is committed before the email so a failed send does not lose
   the fetched data.
6. **Email**, `python3 send_email.py`, with the two heatmaps, the two trends CSVs
   and the run log attached.
7. **Upload artifact**, all of `blast_outputs/` kept 90 days as a fallback.

Two seasonal cron entries bracket the daylight saving change:

```
- cron: '30 20 * * 0'   # 06:30 Monday AEST (winter, UTC+10)
- cron: '30 19 * * 0'   # 06:30 Monday AEDT (summer, UTC+11)
```

The gate job matches the fired schedule against the current UTC offset, so a late
firing cron does not skip the week. A `concurrency` group prevents two
simultaneous runs from clobbering the cache.

---

## Attribution

Please retain the attribution comments in `epirice_model.R` and
`blastam_model.R`, and cite the primary publications if you publish or distribute
results.

**EPIRICE:** Savary, S., Nelson, A., Willocquet, L., Pangga, I., and Aunario, J.
(2012). Modeling and mapping potential epidemics of rice diseases globally. *Crop
Protection* 34: 6 to 17. doi:10.1016/j.cropro.2011.11.009

epicrop package (vendored SEIR engine and leaf blast parameters): Adam H. Sparks
and colleagues (DPIRD, WA; IRRI). Check the epicrop licence before
redistributing the model code. Framework: Zadoks, J.C. (1971). *Phytopathology*
61: 600 to 610.

Hashioka, Y. (1965). Effects of environmental factors on development of causal
fungus, infection, disease development, and epidemiology in rice blast disease.
In: *The Rice Blast Disease*. J Hopkins Press, pp. 153 to 161. [empirical basis
for the 24 to 25 C infection rate optimum]

Kato, H. and Kozaka, T. (1974). Effect of temperature on lesion enlargement and
sporulation of *Pyricularia oryzae* in rice leaves. *Phytopathology* 64: 828 to
830. doi:10.1094/Phyto-64-828 [source of the sporulation figures in the RcT
discussion above]

**BLASTAM:** Koshimizu, Y. (1988). *Bulletin of the Tohoku National Agricultural
Experiment Station* 78: 67 to 121 [in Japanese]. Hayashi, T. and Koshimizu, Y.
(1988). ibid. 78: 123 to 138 [in Japanese].

Barksdale, T.H. and Jones, M.W. (1965). Minimum conditions of temperature and dew
period for infection of rice by *Piricularia oryzae*. *Phytopathology* 55: 1037
to 1040.

Kato, H. (1974). Epidemiology of rice blast disease. *Review of Plant Protection
Research* 7: 1 to 20.

Maehara, H. and Yamada, M. (2025). Annual changes in the timing and frequency of
favorable conditions for rice leaf blast infection estimated by BLASTAM in
Fukushima Prefecture. *Annual Report of the Society of Plant Protection of North
Japan* 76: 41 to 46. doi:10.11455/kitanihon.2025.76_41

**Prior Australian modelling:** Lanoiselet, V., Cother, E.J. and Ash, G.J.
(2002). CLIMEX and DYMEX simulations of the potential occurrence of rice blast
disease in south-eastern Australia. *Australasian Plant Pathology* 31: 1 to 7.
doi:10.1071/AP01070

**Weather:** Open-Meteo historical ERA5 archive, data licensed CC BY 4.0. Non
commercial research use. Hersbach, H. *et al.* (2020). The ERA5 global
reanalysis. *Quarterly Journal of the Royal Meteorological Society* 146: 1999 to
2049. doi:10.1002/qj.3803

ERA5 native resolution is 0.25 degrees. A 0.3 degree model lattice matches this
driver; finer grids buy interpolation rather than information unless
`OPENMETEO_MODEL` is switched to `"era5_land"` (0.1 degrees). Confirm that
`relative_humidity_2m` is served by ERA5-Land before switching, and restate the
map resolution in any publications if you do.

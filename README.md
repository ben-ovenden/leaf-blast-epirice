# Leaf blast risk: EPIRICE + BLASTAM, Open-Meteo

> **Cache schema version 2.** If you are upgrading from an earlier version,
> `run_blast_grid.R` will discard the existing cache and rebuild it. The
> 2026-07-29 run confirmed that every infect/semi/wet\_hours/temp\_wet value
> written under schema version 1 was wrong (the BLASTAM night window was applied
> in UTC rather than local solar time) and cannot be reused. Expect four or five
> daily runs before the 0.3 degree grid is full again.

A self-contained pipeline that runs two complementary leaf blast models each
week from GitHub Actions, using free Open-Meteo ERA5 weather. Two risk maps, a
31-town table and an HTML summary are emailed every Monday morning.

---

## What it does

One hourly weather fetch per grid point feeds both models at no extra API cost.

**EPIRICE** (Savary *et al.* 2012) is a mechanistic SEIR epidemic simulation.
From a rolling crop emergence date it steps day by day, tracking healthy, latent,
infectious and removed leaf sites, and reports **intensity**: the proportion of
leaf tissue diseased (0 to 100%). It answers: *how much disease has the season
built up?* It is cumulative and slow to respond; it reads near zero in cool, dry
conditions regardless of recent infection events.

**BLASTAM** (Koshimizu 1988) is a Japanese infection-warning model. For each
night it asks whether leaf wetness, night temperature and antecedent temperature
together favoured a new *Magnaporthe oryzae* infection event, and counts
**favourable infection days in the last 21 days**. It answers: *where is
infection pressure building right now?* It responds quickly to humid, wet spells
and fires before disease is visible.

The two are complementary. BLASTAM flags when infection windows open, useful for
timing fungicide decisions; EPIRICE estimates the disease that may follow. Both
are weather-driven potentials, not field measurements.

### Relationship to previous Australian work

The most directly comparable prior study is Lanoiselet, Cother and Ash (2002),
which used CLIMEX (climate matching) and a custom DYMEX population model to ask
whether *M. grisea* could establish in the NSW rice belt. That study ran over
four BOM station locations (Finley, Griffith, Hay, Yanco) for the 1988–1999
period and found conditions favourable in 10 of 11 seasons at one or more
locations, with Yanco most at risk and Griffith least, driven primarily by
relative humidity.

This system extends that work in several ways. Coverage expands from four points
in the southern rice belt to approximately 7,700 land cells at 0.3 deg across
the whole continent, including the tropical north and the wild *Oryza* country
most relevant to incursion pathways from Southeast Asia. The weather input is
hourly ERA5 reanalysis rather than daily BOM min/max data, giving a real
diurnal cycle rather than a modelled sine wave. And it runs operationally every
week rather than as a one-time risk assessment.

The biological framing is different too. Lanoiselet's DYMEX model tracked spore
populations explicitly (production rate, UV mortality, latent period, per-lesion
fecundity) and asked how many infection events would occur in a season. This
system uses EPIRICE for the epidemic accumulation question and BLASTAM for the
infection-event question, making the two roles explicit rather than merged.

Both systems share the same fundamental limitation described in Lanoiselet's
paper: weather is measured at ambient height, but in-canopy relative humidity
in irrigated paddocks averages at least 20 percentage points higher than ambient
(measured directly with data loggers at Yanco during the 2000–01 season). Both
systems therefore likely under-count infection-conducive hours, particularly in
the Murrumbidgee and Murray irrigation areas where evapotranspiration from the
flooded surface adds substantially to the rice canopy microclimate. Quantifying
this offset requires in-canopy data loggers deployed alongside ERA5-driven model
runs.

---

## Files

| File | Role |
| --- | --- |
| `blast_config.R` | All configurable settings: sites, thresholds, grid, API, output |
| `epirice_model.R` | Vendored SEIR engine and leaf blast parameters from epicrop |
| `blastam_model.R` | BLASTAM infection-warning model, with full parameter notes |
| `openmeteo_wth.R` | Single-point Open-Meteo adapter (used by run\_blast.R fallback) |
| `openmeteo_batch.R` | Batched Open-Meteo fetcher used by both grid and town runs |
| `run_blast.R` | Town table runner: fetch, model, write CSV, HTML and text summary |
| `run_blast_grid.R` | Continental heatmap runner: fill the cache, model, render maps |
| `send_email.py` | Python stdlib email sender (replaces the broken Node 24 action) |
| `test_offline.R` | Offline regression tests (no network, runs in seconds, in CI) |
| `australia_land.geojson` | Land polygon for masking ocean and clipping the map |
| `australia_rivers.geojson` | River overlay |
| `australia_roads.geojson` | Road overlay |
| `.github/workflows/weekly_blast.yml` | Monday workflow: test → fetch → model → email → commit |

---

## Setup

### Secrets (GitHub → Settings → Secrets → Actions)

| Secret | Value |
| --- | --- |
| `MAIL_USERNAME` | Gmail address used to send (e.g. `wwai.pathology@gmail.com`) |
| `MAIL_PASSWORD` | Gmail App Password (not the account password) |
| `MAIL_TO` | Comma-separated recipient addresses |

The email is sent over SMTP SSL on port 465. App Passwords require 2FA enabled
on the sending account.

### Configuration

Edit `blast_config.R` and commit. The main things to change:

- `MONITOR_TOWNS` — the 31 sites tracked in the town table and trends CSVs.
- `CROP_AGE_DAYS` — the rolling crop age assumed everywhere (default 60 days).
- `INTENSITY_LOW_MAX`, `INTENSITY_MODERATE_MAX` — risk band thresholds. Calibrate
  against field observations; current values are provisional.
- `GRID_EXTENT` — the mapped area (default: Australian continent).

---

## Running locally

```r
source("blast_config.R")
source("epirice_model.R")
source("blastam_model.R")
source("openmeteo_wth.R")
source("openmeteo_batch.R")
Rscript run_blast.R       # town table
Rscript run_blast_grid.R  # continental heatmap
```

Required packages: `data.table`, `jsonlite`, `curl`, `terra`. The workflow
installs them inside the `rocker/geospatial` container; locally use
`install.packages(c("data.table","jsonlite","curl","terra"))`.

Run the offline tests first:

```r
Rscript test_offline.R
```

All 20 tests must pass before a run is meaningful. Each test guards a bug that
was actually shipped.

---

## Outputs

### Town table (`run_blast.R`)

- `blast_outputs/blast_results_latest.csv` — one row per town, EPIRICE intensity
  and BLASTAM infection days for the current run.
- `blast_outputs/blast_summary_latest.txt` — plain text summary for the email.
- `blast_outputs/blast_summary_latest.html` — HTML email body.
- `blast_outputs/town_trends.csv` — wide table, one column per run date,
  EPIRICE intensity (%). Rolling `HISTORY_RUNS` (default 10) runs kept.
- `blast_outputs/blastam_trends.csv` — same layout, BLASTAM favourable days.

**Important:** emergence is computed as `end_date − CROP_AGE_DAYS` and
therefore **moves with each run**. The trends CSVs are a rolling 60-day window
through time, not a cumulative season total.

### Heatmap (`run_blast_grid.R`)

- `blast_outputs/epirice_heatmap_YYYY-MM-DD.png` and `_latest.png` — EPIRICE
  potential risk surface.
- `blast_outputs/blastam_heatmap_YYYY-MM-DD.png` and `_latest.png` — BLASTAM
  infection-day surface.
- `blast_outputs/epirice_heatmap_YYYY-MM-DD.tif` — matching GeoTIFF (set
  `WRITE_GEOTIFF <- FALSE` to suppress).

The date in the filename is the **run date**, not the data date. The data window
ends 6 days before the run (ERA5 archive lag) and the map title states both.

The heatmap colours every land cell as if rice were grown there; it is a
potential risk surface driven by weather, not a map of actual crops or measured
disease. Interpret accordingly.

**Colour scale.** `HEAT_MAX` (default 2%) is the fixed ceiling for EPIRICE;
`BLASTAM_HEAT_MAX` (default 21 days) is the fixed ceiling for BLASTAM. Both
stay fixed week to week so colours are directly comparable across runs. A flat
blue map in July is the correct signal for the Australian winter; the scale is
not wrong. Set `HEAT_STRETCH <- 0.5` for a square-root stretch that makes low
values more legible without changing the anchor.

---

## Weather cache and API

### Open-Meteo weighted calls

Open-Meteo charges in **weighted calls**, not in requests:

```
weight per location = max(1, n_variables / 10) × max(1, n_days / 14)
```

The 14-day floor means a 1-day top-up and a 14-day fetch cost identically.
Three variables over 67 days costs ~4.8 weighted calls per point; a 14-day
refresh costs exactly 1.0.

Free tier limits: 10,000 weighted calls per day, 5,000 per hour, 600 per
minute. This pipeline uses at most 9,000 per day (`DAILY_WEIGHTED_CAP`) and
paces at 80 per minute (`GRID_TARGET_PER_MIN`) to stay under the hourly cap.
`GRID_TARGET_PER_MIN` is in **weighted calls per minute**, not fetches; at ~4.8
weighted per point, 80/min is about 17 points per minute.

Fetching is batched: `OM_BATCH_SIZE` (25) locations per request. Weighted cost
is unchanged but HTTP round trips fall by a factor of 25, which is what governs
wall clock. The `openmeteo_batch.R` fetcher uses `curl` directly so HTTP 429
quota-exceeded responses are distinguishable from transport failures and stop the
run immediately rather than burning three hours confirming the quota is spent.

### Weather cache

`blast_outputs/weather_cache.csv.gz` stores daily EPIRICE inputs and BLASTAM
night judgements for every cached grid point. The cache is committed to the
repository so each run only fetches the latest days, not the full history.

- Points already cached need only an 8-day tail fetched (`REFRESH_TAIL_DAYS`).
  Combined with the 6-day BLASTAM lead-in, the actual fetch is 14 days, costing
  exactly 1.0 weighted calls per point.
- New points need the full window: `CROP_AGE_DAYS` (60) + lead-in (6) = 67
  days, costing ~4.8 weighted per point.
- The spare budget after refreshing all existing points is spent adding new ones.
- The cache is written to a `.tmp.gz` temp path, read back to verify the row
  count and the gzip magic bytes, then renamed atomically, so a cancelled job
  cannot corrupt it.
- `CACHE_KEEP_HISTORY` retains `CACHE_HISTORY_DAYS` (120) of weather beyond the
  modelling window, at roughly 25 MB per year once the grid is full. Do not raise
  this significantly while the cache is committed to git.
- `CACHE_SCHEMA_VERSION` (currently 2): bump this whenever a change alters the
  **values** stored in the cache, not just its columns. On a mismatch the cache
  is discarded entirely rather than partially refreshed, which would mix old and
  new definitions on the same map.

### Failure ledger

`blast_outputs/fetch_failures.csv` tracks points that failed repeatedly. After
`FAIL_LEDGER_MAX_STRIKES` (4) consecutive failures a point is benched with
exponential cooloff, so it is not retried every run when the failure is
structural (offshore cell, bad coordinate).

### Grid fill and resolution

The 7,721 land-cell target at 0.3 deg takes about 8–9 daily runs from cold.
Subsequent runs refresh the whole grid in one pass. Progress:

| Run | Points | Spacing |
| --- | --- | --- |
| 1 (cold) | ~1,800 | 0.60 deg complete |
| 2 | ~3,600 | 0.60 deg complete |
| 4 | ~6,600 | 0.60 complete, 0.30 partial |
| 8–9 | 7,721 | 0.30 deg complete |

Points are added in a bit-reversed Morton (Z-order) sequence within each
resolution level, so any partial run is a spatially uniform sample of the
continent rather than a south-to-north front. Cape York, the Top End and the
Kimberley receive proportional coverage from the first run.

`GRID_WINDOW_MODE = "latest"` (default) ends every cell at the archive edge on
the same date, keeping the map comparable across longitude. Cells that missed
the refresh this run are absent rather than stale.

---

## EPIRICE: parameters and deviations

The SEIR engine and leaf blast parameterisation are vendored from the `epicrop`
package (Sparks and colleagues) into `epirice_model.R`. This means no package
install is needed and the version is fixed.

### How it works

Each day the model moves leaf sites through healthy (H), latent/exposed (E),
infectious (I) and removed (R) states. New infections per day = Rc × H ×
(I / Sx)^a, where Rc is the basic infection rate modulated by four multipliers:

- **RcA**: crop age modifier (young tissue is most susceptible)
- **RcT**: temperature modifier (peaks at the optimum)
- **Wetness switch**: on when daily mean RH ≥ `rhlim` or daily rainfall ≥
  `rainlim`; off otherwise (no wetness = no infection regardless of temperature)
- **Free-site fraction**: infection rate falls as tissue becomes occupied

Daily output is intensity = I / Sx (infectious sites / maximum sites), reported
as a percentage.

### Published parameters (Savary et al. 2012, Table 2)

| Parameter | Value | Description |
| --- | --- | --- |
| `onset` | 15 days | Days after emergence before the epidemic can start |
| `duration` | 120 days | Full season (run here at `CROP_AGE_DAYS` = 60) |
| `rhlim` | 90% | Daily mean RH at or above which leaves count as wet |
| `rainlim` | 5 mm | Daily rainfall at or above which leaves count as wet |
| `H0` | 600 | Initial healthy sites |
| `I0` | 1 | Initial infective site at onset |
| `RcOpt` | 1.14 | Optimum basic infection rate corrected for removals |
| `p` | 5 days | Latent period |
| `i` | 20 days | Infectious period |
| `a` | 1 | Aggregation exponent on the free-sites fraction |
| `Sx` | 30,000 | Maximum leaf sites per unit area |
| `RRS` | 0.01 | Relative senescence rate |
| `RRG` | 0.1 | Relative leaf growth rate |

**Temperature response `RcT`** — relative infection rate vs daily mean
air temperature (as published, Table 2 of Savary *et al.* 2012):

| Temp (°C) | 10 | 15 | 20 | 25 | 30 | 35 | 40 | 45 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RcT | 0 | 0.5 | 0.6 | **1.0** | 0.6 | 0.2 | 0.05 | 0 |

**Crop-age response `RcA`**: 1.0 for the first 10 days, declining to 0.01 by
day 90. Full 25-point table in `epirice_model.R`.

### Deviation from epicrop: RcT peak at 25°C

The `epicrop` package source carries a code comment: *"The optimum temperature
for leaf blast as presented in Table 2 of Savary et al. 2012 has a typo. The
optimal value should be 20°C, not 25°C as shown."* That comment corrects the
table to a 20°C peak.

**This implementation uses the published 25°C peak for the following reasons:**

1. Empirical infection-rate data directly support a 24–25°C optimum. Hashioka
   (1965), used in the Lanoiselet *et al.* (2002) DYMEX model for the Australian
   rice belt, measured the minimum time required for conidial germination and
   penetration: 6 h at 24°C, 8 h at 28°C, 10 h at 32°C. Converting to a rate
   (1/hours), infection is fastest at 24°C and declines on both sides. This is
   the same data used to parameterise the temperature slope in DYMEX and it
   places the infection optimum squarely at 24–25°C, not 20°C.

2. The broader biological literature is consistent. The optimum for spore
   germination, infection, lesion formation and sporulation is reported at 25–28°C
   across multiple sources (Barksdale & Jones 1965; UC IPM California rice blast
   guide; Advances in Rice Blast, ScienceDirect 2025).

3. The consequence for northern Australia is large. At 28°C daily mean (typical
   North Queensland wet season), the 20°C curve gives RcT = 0.36 and the 25°C
   curve gives RcT = 0.76: approximately a 2-fold difference in predicted
   infection rate. Using the 20°C peak would systematically under-predict risk
   at exactly the temperatures and location this tool is built for.

4. The epicrop note is a package author's reading of the paper, not a published
   erratum. It may reflect the original cropsim parameterisation for temperate
   Philippines/Japan rather than a correction to the biological literature.

Note that sporulation peaks at a cooler temperature than infection. Kato and
Kozaka (1974) measured 399 spores per day per lesion at 20°C but only 271 at
25°C and 131 at 32°C. EPIRICE uses a single RcT for both infection and disease
development, so the chosen peak is a compromise. The 25°C value better represents
the infection step, which is the mechanistic basis for the SEIR transition.

To restore the epicrop 20°C version, change the `RcT` coefficient vector in
`epirice_model.R` from `c(0, 0.5, 0.6, 1, 0.6, 0.2, 0.05, 0)` to
`c(0, 0.5, 1, 0.6, 0.2, 0.05, 0.01, 0)`.

### Two differences from a stock epicrop run

1. `duration` is set to `CROP_AGE_DAYS` (60) rather than 120, so each run
   reports "disease built up in a 60-day-old crop", not a full season.
2. Weather is Open-Meteo ERA5 rather than NASA POWER. Column names are mapped
   in `openmeteo_wth.R`.

### Why EPIRICE reads near zero in July, and what that means

The wetness gate (`rhlim = 90%` daily mean) is a much higher bar than BLASTAM's
hourly RH threshold. A day where nights reach 100% RH but afternoons drop to
60% will have a daily mean of perhaps 75–80% and will not open the EPIRICE
wetness gate. BLASTAM detects those nights correctly. Both models are behaving;
they answer different questions.

There is a second, independent reason for low July readings: July is midwinter in
the southern rice belt and early dry season in the tropical north. The
Lanoiselet *et al.* (2002) DYMEX model, run over October–April at four southern
NSW locations, found conditions favourable in 10 of 11 seasons — but that model
ran only during the rice-growing season when temperatures and humidity are both
higher. A flat blue map in July is biologically plausible.

---

## BLASTAM: parameters and deviations

`blastam_model.R` implements the infection-warning model of Koshimizu (1988),
operated across Japan on the AMeDAS automated weather network.

### How it works

For each night the model asks three questions:

1. Was the leaf wetness duration long enough?
2. Was the temperature during the wetness period within the favourable range?
3. Has the background temperature of the past five days been warm enough?

If all three say yes, the night is **favourable**. If the wetness threshold is
met but exactly one temperature criterion fails, the night is
**semi-favourable**. Results are stored per day as `infect` (0/1/NA) and `semi`
(0/1/NA); NA means the night could not be judged (incomplete hourly series or
insufficient lead-in). The BLASTAM score reported in the email and the map is
the count of favourable nights over the past `BLASTAM_WINDOW_DAYS` (21) days.

### Original criteria (Koshimizu 1988)

| Criterion | Original value |
| --- | --- |
| Leaf wetness duration | ≥ 10 hours (fixed) |
| Temperature during wetness | 15–25°C |
| Preceding 5-day mean temperature | 20–25°C |
| Night window | 15:00 to 09:00, local Japan Standard Time |
| Leaf wetness estimation | Energy balance from AMeDAS sunshine, wind, rainfall |
| Heavy rain | Hours ≥ 4 mm/h excluded (Yoshino 1988, via Hayashi & Koshimizu) |

### Our five deviations, and why

**Deviation 1 — Temperature-dependent wetness threshold (replaces fixed 10 h)**

Koshimizu's 10-hour threshold was calibrated empirically for Tohoku, where
wetness-period temperatures are typically 15–20°C. The biology shows the required
wetness duration decreases with temperature (Barksdale & Jones 1965; Kato 1974):

| Wetness-period temp | Required hours |
| --- | --- |
| 15.6°C (60°F) | 12.2 h |
| 18.3°C (65°F) | 10.9 h |
| 21.1°C (70°F) | 9.7 h |
| 23.9°C (75°F) | 8.6 h |
| 26.7°C (80°F) | 7.7 h |

*Source: Barksdale & Jones (1965), lower 95% confidence interval.*

For Tohoku the 10 h approximation is reasonable. For Australia it biases in both
directions: at 17°C nights the required threshold is 11.5 h but we were scoring
against 10 h (over-counting); at 26°C nights the threshold is 7.9 h but we were
demanding 10 h (under-counting by roughly 41% of qualifying nights in the
measured cache data).

The Barksdale & Jones curve is implemented as `blastam_bj_min_hours()`, a linear
interpolation with flat extrapolation outside the range. Set
`BLASTAM_USE_BJ_THRESHOLD <- FALSE` to restore the original fixed 10 h.

**Deviation 2 — Upper temperature bounds raised for tropical Australia**

| Bound | Koshimizu 1988 | This implementation |
| --- | --- | --- |
| Wetness-period mean, upper | 25°C | **32°C** |
| Preceding 5-day mean, upper | 25°C | **30°C** |
| Lower bounds (both) | unchanged | unchanged |

The Japanese 25°C caps were calibrated for temperate Tohoku (July–August mean
temperatures 23–25°C). They exclude the entire tropical wet-season temperature
range (26–32°C), where blast is most active. The shortest required wetness
duration for infection occurs near 25–28°C (Barksdale & Jones 1965; Kato 1974),
exactly where the original upper bound cuts off. These upper bounds are
deliberately raised; the lower bounds are not changed.

To restore the original bounds: `BLASTAM_TWET_MAX <- 25` and
`BLASTAM_PREV5_MAX <- 25` in `blast_config.R` or `blastam_model.R`.

**Deviation 3 — Heavy rain exclusion stated explicitly**

Hours with hourly precipitation ≥ `BLASTAM_RAIN_HEAVY` (4 mm/h, matching
Yoshino 1988) are excluded from the infection-conducive wet count. Heavy rain
washes conidia off the leaf surface, making those hours unfavourable for
infection even though the leaf is physically wet. Koshimizu's energy-balance
estimator handled this implicitly; our RH proxy does not, so it must be stated.
Set `BLASTAM_RAIN_HEAVY <- Inf` to disable.

**Deviation 4 — Leaf wetness estimated from hourly RH**

| Original | This implementation |
| --- | --- |
| Energy balance from AMeDAS sunshine, wind and rain | RH ≥ 90% or rain ≥ 0.2 mm/h, excluding rain ≥ 4 mm/h (Deviation 3) |

ERA5 provides no measured leaf-wetness variable. RH ≥ 90% is the standard
humidity proxy for leaf wetness. The 0.2 mm/h rain threshold picks up drizzle
that the RH sensor may lag. Mean temperature during the wet period (`temp_wet`)
is the mean during infection-conducive wet hours, which differs from the AMeDAS
energy-balance approach. Absolute wet-hour counts are provisional until
calibrated against field leaf-wetness data.

The Lanoiselet *et al.* (2002) DYMEX model required 100% RH for infection, based
on the free-water requirement for conidial germination (Hemmi & Imura 1939; Ou
1985). That threshold appears very strict but was applied to a simulated sine-wave
diurnal cycle, so a day where the minimum RH was 60% and the maximum was 100%
would still produce infection-conducive hours. The BLASTAM hourly RH ≥ 90%
approach, applied to real hourly ERA5 observations, is functionally similar: both
capture the wet portion of the diurnal cycle. Neither uses the published 100%
free-water threshold directly.

**Deviation 5 — Local solar time rather than political timezone**

The 15:00–09:00 window is applied in local solar time (longitude / 15 h) rather
than Japan Standard Time or the local political timezone. Political zone
boundaries jump at state lines and would introduce a step discontinuity in
wet-hour counts across the continental map. Solar time follows the sun, which
governs dew formation and leaf drying.

### Completeness and lead-in

A night can only be judged when both the previous evening (hours ≥ 15:00 local)
and the following morning (hours < 09:00 local) are present, and when five
preceding days of temperature exist for the 5-day mean. Where either is
unavailable, `infect` and `semi` are stored as NA rather than 0.

Callers fetch `BLASTAM_LEADIN_DAYS` (6) extra days and discard them: one day
because the first local solar day is partial, and five days for the preceding
5-day mean. The 8-day refresh tail (`REFRESH_TAIL_DAYS`) plus the 6-day lead-in
totals exactly 14 days, which is the API's minimum charge unit, giving a refresh
cost of exactly 1.0 weighted calls per point.

### BLASTAM scores in the email

The email reports two numbers per town. The first (`BLASTAM days`) is the count
of **favourable** nights in the last 21 days. The second in parentheses
(`7d N`) is the count of favourable nights in the last 7 days, which shows
whether the trend is building. Semi-favourable nights are tracked separately and
reported in `blast_results_latest.csv` but not in the email table.

---

## Known limitations

### Ambient versus in-canopy humidity

Both EPIRICE and BLASTAM are driven by ERA5 ambient air humidity, which
corresponds roughly to a Stevenson screen at field-bank height. Lanoiselet *et
al.* (2002) placed data loggers directly in the rice canopy at Yanco during the
2000–01 season and found in-canopy RH averaged at least 20 percentage points
higher than the ambient reading, with the difference influenced by wind speed,
rainfall and evapotranspiration. The offset is large enough to be practically
significant: at an ambient ERA5 hourly RH of 72%, the in-canopy value would
exceed the BLASTAM 90% threshold. Both EPIRICE's daily mean gate and BLASTAM's
hourly wetness count are therefore likely to under-count infection-conducive
conditions in irrigated paddocks, and model outputs should be read as
conservative lower bounds rather than direct estimates.

Quantifying this offset requires simultaneous in-canopy loggers and ERA5-driven
model runs at the same site. The Yanco data collected by Lanoiselet may still
exist in CSIRO or NSW Agriculture records and would provide a starting point.

### No validation against field outbreaks

Australia remains free of blast in cultivated rice, so direct field validation is
not possible here. Lanoiselet *et al.* (2002) validated their DYMEX model against
the California 1996–1999 outbreak (first recorded blast in California) with
reasonable skill: the model predicted infection events at the outbreak site and
relatively few in a nearby disease-free area, though timing did not match
perfectly. A similar validation of BLASTAM and EPIRICE against documented outbreaks
in climatically similar regions (northern Japan, the Sacramento Valley) would
strengthen confidence in the parameterisation. The 25°C RcT choice is best tested
against tropical incidence data from the Philippines or northern Thailand.

### Rolling emergence and seasonal interpretation

Emergence is computed as `end_date − CROP_AGE_DAYS` and moves with every run.
The EPIRICE output is therefore "disease that would accumulate in a 60-day-old
crop given current weather", not a seasonal total. The BLASTAM window of 21 days
is similarly rolling. Neither output can be read as a season-to-date progression.
For season-total analysis, a fixed emergence date would be needed.

---

## Workflow

The Monday workflow (`.github/workflows/weekly_blast.yml`) runs:

1. **Offline tests** — `Rscript test_offline.R`. 20 tests, no network. Fails the
   run before three hours of fetching are spent on broken code.
2. **Continental heatmaps** — `Rscript run_blast_grid.R`. Fetches weather, runs
   both models on the grid, renders two PNG heatmaps and two GeoTIFFs.
3. **Town table** — `Rscript run_blast.R`. Fetches weather for the 31 monitoring
   sites, runs both models, writes the CSV, text and HTML summary.
4. **Commit** — pushes the updated cache, trends CSVs, failure ledger, version
   file and map stats. The cache is committed first so a failed email step does
   not lose the fetched data.
5. **Email** — `python3 send_email.py`. Sends the HTML summary with four
   attachments (two heatmap PNGs, two trends CSVs) via Gmail SMTP SSL.
6. **Upload artifact** — saves all `blast_outputs/` to the Actions artifact
   store for 90 days as a fallback if email fails.

Two seasonal cron entries bracket the daylight saving change:

```yaml
- cron: '30 20 * * 0'   # 06:30 Monday AEST (winter, UTC+10)
- cron: '30 19 * * 0'   # 06:30 Monday AEDT (summer, UTC+11)
```

The gate job matches the fired schedule against the current UTC offset, so a
late-firing cron (GitHub routinely fires 5–30 minutes late) does not skip the
week.

A `concurrency` group prevents two simultaneous runs from clobbering the cache.

---

## Attribution

Please retain the attribution comments in `epirice_model.R` and
`blastam_model.R`, and cite the primary publications if you publish or
distribute results.

**EPIRICE:**

Savary, S., Nelson, A., Willocquet, L., Pangga, I., and Aunario, J. (2012).
Modeling and mapping potential epidemics of rice diseases globally. *Crop
Protection* 34: 6–17. doi:10.1016/j.cropro.2011.11.009

epicrop package (vendored SEIR engine and leaf blast parameters): Adam H. Sparks
and colleagues (Department of Primary Industries and Regional Development, WA;
IRRI). Check the epicrop licence before redistributing the model code. Framework:
Zadoks, J.C. (1971). A formal definition of host–pathogen interaction in
epidemiological terms. *Phytopathology* 61: 600–610.

Hashioka, Y. (1965). Effects of environmental factors on development of causal
fungus, infection, disease development, and epidemiology in rice blast disease.
In: *The Rice Blast Disease*. J Hopkins Press: Baltimore, pp. 153–161. [empirical
basis for the 24–25°C infection-rate optimum used in the RcT curve]

**BLASTAM:**

Koshimizu, Y. (1988). A forecasting method for occurrence of rice leaf blast with
AMeDAS data. *Bulletin of the Tohoku National Agricultural Experiment Station* 78:
67–121. [in Japanese]

Hayashi, T. and Koshimizu, Y. (1988). Computer program BLASTAM for forecasting
occurrence of rice leaf blast. *Bulletin of the Tohoku National Agricultural
Experiment Station* 78: 123–138. [in Japanese]

Barksdale, T.H. and Jones, M.W. (1965). Minimum conditions of temperature and
dew period for infection of rice by *Piricularia oryzae*. *Phytopathology* 55:
1037–1040. [basis for the temperature-dependent wetness threshold, Deviation 1]

Kato, H. (1974). Epidemiology of blast. *Review of Plant Protection Research*
7: 1–20. [temperature–wetness interaction, Deviations 1 and 2]

Maehara, H. and Yamada, M. (2025). Annual changes in the timing and frequency of
favorable conditions for rice leaf blast infection estimated by BLASTAM in
Fukushima Prefecture. *Annual Report of the Society of Plant Protection of North
Japan* 76: 41–46. doi:10.11455/kitanihon.2025.76\_41

**Prior Australian modelling:**

Lanoiselet, V., Cother, E.J. and Ash, G.J. (2002). CLIMEX and DYMEX simulations
of the potential occurrence of rice blast disease in south-eastern Australia.
*Australasian Plant Pathology* 31: 1–7. doi:10.1071/AP01070

This is the earliest peer-reviewed Australian application of computational blast
risk modelling, using a DYMEX population model validated against the California
1996–1999 outbreak. Key data on in-canopy humidity and temperature–infection
relationships from that study directly inform model choices here. The current
system extends that work from four BOM point locations in the NSW rice belt to a
continental operational grid.

**Weather:**

Open-Meteo historical ERA5 archive, data licensed CC BY 4.0. Non-commercial
research use. Hersbach, H. *et al.* (2020). The ERA5 global reanalysis.
*Quarterly Journal of the Royal Meteorological Society* 146: 1999–2049.
doi:10.1002/qj.3803

ERA5 native resolution is 0.25 deg (~28 km). A 0.3 deg model lattice matches
this driver; finer grids buy interpolation rather than information unless
`OPENMETEO_MODEL` is switched to `"era5_land"` (0.1 deg, ~11 km). Confirm that
`relative_humidity_2m` is served by ERA5-Land before switching, and restate the
map resolution in any publications if you do.

# Weekly blast risk (EPIRICE + BLASTAM, Open-Meteo)

A small, self contained pipeline that runs two complementary blast models each
week from GitHub Actions, using free Open-Meteo weather, and writes two risk
maps, a town table and a summary, and emails them out.

## What this is

Two models, run from one weather fetch:

- **EPIRICE** (Savary *et al.* 2012), a mechanistic SEIR epidemic simulation.
  From an emergence date it steps through the season day by day, tracking
  healthy, latent, infectious and removed leaf sites, and reports leaf blast
  **intensity** (the proportion of leaf tissue diseased). It answers: how much
  disease has the season built up? The SEIR engine and leaf blast parameters are
  vendored into `epirice_model.R` so no package install is needed.
- **BLASTAM** (`blastam_model.R`), the Japanese infection-warning model of
  Koshimizu (1988). For each day it judges whether conditions favoured a new
  infection, and counts the favourable **infection days in the last 21 days** (`BLASTAM_WINDOW_DAYS`). A day is favourable
  when leaf wetness is at least 10 hours, the mean temperature during wetness is
  within bounds, and the preceding 5-day mean temperature is within bounds. The
  original Japanese bounds are 15-25 C and 20-25 C; here the upper bounds are
  raised to 15-32 C and 20-30 C for warmer northern-Australian conditions (a
  deliberate, documented deviation, since blast infects up to about 32 C given
  adequate leaf wetness). Leaf wetness is estimated from hourly humidity (>=90%)
  and rainfall, since ERA5 has no measured leaf-wetness variable. See
  `blastam_model.R` for parameters and refs.

The two are complementary: BLASTAM flags when infection windows open (early
warning), while EPIRICE estimates the disease that may follow. Each run produces
a map and a rolling trends CSV for each model.

`openmeteo_wth.R` fetches the Open-Meteo historical archive. A single hourly
fetch per point feeds BLASTAM (leaf-wetness nights) and is aggregated to the
daily fields EPIRICE expects, so both models share the same weather at no extra
API cost.
A day counts as wet, and therefore infection favourable, when mean relative
humidity is at or above 90 per cent or rainfall is at or above 5 mm, which are
the published thresholds.

## Files

| File | Role |
| --- | --- |
| `epirice_model.R` | Vendored SEIR engine and exact leaf blast parameters |
| `openmeteo_wth.R` | Open-Meteo weather adapter (replaces NASA POWER) |
| `blast_config.R` | Your sites, season, risk bands, output and email settings |
| `run_blast.R` | Runner: fetch weather, run model, write map, table, summary |
| `run_blast_grid.R` | Continental risk **heatmap**: runs the model on a grid |
| `.github/workflows/weekly_blast.yml` | Weekly GitHub Actions schedule |

## Two outputs

`run_blast.R` gives the site level view: a table, a summary and a point map for
your monitoring locations. `run_blast_grid.R` gives a continental **heatmap**: it
runs EPIRICE on a grid of points and renders a continuous risk surface as a PNG
and a GeoTIFF. Both are driven by the same Open-Meteo weather and the same model.

### About the heatmap

The heatmap follows the approach of Savary *et al.* (2012), who mapped potential
epidemics on a grid. Each week it assumes a crop of age `CROP_AGE_DAYS` at every
land cell and colours the modelled leaf blast intensity. It covers all land,
including coastal cells, and is clipped to the coastline so the ocean is blank.
Coastline, rivers, roads, towns and your monitoring sites are drawn on top. Read
it as a potential risk surface: it shades the land as if rice were grown
everywhere, so most of it is a weather driven potential, not actual crop or
measured disease.

Colours run from light blue (low) to dark red (high) on a FIXED scale set by
`HEAT_MAX`, so a given colour means the same intensity every week. That value is
provisional: EPIRICE intensity is near zero in the dry season and rises in the
wet, so run once through a favourable month with `HEAT_MAX <- NULL` to read the
peak, then set `HEAT_MAX` a little above it. Until then, dry-season maps read
low, which is the correct signal.

Output files are date stamped, for example `epirice_heatmap_2026-07-10.png` and `blastam_heatmap_2026-07-10.png` and the
matching GeoTIFF, with a `_latest.png` copy of each for convenience.

Grid settings live in `blast_config.R`:

- Dynamic resolution with a weather cache. Each run keeps a weather cache
  (`weather_cache.csv.gz`, committed to the repo, pruned to the 60 day window).
  Points already in the cache only need their newest days fetched, which is cheap,
  so the spare API budget is spent ADDING new points. The map therefore fills
  coarse to fine over successive runs, reaching about 0.5 degree (~2,800 land
  points over the continent) in roughly six weeks, then holds there. This keeps
  every run inside the free Open-Meteo hourly limit (5,000 weighted calls).
- `GRID_RES_FINEST`: the finest resolution the map refines toward (default 0.5).
- `GRID_RES_LEVELS`: the coarse-to-fine fill order (each a whole multiple of
  `GRID_RES_FINEST`).
- `TARGET_CALLS_PER_RUN`: weighted-call budget per run (default 4,000, under the
  5,000 per hour cap). Larger fills faster but risks rate limits.
- `GRID_EXTENT`: the mapped area. Default is the Australian continent.
- `LAND_ONLY`: keeps land points and clips the surface to the coast, using the
  bundled `australia_land.geojson` (read with terra, no `sf`).
- `SHOW_COAST`, `SHOW_TOWNS`, `SHOW_ROADS`, `SHOW_RIVERS`: overlay toggles. Coast,
  rivers and roads come from the bundled GeoJSON files. `MONITOR_TOWNS` is the
  editable list of towns, highlighted on the map and tracked in the trends CSV.
- `CROP_AGE_DAYS`: the trailing crop age and weather window (default 60).
- `HEAT_MAX`: fixes the colour scale (deepest red at this %); `NULL` auto-scales.

A note on magnitude: EPIRICE intensity values are often small in absolute terms,
especially early in a crop. The heatmap is most useful for the spatial pattern
and its change over the season, rather than the absolute percentage.

## Setup

1. Put these files in a new repository.
2. Edit `blast_config.R`:
   - `SITES`: name, latitude, longitude for each monitoring location.
   - `EMERGENCE_DATE`: transplanting or emergence date for the current crop.
     Optionally add an `emergence` column to `SITES` for per site dates.
   - Risk band cut points, if you have grounds to change them.
3. Commit and push. In the repository, open Settings, Actions, General, and
   under Workflow permissions choose read and write.
4. Trigger a first run from the Actions tab (Run workflow), or wait for the
   Monday schedule.

Each run writes to `blast_outputs/`: a dated map, table and summary, plus
`*_latest` copies, and commits them back to the repository.

## Run it locally first

```bash
Rscript -e 'install.packages(c("data.table","jsonlite"))'
Rscript run_blast.R
```

The Open-Meteo archive lags real time by roughly five days, so the run reports
up to about five days ago. That is set by `ARCHIVE_LAG_DAYS` in the config.

## Email (required)

Email is how the two heatmaps leave the runner: the PNGs are NOT committed to the
repo, so email is their delivery path (with the run artifact as a 90-day fallback
copy). Email is sent by the workflow, not by R, which avoids a Java dependency.

Add three repository secrets: `MAIL_USERNAME` (a Gmail address), `MAIL_PASSWORD`
(a Gmail app password, not the account password), and `MAIL_TO` (recipient, or a
comma separated list). If any is missing the email step fails loudly (a red run)
rather than skipping silently, so a delivery problem is visible. The maps are
still recoverable for 90 days from the run's uploaded artifact.

Each email carries: the combined town table (EPIRICE and BLASTAM), the two
heatmap PNGs, and the two trends CSVs.

## Reading the output

Intensity is the modelled proportion of diseased leaf sites. The risk bands
(low, moderate, high) are provisional cut points in the config. They only change
the labelling and the map colours, not the model. Before relying on the bands,
compare a season or two of modelled intensity against your own field
observations and adjust the cut points so the labels match what you see. The
`trend7` figure in the summary is the change in intensity over the last seven
days, which indicates whether the epidemic is building or easing.

## Scope and limits

- The model estimates unmanaged disease risk from weather. It does not account
  for fungicide programmes, cultivar resistance, or nitrogen status, all of
  which change actual disease.
- Weather is ERA5 reanalysis on a grid of roughly 10 to 25 km, not a station in
  your canopy, so it is a regional signal rather than a paddock measurement.
- The wet day proxy uses daily mean humidity and rainfall. It does not resolve
  the length of overnight leaf wetness directly.

## Attribution

Please retain the attribution in `epirice_model.R` and `blastam_model.R`, and
cite the models if you publish or distribute results.

EPIRICE:

Savary, S., Nelson, A., Willocquet, L., Pangga, I., and Aunario, J. (2012).
Modeling and mapping potential epidemics of rice diseases globally. *Crop
Protection* 34: 6 to 17. doi:10.1016/j.cropro.2011.11.009.

epicrop package (vendored SEIR engine and leaf blast parameters): Adam H. Sparks
and colleagues. Check the epicrop licence before redistributing the model code
itself. Framework: Zadoks, J.C. (1971). *Phytopathology* 61: 600 to 610.

BLASTAM:

Koshimizu, Y. (1988). A forecasting method for occurrence of rice leaf blast
with AMeDAS data. *Bulletin of the Tohoku National Agricultural Experiment
Station* 78: 67 to 121. [in Japanese]

Hayashi, T. and Koshimizu, Y. (1988). Computer program BLASTAM for forecasting
occurrence of rice leaf blast. *Bulletin of the Tohoku National Agricultural
Experiment Station* 78: 123 to 138. [in Japanese]

Maehara, H. and Yamada, M. (2025). Annual changes in the timing and frequency of
favorable conditions for rice leaf blast infection estimated by BLASTAM in
Fukushima Prefecture. *Annual Report of the Society of Plant Protection of North
Japan* 76: 41 to 46. doi:10.11455/kitanihon.2025.76_41.

The infection criteria follow the operational Japanese BLASTAM (leaf wetness at
least 10 hours; mean temperature during wetness and preceding 5-day mean within
bounds), with the upper temperature bounds deliberately raised from the Japanese
25 C to 32 C (wetness period) and 30 C (5-day mean) for warmer
northern-Australian conditions. Set both maxima back to 25 in `blastam_model.R`
to reproduce the original. Leaf wetness is estimated from hourly humidity and
rainfall, since ERA5 has no measured leaf-wetness variable.

Weather:

Open-Meteo historical (ERA5) archive, data licensed CC BY 4.0.
Hersbach, H. et al. (2020). The ERA5 global reanalysis. *Quarterly Journal of
the Royal Meteorological Society* 146: 1999 to 2049. doi:10.1002/qj.3803.

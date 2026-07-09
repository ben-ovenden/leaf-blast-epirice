# Weekly rice leaf blast risk (EPIRICE + Open-Meteo)

A small, self contained pipeline that runs the EPIRICE leaf blast model each
week from GitHub Actions, using free Open-Meteo weather, and writes a risk map,
a table and a summary. Email delivery is optional.

## What this is

The disease model is the EPIRICE leaf blast model of Savary *et al.* (2012),
implemented in the epicrop R package by Adam Sparks and colleagues. The core
SEIR engine and the exact leaf blast parameters are vendored into
`epirice_model.R` so the pipeline runs without installing the package. The only
substitution is the weather source: `openmeteo_wth.R` replaces the package's
NASA POWER downloader with the Open-Meteo historical archive, returning the same
daily fields the model expects (mean temperature, mean relative humidity, and
rainfall).

The model is a mechanistic epidemic simulation. From an emergence date it steps
through the season day by day, tracking healthy, latent, infectious and removed
leaf sites, and reports leaf blast intensity (the proportion of diseased sites).
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
grid cell and colours the modelled leaf blast intensity. Read it as a potential
risk surface: it shades the whole extent as if rice were grown everywhere, so
most of the continent is a weather driven potential, not actual crop or measured
disease. The monitoring sites are drawn on top for reference.

Grid settings live in `blast_config.R`:

- `GRID_EXTENT` and `GRID_RES`: area and resolution. The default is the continent
  at 0.75 degree with the ocean masked out, about 1,240 land cells, which uses
  roughly 80 per cent of the free Open-Meteo allowance (10,000 calls per day) for
  a 60 day window. The runner prints an estimated call budget each time and warns
  if a change pushes it over. Finer than about 0.70 degree over the whole
  continent exceeds the free tier, so go finer only over a smaller extent.
- `LAND_ONLY`: drops ocean cells before fetching using an Australia polygon (needs
  `ozmaps`). This spends no budget on water and roughly halves the cell count.
- `CROP_AGE_DAYS`: the trailing crop age and weather window. Larger values give a
  fuller season potential but cost more requests.
- `SMOOTH_FACTOR`: display only smoothing; the model still runs at `GRID_RES`.
- `HEAT_MAX`: fixes the colour scale for week to week comparability, or leave
  `NULL` to auto scale each week.

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

## Email (optional)

Email is sent by the workflow, not by R, which avoids a Java dependency. Add
three repository secrets: `MAIL_USERNAME` (a Gmail address), `MAIL_PASSWORD` (a
Gmail app password, not the account password), and `MAIL_TO` (recipient, or a
comma separated list). If the secrets are absent the email step is skipped.

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

Please retain the attribution in `epirice_model.R` and cite the model if you
publish or distribute results:

Savary, S., Nelson, A., Willocquet, L., Pangga, I., and Aunario, J. (2012).
Modeling and mapping potential epidemics of rice diseases globally. *Crop
Protection* 34: 6 to 17. doi:10.1016/j.cropro.2011.11.009.

epicrop package: Adam H. Sparks and colleagues. Check the epicrop licence
before redistributing the model code itself.

# IFS trial

A second, separate run of the leaf blast pipeline on ECMWF IFS weather. It is
sent beside the regular (ERA5) email from November 2026 so the maintainer can
compare the two and keep one.

## Completely separate from the product

The product is the set of scripts in the repository root. **Nothing in this
directory is shared with it, and nothing in the root was changed for it.**

| | Product (repository root) | Trial (this directory) |
| --- | --- | --- |
| Scripts, config, tests | the root's own files | its own copy of every one, here |
| Weather | ERA5, 0.25 degrees, 6 days behind | ECMWF IFS, 9 km, to the newest settled day |
| Outputs | `blast_outputs/` | `ifs/blast_outputs_ifs/` |
| Cache branch | `cache-data` | `cache-data-ifs` |
| Workflows | `weekly_blast.yml`, `midweek_topup.yml` | `ifs_weekly.yml`, `ifs_topup.yml` |
| Grid fetched (UTC days) | Sunday; Wednesday, Saturday | Monday, Tuesday, Thursday, Friday |
| Weekly email run | Sunday UTC | Sunday UTC (towns only; map from the cache) |
| Concurrency group | `blast-grid` | `blast-grid-ifs` |
| Email | Monday by 7 am, towns and maps to the previous Sunday | Monday by 7 am, towns to Saturday, maps to Thursday, subject `[IFS trial] ...` |

The two workflows run every step with `working-directory: ifs`. The root scripts
contain no switch and no mention of the trial; test 31 here fails if that
changes. What the two still share is the repository itself, its `main` branch
(each commits its own run state there, in its own directory), the mail secrets,
and the weather API's free daily allowance (see the next section).

The first version of the trial (15 November 2026, PR 15) was a
`BLAST_SOURCE=ifs` switch inside the product's scripts. The maintainer asked the
same day for it to be completely separate, so the product's files were put back
exactly as they were at the tag `era5-stable-2026-11-15` and the trial became
this directory (PR 16).

## The schedule, and why the map is two days older than the town table

The free weather API allows 10,000 weighted calls a day. Refreshing the whole map
grid (7,721 cells) costs about 7,700, and the product does that on Sunday,
Wednesday and Saturday UTC. So the trial cannot refresh its own grid on those
days without the two together going over the allowance.

- **Monday, Tuesday, Thursday, Friday UTC** (`ifs_topup.yml`, 07:07): the trial
  fetches its grid, a whole day's budget each time.
- **Sunday UTC** (`ifs_weekly.yml`, 06:47): the email run. It fetches only the 31
  towns (150 to 270 calls, which fit beside the product's run even at the
  product's own ceiling of 9,500) and draws the map from the cache as Friday's
  top-up left it. `map_window_edge.R` names that day for the grid step.

GitHub starts this repository's schedules five to seven hours late, so the
Sunday run starts at about midday to 2 pm UTC and the email arrives at about
midnight to 2 am on Monday, Sydney time.

The result, in the Monday email: **town table to Saturday, maps to Thursday**.
The health check allows the map to be two days behind the table
(`HEALTH_MAX_MAP_BEHIND_DAYS`); three means Friday's top-up did not run, and the
email says so. The subject carries both dates.

To have the maps to Saturday as well, the trial's grid would have to be fetched
on Sunday, which needs either the product's Sunday allowance (if IFS is kept and
ERA5 switched off) or a paid API plan.

## The latest weather that has happened

The IFS archive is stitched from forecast runs started at 00, 06, 12 and 18 UTC.
Each reaches the archive about six hours after it starts, and the archive then
runs on for days into the future, so "has data" does not mean "has happened".

A UTC day may be the last day fetched once a forecast run started on that day
has landed: the only hours of it the models use are its first three, and those
are then a same-day run's own first hours. The archive publishes which run it
holds, and `probe_archive_edge.R` asks. If it cannot tell, or the day's first
run has not landed after the wait, the window ends a day earlier. A run that
has not asked at all (a local run) never uses today.

So a Sunday run models to Saturday, where the product, six days behind, models
to the Sunday before.

## Filling the grid, and when the second email starts

From an empty cache every cell costs 4.86 weighted calls. While any cell remains
to be added the grid runner spends its budget on new cells first
(`IFS_ADD_FIRST`) and refreshes cached ones only with what is left; the
product's order, refresh first, is right for a full grid and takes about twelve
fetch days to fill an empty one. Add-first takes about five, then one more to
bring every cell to the same day.

**The second email starts only when the data is sufficient**: the map complete
at 0.30 degrees AND the run healthy (every town modelled, the map current,
almost no cell grey). The run that first qualifies writes
`blast_outputs_ifs/email_started.txt`, and every weekly run after that emails,
as the product does. Until then the weekly run still models the towns, draws
what map it can and keeps the results as a workflow artifact, and an incomplete
map is reported as a notice, not a failed run.

From the first scheduled top-up on Monday 16 November 2026: grid full on about
Monday 23 November, every cell on one day by Tuesday 24, and the first email on
**Monday 30 November**, after the first Friday top-up of a full grid.

It goes to the `MAIL_TO` recipients unless a `MAIL_TO_IFS` secret names a
different list.

## What differs from the product's code

Seven files are the product's, line for line: `blastam_model.R`,
`epirice_model.R`, `grid_window.R`, `map_labels.R`, `openmeteo_batch.R`,
`openmeteo_wth.R`, `run_health.R`, and the three `.geojson` overlays. The test
suite prints a note if they drift apart; a comparison of two weather sources is
only that while the models are the same.

The rest differ only as follows.

- `blast_config.R`, section 6b: the IFS archive (`ecmwf_ifs`), the window rule
  above, the output directory, add-first, the map and top-up age limits, and the
  wording that names the source (`WEATHER_NAME`, the "Trial version" banner).
- `run_blast_grid.R`: add-first while the grid is filling; "IFS:" in map titles.
- `run_blast.R`: prints the trial's wording (email header and banner).
- `probe_archive_edge.R`: asks the archive which forecast run it holds.
- `map_window_edge.R`: new, the trial's only; the day the weekly map is drawn at.
- `send_email.py`: subject starts `[IFS trial]`, attachments are named `ifs_...`.
- Every entry script stops unless the config it loaded is the trial's, so none
  can be run against the product's directory by mistake.
- `test_offline.R`: the product's suite with the expectations that differ for
  IFS changed, plus test 31 (the separation, the schedule, the fill).

## Running it

Always from this directory:

```
cd ifs
Rscript test_offline.R
```

The weather cache is not on `main`. To get it into a local checkout, see the
comment in `.gitignore` here.

A capped top-up proves the workflow in about three minutes and thirty calls:
`gh workflow run ifs_topup.yml --ref main -f weighted_cap=30`.

## The two are not interchangeable

Over the same 60 days to 7 November 2026 at the 31 towns, IFS had 190 days of
5 mm or more against ERA5's 144, and 163 infection nights against 144; the
BLASTAM count at Home Hill was 10 against 1. The gap is in the tropics and on
the coast; in NSW the two agree. Which is closer to the ground has not been
established, which is what the trial is for. The trial's trends files are a
separate record and must not be joined to the product's.

## Ending the trial

- **Keep ERA5:** delete this directory and the two `ifs_*.yml` workflows; the
  `cache-data-ifs` branch can then be deleted too. Nothing in the root needs
  changing.
- **Keep IFS:** disable the two ERA5 workflows. Sunday's allowance is then free,
  so the weekly run can fetch its own grid again (remove the map-window step and
  the zero budget from `ifs_weekly.yml`) and the maps reach Saturday too.

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
| Weather | ERA5, 0.25 degrees, 6 days behind | ECMWF IFS, 9 km, held back 2 days |
| Outputs | `blast_outputs/` | `ifs/blast_outputs_ifs/` |
| Cache branch | `cache-data` | `cache-data-ifs` |
| Workflows | `weekly_blast.yml`, `midweek_topup.yml` | `ifs_weekly.yml`, `ifs_topup.yml` |
| UTC days | Sunday; Wednesday, Saturday | Monday; Tuesday, Thursday, Friday |
| Concurrency group | `blast-grid` | `blast-grid-ifs` |
| Email | Monday by 7 am, weather to the previous Sunday | Monday evening, weather to Friday, subject `[IFS trial] ...` |

The two workflows run every step with `working-directory: ifs`. The root scripts
contain no switch and no mention of the trial; test 31 here fails if that
changes. What the two still share is the repository itself, its `main` branch
(each commits its own run state there, in its own directory) and the mail
secrets.

The first version of the trial (15 November 2026, PR 15) was a
`BLAST_SOURCE=ifs` switch inside the product's scripts. The maintainer asked the
same day for it to be completely separate, so the product's files were put back
exactly as they were at the tag `era5-stable-2026-11-15` and the trial became
this directory.

## What differs from the product's code

Seven files are the product's, line for line: `blastam_model.R`,
`epirice_model.R`, `grid_window.R`, `map_labels.R`, `openmeteo_batch.R`,
`openmeteo_wth.R`, `run_health.R`, and the three `.geojson` overlays. The test
suite prints a note if they drift apart; a comparison of two weather sources is
only that while the models are the same.

The rest differ only as follows.

- `blast_config.R`, section 6b: the IFS archive (`ecmwf_ifs`), a two day lag,
  the output directory, the top-up age limit, and the wording that names the
  source (`WEATHER_NAME`, the "Trial version" banner).
- `run_blast.R`, `run_blast_grid.R`: print that wording (email header and
  banner, "IFS:" in map titles).
- `send_email.py`: subject starts `[IFS trial]`, attachments are named `ifs_...`.
- `run_blast.R`, `run_blast_grid.R`, `probe_archive_edge.R`: stop unless the
  config they loaded is the trial's, so they cannot be run against the
  product's directory by mistake.
- `test_offline.R`: the product's suite with the expectations that differ for
  IFS changed, plus test 31 (the separation).

## Running it

Always from this directory:

```
cd ifs
Rscript test_offline.R
```

The weather cache is not on `main`. To get it into a local checkout, see the
comment in `.gitignore` here.

## Why the days differ

Each grid run spends most of a day's free API allowance, so the two must not run
on the same UTC day. The product keeps its days; the trial takes the other four
and fires early in the UTC day (00:17) so that a late start cannot spill into a
product day. That is why the trial's email arrives on Monday evening. If IFS is
kept, it takes over the Sunday slot and the 7 am delivery.

## Why IFS is held back two days

Its archive is updated every six hours and runs to the end of today, so its
newest hours are forecast (on 14 November 2026 at 02:41 UTC it returned every
hour of that day). Two UTC days keeps every modelled hour on settled data.

## Accumulation, and when the second email starts

The trial started with an empty cache and fills it coarse to fine, about twelve
fetch days: at four a week, about three weeks. Its weekly run models the towns
and draws whatever map it can from the first week, and keeps the results as a
workflow artifact, but **the second email starts only when the map is complete
at 0.30 degrees**. The run that first reaches it writes
`blast_outputs_ifs/email_started.txt`, and every weekly run after that emails.
It goes to the `MAIL_TO` recipients unless a `MAIL_TO_IFS` secret names a
different list.

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
- **Keep IFS:** disable the two ERA5 workflows and move `ifs_weekly.yml` to the
  Sunday slot (06:17 UTC), so its email arrives by 7 am Monday.

#!/usr/bin/env python3
"""Send the weekly blast summary email (HTML body + attachments) over SMTP.

This is a Node-independent replacement for the send-mail GitHub Action, which
stopped delivering once GitHub forced JavaScript actions onto Node 24. It uses
only the Python standard library, so no pip install is needed.

Reads from the environment: MAIL_USERNAME, MAIL_PASSWORD (Gmail app password),
MAIL_TO (comma-separated), and the run date. Exits non-zero with a clear
::error:: message if anything required is missing, so a failure is loud in the
Actions log rather than silent.

The run date is taken from BLAST_RUN_DATE or RUN_DATE, both exported by the
workflow's "Resolve run date" step, and falls back to blast_outputs/run_date.txt.
All three now agree by construction; previously the R scripts each called
Sys.Date() separately and a run straddling local midnight produced a map and a
table dated differently.
"""
import os
import sys
import ssl
import smtplib
import mimetypes
from email.message import EmailMessage

OUT = "blast_outputs"
SMTP_HOST = "smtp.gmail.com"
SMTP_PORT = 465
FROM_NAME = "WWAI Cereal Pathology: blast models"


def require(name):
    value = os.environ.get(name, "").strip()
    if not value:
        sys.exit(f"::error::{name} is empty; cannot send email.")
    return value


def resolve_run_date():
    for key in ("BLAST_RUN_DATE", "RUN_DATE"):
        v = os.environ.get(key, "").strip()
        if v:
            return v
    path = os.path.join(OUT, "run_date.txt")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            v = f.read().strip()
            if v:
                return v
    sys.exit("::error::no run date in BLAST_RUN_DATE, RUN_DATE or run_date.txt.")


def read_first_line(path):
    if not os.path.exists(path):
        return ""
    with open(path, encoding="utf-8") as f:
        return f.readline().strip()


def read_run_status(path):
    """key=value lines from run_blast.R; {} if the file is absent."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf-8") as f:
        for line in f:
            if "=" in line:
                k, v = line.rstrip("\n").split("=", 1)
                out[k.strip()] = v.strip()
    return out


def build_subject(run_date, status, map_end):
    """'[DEGRADED] Blast risk summary 2026-09-21 (weather to 2026-09-14; maps to 2026-08-29)'

    The window in the subject is the TOWN table's, the core product; the map's
    window is added only when it differs. The subject used to carry the map
    window alone, so three September 2026 emails said "weather to 2026-08-29"
    above a town table modelled to a later date, and nothing in the subject said
    the run was degraded.
    """
    town_end = status.get("town_window_end", "")
    window = town_end or map_end
    subject = f"Blast risk summary {run_date}"
    if window:
        subject += f" (weather to {window}"
        if town_end and map_end and map_end != town_end:
            subject += f"; maps to {map_end}"
        subject += ")"
    if status.get("degraded") == "1":
        subject = "[DEGRADED] " + subject
    return subject


def attach(msg, path):
    ctype, _ = mimetypes.guess_type(path)
    maintype, subtype = (ctype or "application/octet-stream").split("/", 1)
    with open(path, "rb") as f:
        msg.add_attachment(f.read(), maintype=maintype, subtype=subtype,
                           filename=os.path.basename(path))


def main():
    user = require("MAIL_USERNAME")
    password = require("MAIL_PASSWORD")
    recipients = [a.strip() for a in require("MAIL_TO").split(",") if a.strip()]
    run_date = resolve_run_date()

    # Field 8 of map_stats.txt is the date the maps were modelled to. The town
    # table's window and the run's health verdict come from run_status.txt,
    # written by run_blast.R (see run_health.R).
    stats = read_first_line(os.path.join(OUT, "map_stats.txt")).split("|")
    map_end = stats[7] if len(stats) > 7 and stats[7] else ""
    status = read_run_status(os.path.join(OUT, "run_status.txt"))
    if status.get("degraded") == "1":
        print(f"::warning title=Degraded blast run::{status.get('reasons', '')}")

    msg = EmailMessage()
    msg["Subject"] = build_subject(run_date, status, map_end)
    msg["From"] = f"{FROM_NAME} <{user}>"
    msg["To"] = ", ".join(recipients)

    txt_path = os.path.join(OUT, "blast_summary_latest.txt")
    html_path = os.path.join(OUT, "blast_summary_latest.html")
    for p in (txt_path, html_path):
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            sys.exit(f"::error::email body missing or empty: {p}")
    with open(txt_path, encoding="utf-8") as f:
        msg.set_content(f.read())
    with open(html_path, encoding="utf-8") as f:
        msg.add_alternative(f.read(), subtype="html")

    # The town table and the trends are the core product and must be present.
    required = [
        os.path.join(OUT, "town_trends.csv"),
        os.path.join(OUT, "blastam_trends.csv"),
    ]
    # The heatmaps are OPTIONAL, deliberately. A run that could not refresh the
    # grid, typically because the daily weather-API quota was already spent, may
    # legitimately have too few current cells to render. Failing the whole email
    # on a missing PNG means the one run that most needs explaining is the one
    # nobody hears about: the 2026-07-31 run died here with 1,874 usable cached
    # points sitting in the repository. The R side writes a "Degraded run"
    # banner into the body naming exactly what is missing and why.
    #
    # run_log.csv records the schema version and model parameters behind each
    # trends column. Optional: an older cache will not have one yet.
    optional = [
        os.path.join(OUT, f"epirice_heatmap_{run_date}.png"),
        os.path.join(OUT, f"blastam_heatmap_{run_date}.png"),
        os.path.join(OUT, "run_log.csv"),
    ]

    n = 0
    for p in required:
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            sys.exit(f"::error::attachment missing or empty: {p} (run date {run_date})")
        attach(msg, p)
        n += 1
    skipped = []
    for p in optional:
        if os.path.exists(p) and os.path.getsize(p) > 0:
            attach(msg, p)
            n += 1
        else:
            skipped.append(os.path.basename(p))
    if skipped:
        # A warning, not an error: the email still goes out and explains itself.
        print(f"::warning::not attached (absent or empty): {', '.join(skipped)}")

    context = ssl.create_default_context()
    with smtplib.SMTP_SSL(SMTP_HOST, SMTP_PORT, context=context, timeout=60) as server:
        server.login(user, password)
        server.send_message(msg)

    print(f"Email sent to {', '.join(recipients)} with {n} attachments "
          f"(run date {run_date}, weather to {window_end or 'unknown'})."
          + (f" Skipped: {', '.join(skipped)}." if skipped else ""))


if __name__ == "__main__":
    main()

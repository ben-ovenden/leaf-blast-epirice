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

    # Field 8 of map_stats.txt is the date the maps were modelled to. Putting it
    # in the subject line makes the data window visible without opening the mail,
    # and makes a stale run obvious.
    stats = read_first_line(os.path.join(OUT, "map_stats.txt")).split("|")
    window_end = stats[7] if len(stats) > 7 and stats[7] else ""

    msg = EmailMessage()
    subject = f"Blast risk summary {run_date}"
    if window_end:
        subject += f" (weather to {window_end})"
    msg["Subject"] = subject
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

    required = [
        os.path.join(OUT, f"epirice_heatmap_{run_date}.png"),
        os.path.join(OUT, f"blastam_heatmap_{run_date}.png"),
        os.path.join(OUT, "town_trends.csv"),
        os.path.join(OUT, "blastam_trends.csv"),
    ]
    # run_log.csv records the schema version and model parameters behind each
    # trends column, so a change of method is not read as a change in the
    # weather. Optional: an older cache will not have one yet.
    optional = [os.path.join(OUT, "run_log.csv")]

    n = 0
    for p in required:
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            sys.exit(f"::error::attachment missing or empty: {p} (run date {run_date})")
        attach(msg, p)
        n += 1
    for p in optional:
        if os.path.exists(p) and os.path.getsize(p) > 0:
            attach(msg, p)
            n += 1

    context = ssl.create_default_context()
    with smtplib.SMTP_SSL(SMTP_HOST, SMTP_PORT, context=context, timeout=60) as server:
        server.login(user, password)
        server.send_message(msg)

    print(f"Email sent to {', '.join(recipients)} with {n} attachments "
          f"(run date {run_date}, weather to {window_end or 'unknown'}).")


if __name__ == "__main__":
    main()

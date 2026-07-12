#!/usr/bin/env python3
"""Send the weekly blast summary email (HTML body + attachments) over SMTP.

This is a Node-independent replacement for the send-mail GitHub Action, which
stopped delivering once GitHub forced JavaScript actions onto Node 24. It uses
only the Python standard library, so no pip install is needed.

Reads from the environment: MAIL_USERNAME, MAIL_PASSWORD (Gmail app password),
MAIL_TO (comma-separated), and RUN_DATE (YYYY-MM-DD, used to pick the dated
heatmap files). Exits non-zero with a clear ::error:: message if anything is
missing, so a failure is loud in the Actions log rather than silent.
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
FROM_NAME = "WWAI Cereal Pathology - Blast Models"


def require(name):
    value = os.environ.get(name, "").strip()
    if not value:
        sys.exit(f"::error::{name} is empty; cannot send email.")
    return value


def main():
    user = require("MAIL_USERNAME")
    password = require("MAIL_PASSWORD")
    recipients = [a.strip() for a in require("MAIL_TO").split(",") if a.strip()]
    run_date = os.environ.get("RUN_DATE", "").strip()

    msg = EmailMessage()
    msg["Subject"] = f"Blast risk summary {run_date}".strip()
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

    attachments = [
        os.path.join(OUT, f"epirice_heatmap_{run_date}.png"),
        os.path.join(OUT, f"blastam_heatmap_{run_date}.png"),
        os.path.join(OUT, "town_trends.csv"),
        os.path.join(OUT, "blastam_trends.csv"),
    ]
    for p in attachments:
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            sys.exit(f"::error::attachment missing or empty: {p} (RUN_DATE={run_date})")
        ctype, _ = mimetypes.guess_type(p)
        maintype, subtype = (ctype or "application/octet-stream").split("/", 1)
        with open(p, "rb") as f:
            msg.add_attachment(f.read(), maintype=maintype, subtype=subtype,
                               filename=os.path.basename(p))

    context = ssl.create_default_context()
    with smtplib.SMTP_SSL(SMTP_HOST, SMTP_PORT, context=context, timeout=60) as server:
        server.login(user, password)
        server.send_message(msg)

    print(f"Email sent to {', '.join(recipients)} with "
          f"{len(attachments)} attachments (RUN_DATE={run_date}).")


if __name__ == "__main__":
    main()

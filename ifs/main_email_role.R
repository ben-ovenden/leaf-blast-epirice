#!/usr/bin/env Rscript
################################################################################
# main_email_role.R   (IFS only)
#
# WHICH EMAIL IS THE MAIN ONE THIS WEEK. Asks GitHub whether the product's weekly
# workflow is switched on, and prints ROLE=trial (it is: this run is the trial
# beside it) or ROLE=main (it is off or gone: this run is the email). The rule
# and the reasons are in section 6b of blast_config.R; ifs_role_from_state() is
# the pure part and is tested offline.
#
# Reads from the environment, all set by GitHub Actions or by the workflow step:
#   GITHUB_REPOSITORY, GITHUB_API_URL   which repository to ask about
#   GH_TOKEN                            the job's own token (actions: read)
#   IFS_PRODUCT_LEDGER                  optional: the product's spend ledger, to
#                                       refuse "main" on a day it has fetched
# Never fails: anything it cannot find out reads as "trial".
################################################################################
SCRIPT_DIR <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile), winslash = "/"),
  error = function(e) normalizePath(getwd(), winslash = "/"))
suppressPackageStartupMessages(library(data.table))
source(file.path(SCRIPT_DIR, "blast_config.R"))
# THE IFS TRIAL'S COPY (see README.md in this directory). It must be run from
# ifs/, where blast_config.R is the trial's; refuse anywhere else.
if (!isTRUE(get0("IFS_TRIAL")))
  stop("this is the IFS trial's copy: run it from the ifs/ directory", call. = FALSE)
say <- function(...) cat(..., "\n", sep = "", file = stderr())

api  <- Sys.getenv("GITHUB_API_URL", "https://api.github.com")
repo <- Sys.getenv("GITHUB_REPOSITORY", "")
tok  <- Sys.getenv("GH_TOKEN", "")
status <- NA_integer_; state <- NA_character_
if (nzchar(repo)) {
  r <- tryCatch({
    h <- curl::new_handle(timeout = 30, connecttimeout = 20, useragent = "curl")
    hdr <- list(Accept = "application/vnd.github+json")
    if (nzchar(tok)) hdr$Authorization <- paste("Bearer", tok)
    do.call(curl::handle_setheaders, c(list(h), hdr))
    curl::curl_fetch_memory(sprintf("%s/repos/%s/actions/workflows/%s", api, repo,
                                    IFS_PRODUCT_WEEKLY_WORKFLOW), handle = h)
  }, error = function(e) NULL)
  if (!is.null(r)) {
    status <- as.integer(r$status_code)
    if (identical(status, 200L))
      state <- tryCatch(as.character(jsonlite::fromJSON(rawToChar(r$content))$state)[1],
                        error = function(e) NA_character_)
  }
}

# What the product has spent today, if the workflow handed over its ledger.
spent <- 0
lf <- Sys.getenv("IFS_PRODUCT_LEDGER", "")
if (nzchar(lf) && file.exists(lf)) spent <- tryCatch({
  d <- fread(lf)
  sum(as.numeric(d[as.character(utc_day) == format(blast_utc_date()), spent]), na.rm = TRUE)
}, error = function(e) 0)

role <- ifs_role_from_state(status, state, spent)
say(sprintf("The regular email's workflow (%s) is '%s' (HTTP %s); it has spent %.0f weighted calls today. This run is the %s.",
            IFS_PRODUCT_WEEKLY_WORKFLOW, if (is.na(state)) "unknown" else state,
            if (is.na(status)) "no answer" else as.character(status), spent,
            if (identical(role, "main")) "MAIN email" else "trial beside it"))
cat("STATE=", if (is.na(state)) "unknown" else state, "\n", sep = "")
cat("ROLE=", role, "\n", sep = "")

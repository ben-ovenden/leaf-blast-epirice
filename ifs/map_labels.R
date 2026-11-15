################################################################################
# map_labels.R
#
# Where each town label goes on the heatmaps. A pure function of the marker
# positions and the label sizes, so it can be tested without a graphics device.
#
# WHY. Every label used to go on one fixed side of its marker (right, or left
# near the eastern edge), and the only defence against collisions was to DROP any
# label whose marker was within LABEL_MIN_SEP_DEG of one already placed. That
# tested marker to marker distance, not what the text covers: Humpty Doo and
# Jabiru are 1.6 degrees apart on the same latitude, so both labels were kept and
# "Humpty Doo" ran straight through Jabiru's marker and name; Kununurra did the
# same to Timber Creek, and Goondiwindi to Warwick. Meanwhile five towns that
# could have been labelled on another side were not labelled at all.
#
# Each label now tries twelve positions around its marker and takes the first
# that stays inside the plot, touches no other label and covers no town marker.
#
#   Pass 1  most constrained first: the label with the fewest positions still
#           open goes next, which is what lets the crowded clusters resolve.
#   Pass 2  repair: a label left with nowhere to go may still fit if exactly ONE
#           label already placed steps aside to another free position.
#
# A label that still has no position is suppressed and counted; its marker is
# drawn regardless.
################################################################################

# Candidate positions: name, direction of the offset, and text adjustment. The
# order is the order of preference: beside the marker first, then beside it and
# nudged half a line up or down, then the diagonals, then above and below.
#
# The four nudged positions are there for robustness, not looks. Warwick sits
# 0.59 degrees of latitude from Lismore, which is almost exactly a label's half
# height plus a marker's radius, so whether "Warwick" fitted beside its marker
# came down to a pixel of font metrics: it did on one of two maps drawn in the
# same run and not on the other. Fonts differ again between Windows and the Linux
# runner. A test places the labels at 85%, 100% and 115% of the measured size.
.LABEL_POS <- data.frame(
  where = c("E", "W", "ENE", "ESE", "WNW", "WSW", "NE", "SE", "NW", "SW", "N", "S"),
  dx    = c( 1,  -1,  0.95,  0.95, -0.95, -0.95, 0.7,  0.7, -0.7, -0.7,  0,   0),
  dy    = c( 0,   0,  0.45, -0.45,  0.45, -0.45, 0.7, -0.7,  0.7, -0.7,  1,  -1),
  adjx  = c( 0,   1,  0,     0,     1,     1,    0,    0,    1,    1,    0.5, 0.5),
  adjy  = c( 0.5, 0.5, 0.2,  0.8,   0.2,   0.8,  0,    1,    0,    1,    0,   1),
  stringsAsFactors = FALSE)

#   x, y        marker positions (map units)
#   w, h        label width and height, same units as x and y respectively
#   offx, offy  gap between the marker centre and the text, in x and y units
#   xlim, ylim  the plot region; a label may not leave it
#   mrx, mry    half-size of a town marker
#   pad         extra clearance between labels, as a fraction of label height
#   prefer      positions to try, in order
#
# Returns a data.frame with one row per town: tx, ty (where to draw the text),
# adjx, adjy (its adjustment) and where (the position chosen, NA if suppressed).
place_labels <- function(x, y, w, h, offx, offy, xlim, ylim,
                         mrx = offx * 0.6, mry = offy * 0.6, pad = 0.15,
                         prefer = .LABEL_POS$where) {
  n <- length(x)
  out <- data.frame(tx = rep(NA_real_, n), ty = NA_real_, adjx = NA_real_, adjy = NA_real_,
                    where = NA_character_, stringsAsFactors = FALSE)
  if (n == 0L) return(out)
  w <- rep_len(w, n); h <- rep_len(h, n)
  cand <- .LABEL_POS[match(prefer, .LABEL_POS$where), , drop = FALSE]
  cand <- cand[!is.na(cand$where), , drop = FALSE]
  nk <- nrow(cand)

  # c(x0, x1, y0, y1, tx, ty): the box a label would occupy at candidate k, and
  # the point its text is drawn at.
  box <- function(i, k) {
    tx <- x[i] + cand$dx[k] * offx; ty <- y[i] + cand$dy[k] * offy
    c(tx - cand$adjx[k] * w[i], tx + (1 - cand$adjx[k]) * w[i],
      ty - cand$adjy[k] * h[i], ty + (1 - cand$adjy[k]) * h[i], tx, ty)
  }
  hit <- function(a, b, px = 0, py = 0)
    a[1] < b[2] + px && b[1] < a[2] + px && a[3] < b[4] + py && b[3] < a[4] + py
  markers <- lapply(seq_len(n), function(j) c(x[j] - mrx, x[j] + mrx, y[j] - mry, y[j] + mry))

  boxes <- matrix(NA_real_, n, 4L)      # the box each placed label occupies
  kpos  <- rep(NA_integer_, n)          # the candidate each placed label took
  inside <- function(b) b[1] >= xlim[1] && b[2] <= xlim[2] && b[3] >= ylim[1] && b[4] <= ylim[2]
  on_marker <- function(b, i)
    any(vapply(seq_len(n)[-i], function(j) hit(b, markers[[j]]), logical(1)))
  blockers <- function(b, i) {
    js <- setdiff(which(!is.na(kpos)), i)
    js[vapply(js, function(j) hit(b, boxes[j, ], pad * h[i], pad * h[i]), logical(1))]
  }
  free <- function(i, k) {
    b <- box(i, k)
    inside(b) && !on_marker(b, i) && length(blockers(b, i)) == 0L
  }
  put <- function(i, k) { boxes[i, ] <<- box(i, k)[1:4]; kpos[i] <<- k }

  # Pass 1: most constrained first. Ties go south to north, then west to east, so
  # the result does not depend on the order the towns are listed in.
  todo <- seq_len(n)
  while (length(todo) > 0L) {
    opts <- lapply(todo, function(i) vapply(seq_len(nk), function(k) free(i, k), logical(1)))
    left <- vapply(opts, sum, integer(1))
    pick <- order(left, y[todo], x[todo])[1]
    if (left[pick] > 0L) put(todo[pick], which(opts[[pick]])[1])
    todo <- todo[-pick]
  }

  # Pass 2: repair. For a label with nowhere to go, look for a position blocked by
  # exactly one placed label that can itself move to another free position.
  stuck <- which(is.na(kpos))
  for (i in stuck[order(y[stuck], x[stuck])]) {
    moved <- FALSE
    for (k in seq_len(nk)) {
      b <- box(i, k)
      if (!inside(b) || on_marker(b, i)) next
      j <- blockers(b, i)
      if (length(j) != 1L) next
      for (k2 in setdiff(seq_len(nk), kpos[j])) {
        b2 <- box(j, k2)
        if (!inside(b2) || on_marker(b2, j) || length(blockers(b2, j)) > 0L) next
        if (hit(b2, b, pad * h[j], pad * h[j])) next
        put(j, k2); put(i, k); moved <- TRUE
        break
      }
      if (moved) break
    }
  }

  for (i in which(!is.na(kpos))) {
    b <- box(i, kpos[i])
    out$tx[i] <- b[5]; out$ty[i] <- b[6]
    out$adjx[i] <- cand$adjx[kpos[i]]; out$adjy[i] <- cand$adjy[kpos[i]]
    out$where[i] <- cand$where[kpos[i]]
  }
  out
}

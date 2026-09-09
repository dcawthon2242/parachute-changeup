## Batter silhouette for catcher's-view strike-zone plots.
##
## Follows the approach in dcawthon2242/CSUF-Pitch-Trajectory-Plots: drop in a batter PNG
## and place it beside the plate, flipping it horizontally for the other box. Put the
## artwork at assets/BatterPerspective.png (or point BATTER_PNG at it) and every figure
## picks it up automatically. Without that file we fall back to a drawn silhouette.
##
## These plots look down the pitcher-home axis, so the batter -- who stands beside the
## plate with his chest facing it -- is seen in profile. Geometry is in feet with the
## middle of the plate at x = 0 and the ground at z = 0.
##
## Side convention, verified against 2026 hit-by-pitches (n = 1,420) in the mirrored
## frame these scripts use (+x = pitcher's glove side): a same-handed hitter sits at mean
## x = -1.89 ft and an opposite-handed hitter at +2.07 ft. Those are contact points on the
## near edge of the body, so the centreline sits a little further out.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(grid) })

BATTER_X      <- 2.35    # body centreline (feet), ft from the middle of the plate
BATTER_FILL   <- "grey76"
BATTER_TOP    <- 6.10    # ft from the ground to the top of the artwork (the bat tip)

.batter_png_path <- function() {
  cand <- c(Sys.getenv("BATTER_PNG", ""),
            file.path("assets", "BatterPerspective.png"),
            file.path("data", "assets", "BatterPerspective.png"))
  cand <- cand[nzchar(cand)]
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else NA_character_
}

# Read the artwork, reduce it to an alpha mask, trim the margins, and note where the feet
# sit across the width so the figure can be pinned to the batter's box rather than to the
# image's centre (the bat throws that centre off). Works with or without a real alpha
# channel: art exported on a white background gets its mask derived from luminance.
.batter_art <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    pth <- .batter_png_path()
    if (is.na(pth) || !requireNamespace("png", quietly = TRUE)) return(NA)
    img <- png::readPNG(pth); d <- dim(img)
    a <- if (d[3] == 4 && diff(range(img[, , 4])) > 0.1) {
      img[, , 4]
    } else {
      lum <- if (d[3] >= 3) (img[,,1] + img[,,2] + img[,,3])/3 else img[,,1]
      matrix(pmin(1, pmax(0, (0.90 - lum) / 0.30)), nrow = nrow(lum))
    }
    rr <- which(rowSums(a > 0.5) > 0); cc <- which(colSums(a > 0.5) > 0)
    if (!length(rr) || !length(cc)) return(NA)
    a <- a[min(rr):max(rr), min(cc):max(cc), drop = FALSE]
    H <- nrow(a); W <- ncol(a)
    foot <- a[round(0.90*H):H, , drop = FALSE]
    fc <- sum(col(foot) * (foot > 0.5)) / sum(foot > 0.5)
    cache <<- list(alpha = a, aspect = W/H, foot_frac = (fc - 1)/max(W - 1, 1))
    cache
  }
})

# Recolour to a flat BATTER_FILL and mirror for the opposite box, as the CSUF script does.
.batter_grob <- function(side) {
  art <- .batter_art(); if (!is.list(art)) return(NULL)
  a <- art$alpha
  if (side > 0) a <- a[, rev(seq_len(ncol(a))), drop = FALSE]
  rgb0 <- grDevices::col2rgb(BATTER_FILL) / 255
  arr <- array(0, c(nrow(a), ncol(a), 4))
  arr[, , 1] <- rgb0[1]; arr[, , 2] <- rgb0[2]; arr[, , 3] <- rgb0[3]; arr[, , 4] <- a
  # fill the placement box exactly; the box is built to the art's own aspect, so nothing
  # is stretched, but leaving these NULL would letterbox and float the feet off the ground
  rasterGrob(arr, interpolate = TRUE, width = unit(1, "npc"), height = unit(1, "npc"))
}

## ---- drawn fallback ---------------------------------------------------------
# Local coords: u = horizontal offset from the centreline, positive TOWARD the plate;
# v = height. Control points for a ~6'0" hitter in a loaded stance, traced as one closed
# loop -- up the back, over the helmet, down the chest and shin, then back along the
# inside of the legs. Points repeated 3x stay sharp through the corner-cutting below.
.rep3 <- function(u, v) list(u = rep(u, 3), v = rep(v, 3))
.CTRL <- local({
  pts <- list(
    list(-0.40, 0.00, TRUE ),  # back heel
    list(-0.36, 0.22, FALSE),  list(-0.33, 0.62, FALSE),  # achilles, calf
    list(-0.28, 1.10, FALSE),  list(-0.27, 1.35, FALSE),  # knee
    list(-0.33, 1.75, FALSE),  list(-0.40, 2.20, FALSE),  # hamstring, glute
    list(-0.41, 2.55, FALSE),  list(-0.35, 2.95, FALSE),  # waist taper
    list(-0.38, 3.45, FALSE),  list(-0.43, 3.95, FALSE),  # back, lat
    list(-0.40, 4.30, FALSE),  list(-0.28, 4.55, FALSE),  # rear delt, trap
    list(-0.15, 4.68, FALSE),  list(-0.20, 4.95, FALSE),  # neck, helmet back
    list(-0.13, 5.25, FALSE),  list( 0.07, 5.34, FALSE),  # helmet crown
    list( 0.24, 5.16, FALSE),  list( 0.28, 4.92, FALSE),  # helmet front
    list( 0.50, 4.88, TRUE ),  list( 0.48, 4.78, TRUE ),  # brim tip
    list( 0.29, 4.76, FALSE),  list( 0.24, 4.60, FALSE),  # cheek, jaw
    list( 0.20, 4.48, FALSE),  list( 0.30, 4.38, FALSE),  # neck, front delt
    list( 0.33, 4.00, FALSE),  list( 0.27, 3.45, FALSE),  # chest, ribs
    list( 0.21, 2.95, FALSE),  list( 0.26, 2.45, FALSE),  # waist, hip
    list( 0.25, 1.95, FALSE),  list( 0.21, 1.35, FALSE),  # quad, knee
    list( 0.19, 1.05, FALSE),  list( 0.22, 0.55, FALSE),  # shin
    list( 0.25, 0.14, FALSE),  list( 0.44, 0.00, TRUE ),  # front toe
    list( 0.12, 0.00, TRUE ),  list( 0.02, 0.85, FALSE),  # inner front, crotch
    list(-0.12, 0.00, TRUE ))                             # inner back
  u <- unlist(lapply(pts, function(p) if (p[[3]]) rep(p[[1]], 3) else p[[1]]))
  v <- unlist(lapply(pts, function(p) if (p[[3]]) rep(p[[2]], 3) else p[[2]]))
  data.table(u = u, v = v)
})

# Chaikin corner cutting: rounds a closed outline without the overshoot a spline gives.
.chaikin <- function(u, v, iter = 3) {
  for (k in seq_len(iter)) {
    j <- c(seq_along(u)[-1], 1L)
    u <- as.vector(rbind(0.75*u + 0.25*u[j], 0.25*u + 0.75*u[j]))
    v <- as.vector(rbind(0.75*v + 0.25*v[j], 0.25*v + 0.75*v[j]))
  }
  data.table(u = u, v = v)
}

.thicken <- function(u, v, w) {
  du <- diff(u); dv <- diff(v)
  seg <- cbind(du, dv) / sqrt(du^2 + dv^2)
  dir <- rbind(seg[1, ], (seg[-1, , drop = FALSE] + seg[-nrow(seg), , drop = FALSE]) / 2,
               seg[nrow(seg), ])
  dir <- dir / sqrt(rowSums(dir^2))
  nx <- -dir[, 2] * w/2; nz <- dir[, 1] * w/2
  data.table(u = c(u + nx, rev(u - nx)), v = c(v + nz, rev(v - nz)))
}

.PARTS <- local({
  body <- .chaikin(.CTRL$u, .CTRL$v, iter = 3)
  arms <- .thicken(c( 0.16, -0.10, -0.50), c(4.30, 3.86, 4.42), w = 0.16)
  bat  <- .thicken(c(-0.48, -0.78),        c(4.44, 5.46),       w = 0.085)
  th   <- seq(0, 2*pi, length.out = 28)
  knob <- data.table(u = -0.48 + 0.085*cos(th), v = 4.44 + 0.085*sin(th))
  list(body = body, arms = arms, bat = bat, knob = knob)
})

#' Silhouette polygons for a batter standing on one side of the plate.
#' @param side -1 for the negative-x box (same-handed in the mirrored frame), +1 for the
#'   positive-x box (opposite-handed).
batter_silhouette <- function(side) {
  cx <- BATTER_X * side
  rbindlist(lapply(names(.PARTS), function(nm) {
    p <- .PARTS[[nm]]
    # the batter faces the plate, so local +u points opposite the side he stands on
    data.table(x = cx - side * p$u, z = p$v,
               grp = paste0(ifelse(side < 0, "L", "R"), "_", nm))
  }))
}

#' ggplot layers placing a batter on one or both sides of the plate.
#' @param specs list of lists, each with `side` (-1/+1) and optional `data`, a one-row
#'   data.table of facet keys restricting the batter to those panels.
#' @return list of layers, usable directly inside a ggplot() chain.
batter_layers <- function(specs) {
  lapply(specs, function(s) {
    g <- .batter_grob(s$side)
    if (!is.null(g)) {
      art <- .batter_art()
      cx <- BATTER_X * s$side
      wft <- BATTER_TOP * art$aspect
      # pin the feet, not the image centre; the mirrored copy has its feet mirrored too
      ff <- if (s$side > 0) 1 - art$foot_frac else art$foot_frac
      layer(data = if (is.null(s$data)) data.frame(x = NA) else s$data,
            stat = StatIdentity, position = PositionIdentity, geom = GeomCustomAnn,
            inherit.aes = FALSE, show.legend = FALSE,
            params = list(grob = g, xmin = cx - ff*wft, xmax = cx + (1-ff)*wft,
                          ymin = 0, ymax = BATTER_TOP))
    } else {
      poly <- batter_silhouette(s$side)
      if (!is.null(s$data)) poly <- cbind(poly, s$data[rep(1, nrow(poly))])
      geom_polygon(data = poly, inherit.aes = FALSE, aes(x, z, group = grp),
                   fill = BATTER_FILL, colour = NA)
    }
  })
}

batter_png_in_use <- function() !is.na(.batter_png_path())

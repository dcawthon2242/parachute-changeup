#!/usr/bin/env Rscript

# WITHIN A GRADE BAND, WHAT SEPARATES THE ARMS THAT BEAT IT FROM THE ARMS THAT DO NOT?
#
# The quintile table in the last run was whiff, and it split on the stuff model's own predicted
# whiff rate. That split is the useful control: inside a quintile every changeup has roughly the
# same raw stuff, so any trait that still predicts the residual is describing deception rather than
# quality. A trait that only looks good because good pitchers have it will flatten out here.
#
# Run separately for whiff and for run value, on their own predicted quintiles, because the two
# have already disagreed once - the matched-axis group misses by eleven whiff points and by nothing
# at all in runs.
#
# Two views of the same thing:
#   CONTRAST  top third versus bottom third of the residual inside each quintile, trait by trait
#   SLOPE     standardised regression of residual on each trait inside each quintile, so the
#             direction and size are comparable across traits and bands
#
# The question that matters more than any single coefficient is whether the signs hold across
# quintiles. A trait that helps low-stuff changeups and hurts high-stuff ones is a different
# finding from one that helps everywhere, and only the second is worth putting in a model.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
set.seed(11); options(width = 230); MDIR <- "data/statcast_model"

OUT <- readRDS(file.path(MDIR, "absorb_preds.rds"))
P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(ax) & is.finite(az) & is.finite(sax) & is.finite(cax) & is.finite(speed_diff) &
       is.finite(axis_diff) & is.finite(arm_angle) & is.finite(arm_diff) &
       is.finite(release_extension) & is.finite(release_pos_z) & is.finite(release_spin_rate)]
P[, lh := p_throws == "L"]
P[, `:=`(mx = fifelse(lh, -ax, ax), mz = az + 32.174, sep = -speed_diff)]
P[, spin_axis := (atan2(sax, cax) * 180/pi) %% 360]
P[, spin_axis := fifelse(lh, (360 - spin_axis) %% 360, spin_axis)]
wrap <- function(d) ((d + 180) %% 360) - 180
P[, dev_raw := wrap(atan2(mx, mz) * 180/pi + (spin_axis - 180))]
P[, dev := wrap(dev_raw - median(dev_raw, na.rm = TRUE))]
P[, id := as.character(pitcher)]

TR <- P[, .(velo = mean(release_speed), sep = mean(sep), axis = mean(axis_diff),
            seamdev = mean(abs(dev)), signed = mean(dev), spin = mean(release_spin_rate),
            slot = mean(arm_angle), armgap = mean(arm_diff),
            ivb = mean(mz)/32.174*12, hb = mean(mx)/32.174*12,
            ext = mean(release_extension), relz = mean(release_pos_z),
            bu = mean(release_spin_rate)/mean(release_speed)), by = .(id, season)]
NM <- unique(P[, .(id, season, nm = player_name)])
TR <- merge(TR, NM, by = c("id","season"))

VARS <- c("velo","sep","axis","seamdev","signed","spin","bu","slot","armgap","ivb","hb","ext","relz")
NICE <- c(velo = "Velocity (mph)", sep = "Separation (mph)", axis = "Spin-axis gap (deg)",
          seamdev = "Seam deviation (deg)", signed = "Signed seam dev (deg)", spin = "Spin rate (rpm)",
          bu = "Bauer units (rpm/mph)", slot = "Arm slot (deg)", armgap = "Arm-angle gap (deg)",
          ivb = "Induced vertical break (in)", hb = "Horizontal break (in)",
          ext = "Extension (ft)", relz = "Release height (ft)")

agg <- function(X) X[, .(n = .N, act = mean(y), pred = mean(p)), by = .(id, season)]

clus_t <- function(D, yv, xv) {                       # slope of y on standardised x, arm-clustered
  D <- D[is.finite(get(yv)) & is.finite(get(xv))]
  if (nrow(D) < 25 || uniqueN(D$id) < 10) return(c(NA, NA, NA))
  x <- scale(D[[xv]])[,1]; y <- D[[yv]]; w <- as.numeric(D$n)
  m <- lm(y ~ x, weights = w)
  u <- residuals(m)*sqrt(w); X <- model.matrix(m)*sqrt(w)
  nc <- uniqueN(D$id); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[2]); s <- unname(sqrt(diag(V))[2]); c(e, s, 2*pt(-abs(e/s), nc-1))
}

run_outcome <- function(key, floor_n, scale, label, unit) {
  A <- agg(OUT$STUFF[[key]])[n >= floor_n]
  A <- merge(A, TR, by = c("id","season"))
  A[, r := scale*(act - pred)]
  A[, qq := cut(pred, quantile(pred, 0:5/5), include.lowest = TRUE, labels = paste0("Q", 1:5))]

  cat(sprintf("\n\n################ %s ################\n", label))
  cat(sprintf("%d pitcher-seasons, %d arms, minimum %d %s. Residual is %s.\n",
              nrow(A), uniqueN(A$id), floor_n, if (key == "w") "swings" else "pitches", unit))
  Q <- A[, .(seasons = .N, pred = mean(scale*pred), actual = mean(scale*act),
             bias = mean(r), sd_resid = sd(r)), by = qq][order(qq)]
  print(Q[, .(quintile = qq, seasons, predicted = round(pred,2), actual = round(actual,2),
              bias = round(bias,2), spread = round(sd_resid,2))], row.names = FALSE)

  ## ---- CONTRAST: top third vs bottom third of the residual, inside each quintile ----
  cat(sprintf("\n--- trait means: top third of the residual minus bottom third, within quintile ---\n"))
  CN <- rbindlist(lapply(levels(A$qq), function(qv) {
    D <- A[qq == qv]; lo <- quantile(D$r, 1/3); hi <- quantile(D$r, 2/3)
    O <- D[r >= hi]; U <- D[r <= lo]
    o <- list(quintile = qv, over = nrow(O), under = nrow(U))
    for (v in VARS) {
      d <- mean(O[[v]]) - mean(U[[v]])
      p <- tryCatch(t.test(O[[v]], U[[v]])$p.value, error = function(e) NA_real_)
      o[[v]] <- d; o[[paste0("p_", v)]] <- p
    }
    as.data.table(o)
  }))
  M <- as.matrix(CN[, ..VARS]); PM <- as.matrix(CN[, paste0("p_", VARS), with = FALSE])
  disp <- matrix(sprintf("%+.2f%s", M, ifelse(PM < .01, "**", ifelse(PM < .05, "*", " "))),
                 nrow = nrow(M), dimnames = list(CN$quintile, VARS))
  print(as.data.frame(disp))
  cat("  (* p<.05, ** p<.01, two-sample t-test; sign is overperformers minus underperformers)\n")

  ## ---- SLOPE: standardised regression inside each quintile, plus overall ----
  cat(sprintf("\n--- residual change per 1 SD of the trait, within quintile (arm-clustered) ---\n"))
  SL <- rbindlist(lapply(c("ALL", levels(A$qq)), function(qv) {
    D <- if (qv == "ALL") A else A[qq == qv]
    o <- list(band = qv, seasons = nrow(D))
    for (v in VARS) { e <- clus_t(D, "r", v); o[[v]] <- e[1]; o[[paste0("p_", v)]] <- e[3] }
    as.data.table(o)
  }))
  M2 <- as.matrix(SL[, ..VARS]); P2 <- as.matrix(SL[, paste0("p_", VARS), with = FALSE])
  d2 <- matrix(sprintf("%+.2f%s", M2, ifelse(P2 < .01, "**", ifelse(P2 < .05, "*", " "))),
               nrow = nrow(M2), dimnames = list(SL$band, VARS))
  print(as.data.frame(d2))

  ## ---- consistency: does the sign hold across all five bands? ----
  QM <- M2[-1, , drop = FALSE]; QP <- P2[-1, , drop = FALSE]
  CONS <- data.table(trait = NICE[VARS], overall = M2[1, ], p_overall = P2[1, ],
                     same_sign = apply(QM, 2, function(z) max(sum(z > 0), sum(z < 0))),
                     n_sig = colSums(QP < .05, na.rm = TRUE),
                     range = apply(QM, 2, function(z) diff(range(z))))
  setorder(CONS, p_overall)
  cat("\n--- which traits work everywhere, not just on average ---\n")
  print(CONS[, .(trait, overall = round(overall,2), p = round(p_overall,4),
                 `bands same sign` = paste0(same_sign, "/5"), `bands sig` = n_sig,
                 `spread across bands` = round(range,2))], row.names = FALSE)
  invisible(A)
}

AW <- run_outcome("w", 75, 100, "WHIFF  -  quintiles of the model's predicted whiff rate",
                  "actual minus predicted whiff, in percentage points")
AR <- run_outcome("r", 200, 100, "RUN VALUE  -  quintiles of the model's predicted run value",
                  "actual minus predicted run value, per 100 pitches, pitcher-positive")

## =============================================================================================
cat("\n\n################ DO THE TWO OUTCOMES AGREE? ################\n\n")
W <- agg(OUT$STUFF$w)[n >= 75][, .(id, season, nw = n, rw = 100*(act - pred))]
R <- agg(OUT$STUFF$r)[n >= 200][, .(id, season, nr = n, rr = 100*(act - pred))]
B <- merge(merge(W, R, by = c("id","season")), TR, by = c("id","season"))
cat(sprintf("%d pitcher-seasons with both. correlation of the two residuals: %+.3f\n\n",
            nrow(B), cor(B$rw, B$rr)))
cat(sprintf("  %-30s %10s %10s %10s\n", "trait", "whiff", "run value", "agree?"))
for (v in VARS) {
  a <- clus_t(copy(B)[, n := nw], "rw", v); b <- clus_t(copy(B)[, n := nr], "rr", v)
  st <- function(e, p) sprintf("%+.2f%s", e, ifelse(p < .01, "**", ifelse(p < .05, "*", " ")))
  ag <- if (is.na(a[1]) || is.na(b[1])) "-" else
        if (a[3] < .05 && b[3] < .05) (if (sign(a[1]) == sign(b[1])) "both, same" else "both, OPPOSITE") else
        if (a[3] < .05) "whiff only" else if (b[3] < .05) "run value only" else ""
  cat(sprintf("  %-30s %10s %10s %10s\n", NICE[[v]], st(a[1], a[3]), st(b[1], b[3]), ag))
}
saveRDS(list(whiff = AW, rv = AR, both = B), file.path(MDIR, "quintile_traits.rds"))

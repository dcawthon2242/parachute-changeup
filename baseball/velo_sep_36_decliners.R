#!/usr/bin/env Rscript

# WHAT THE DOWNGRADED CHANGEUPS HAVE IN COMMON.
#
# Part of the answer is trivially "the model was over-rating them", since the move correlates +0.40
# with the baseline miss. That is the compression doing its job and says nothing about pitch type.
# The interesting question is what remains after that is taken out: given two changeups the baseline
# model over-rated by the same amount, which one gets cut harder?
#
# Four passes.
#   1. Trait means, decliners against risers, standardised so magnitudes compare.
#   2. Multivariate regression of the move on all traits at once, since several are correlated.
#   3. The same regression with the baseline miss included, isolating what the compression keys on
#      beyond the error it is correcting.
#   4. Within-season variability. The compression adds a season AVERAGE, so the arms it should hurt
#      are the ones whose individual pitches read better than their average - inconsistent arsenals
#      where the per-pitch model was being flattered by the good ones.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 215); MDIR <- "data/statcast_model"

D <- readRDS(file.path(MDIR, "smean_movers.rds")); setDT(D)
P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(ax) & is.finite(az) & is.finite(sax) & is.finite(cax) & is.finite(speed_diff) &
       is.finite(axis_diff) & is.finite(arm_angle) & is.finite(arm_diff) &
       is.finite(release_extension) & is.finite(release_pos_z) & is.finite(release_spin_rate)]
P[, lh := p_throws == "L"]
P[, `:=`(mx = fifelse(lh, -ax, ax), mz = az + 32.174, sep = -speed_diff)]
P[, spin_axis := (atan2(sax, cax) * 180/pi) %% 360]
P[, spin_axis := fifelse(lh, (360 - spin_axis) %% 360, spin_axis)]
wrap <- function(d) ((d + 180) %% 360) - 180
P[, dev := wrap(wrap(atan2(mx, mz) * 180/pi + (spin_axis - 180)) -
                median(wrap(atan2(mx, mz) * 180/pi + (spin_axis - 180)), na.rm = TRUE))]
P[, id := as.character(pitcher)]

TR <- P[, .(velo = mean(release_speed), sep = mean(sep), axis = mean(axis_diff),
            seamdev = mean(abs(dev)), signed = mean(dev), spin = mean(release_spin_rate),
            slot = mean(arm_angle), armgap = mean(arm_diff), ivb = mean(mz)/32.174*12,
            hb = mean(mx)/32.174*12, ext = mean(release_extension), relz = mean(release_pos_z),
            # within-season spread of the same measurements, pitch to pitch
            sd_velo = sd(release_speed), sd_sep = sd(sep), sd_axis = sd(axis_diff),
            sd_dev = sd(dev), sd_spin = sd(release_spin_rate), sd_ivb = sd(mz)/32.174*12,
            sd_hb = sd(mx)/32.174*12, sd_relz = sd(release_pos_z)), by = .(id, season)]
D <- merge(D, TR, by = c("id","season"))

LVL <- c("velo","sep","axis","seamdev","signed","spin","slot","armgap","ivb","hb","ext","relz")
# sd_sep is identical to sd_velo by construction (the fastball anchor is a season constant), dropped
VAR <- c("sd_velo","sd_axis","sd_dev","sd_spin","sd_ivb","sd_hb","sd_relz")
NICE <- c(velo = "Velocity (mph)", sep = "Separation (mph)", axis = "Spin-axis gap (deg)",
          seamdev = "Seam deviation (deg)", signed = "Signed seam dev (deg)", spin = "Spin rate (rpm)",
          slot = "Arm slot (deg)", armgap = "Arm-angle gap (deg)", ivb = "Ind. vertical break (in)",
          hb = "Horizontal break (in)", ext = "Extension (ft)", relz = "Release height (ft)",
          sd_velo = "Velocity, pitch to pitch", sd_sep = "Separation, pitch to pitch",
          sd_axis = "Axis gap, pitch to pitch", sd_dev = "Seam dev, pitch to pitch",
          sd_spin = "Spin rate, pitch to pitch", sd_ivb = "IVB, pitch to pitch",
          sd_hb = "HB, pitch to pitch", sd_relz = "Release height, pitch to pitch")

D[, down := mv < 0]
cat(sprintf("%d pitcher-seasons. %d declined (%.0f%%), %d rose. mean move %+.1f down, %+.1f up.\n",
            nrow(D), sum(D$down), 100*mean(D$down), sum(!D$down),
            D[down == TRUE, mean(mv)], D[down == FALSE, mean(mv)]))

## =============================================================================================
cat("\n=== 1. TRAIT MEANS: DECLINERS VERSUS RISERS ===\n\n")
CMP <- rbindlist(lapply(c(LVL, VAR), function(v) {
  a <- D[down == TRUE][[v]]; b <- D[down == FALSE][[v]]
  s <- sd(D[[v]], na.rm = TRUE); tt <- t.test(a, b)
  data.table(trait = NICE[[v]], kind = if (v %in% LVL) "level" else "spread",
             decliners = mean(a, na.rm = TRUE), risers = mean(b, na.rm = TRUE),
             diff = mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE),
             d = (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE))/s, p = tt$p.value)
}))
CMP[, ad := abs(d)]; setorder(CMP, -ad)
print(CMP[, .(trait, kind, decliners = round(decliners,2), risers = round(risers,2),
              diff = round(diff,2), `effect (SD)` = round(d,2), p = signif(p,3))], row.names = FALSE)

## =============================================================================================
cat("\n=== 2. WHAT INDEPENDENTLY DRIVES THE MOVE ===\n\n")
Z <- copy(D); for (v in c(LVL, VAR)) Z[, (v) := scale(get(v))[,1]]
crob <- function(f, dat) {
  environment(f) <- environment()
  dat <- dat[complete.cases(dat[, c(all.vars(f), "id", "nsw"), with = FALSE])]
  w <- as.numeric(dat$nsw); m <- lm(f, dat, weights = w)
  u <- residuals(m)*sqrt(w); X <- model.matrix(m)*sqrt(w)
  nc <- uniqueN(dat$id); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X*u, dat$id)) %*% b * (nc/(nc-1))
  e <- coef(m); s <- sqrt(diag(V))
  data.table(term = names(e), est = unname(e), se = unname(s),
             p = unname(2*pt(-abs(e/s), nc-1)))[term != "(Intercept)"]
}
f1 <- as.formula(paste("mv ~", paste(c(LVL, VAR), collapse = " + ")))
R1 <- crob(f1, Z)[order(p)]
f2 <- as.formula(paste("mv ~ e0 +", paste(c(LVL, VAR), collapse = " + ")))
R2 <- crob(f2, Z)
R <- merge(R1[, .(term, raw = est, p_raw = p)], R2[, .(term, adj = est, p_adj = p)], by = "term")
R[, trait := NICE[term]][, kind := fifelse(term %in% VAR, "spread", "level")]
R <- R[!is.na(trait)]; setorder(R, p_adj)
cat("  index points of move per 1 SD of the trait. 'raw' is traits only; 'adjusted' also\n")
cat("  controls for the baseline miss, so it shows what the compression keys on beyond the\n")
cat("  error it is already correcting. Arm-clustered, swing-weighted.\n\n")
print(R[, .(trait, kind, raw = round(raw,2), p_raw = signif(p_raw,3),
            adjusted = round(adj,2), p_adjusted = signif(p_adj,3))], row.names = FALSE)
cat(sprintf("\n  the baseline miss itself: %+.2f index points per whiff point over-rated (p %s)\n",
            R2[term == "e0", est], signif(R2[term == "e0", p], 3)))

## =============================================================================================
cat("\n=== 3. THE CAREER DECLINERS, TRAIT BY TRAIT AGAINST LEAGUE AVERAGE ===\n\n")
CA <- D[, .(seasons = .N, swings = sum(nsw), move = weighted.mean(mv, nsw)), by = .(id, nm)][seasons >= 3]
worst <- CA[order(move)][1:10]$id; best <- CA[order(-move)][1:10]$id
prof <- function(ids, lab) {
  X <- D[id %in% ids]
  c(list(group = lab, seasons = nrow(X)),
    setNames(lapply(c(LVL, VAR), function(v) mean(X[[v]], na.rm = TRUE)), c(LVL, VAR)))
}
PF <- rbindlist(list(as.data.table(prof(unique(D$id), "league")),
                     as.data.table(prof(worst, "ten biggest career declines")),
                     as.data.table(prof(best, "ten biggest career gains"))))
setnames(PF, c(LVL, VAR), NICE[c(LVL, VAR)])
print(t(PF[, -"seasons"]), quote = FALSE)

cat("\n  the ten arms, individually:\n")
W <- merge(CA[id %in% worst], D[, lapply(.SD, function(z) mean(z)), by = id,
                                .SDcols = c(LVL, "sd_axis","sd_dev","sd_ivb")], by = "id")
print(W[order(move), .(pitcher = nm, seasons, swings, move = round(move,1),
                       velo = round(velo,1), sep = round(sep,1), axis = round(axis,1),
                       seamdev = round(seamdev,1), spin = round(spin), slot = round(slot,1),
                       ivb = round(ivb,2), hb = round(hb,2), ext = round(ext,2),
                       sd_axis = round(sd_axis,1), sd_ivb = round(sd_ivb,2))], row.names = FALSE)

## =============================================================================================
cat("\n=== 4. THE MECHANISM: DOES THE SEASON AVERAGE FLATTER OR PUNISH? ===\n\n")
# a changeup whose per-pitch readings are widely scattered can look good on its best pitches while
# its average is ordinary. quantify by comparing the model's per-pitch spread against the mean.
PS <- P[, .(id, season, dev, axis_diff, sep)][, `:=`(id = as.character(id))]
Q <- D[, .(id, season, mv, nsw, sd_axis, sd_dev, sd_ivb, sd_velo, axis, seamdev)]
Q[, `:=`(cv_axis = sd_axis/pmax(axis, 1), cv_dev = sd_dev/pmax(seamdev, 1))]
LBL <- c(sd_axis = "Axis gap, pitch to pitch", sd_dev = "Seam dev, pitch to pitch",
         sd_ivb = "IVB, pitch to pitch", sd_velo = "Velocity, pitch to pitch",
         cv_axis = "Axis gap, spread / level", cv_dev = "Seam dev, spread / level")
for (v in names(LBL)) {
  z <- Q[[v]]; k <- cut(z, quantile(z, 0:4/4, na.rm = TRUE), include.lowest = TRUE,
                        labels = c("Q1 tightest","Q2","Q3","Q4 loosest"))
  m <- Q[, .(seasons = .N, move = mean(mv)), by = .(k = k)][order(k)]
  cat(sprintf("  %-28s %s\n", LBL[[v]], paste(sprintf("%s %+6.1f", m$k, m$move), collapse = "   ")))
}
cat("\n  correlation of the move with each spread measure:\n")
for (v in VAR) cat(sprintf("    %-30s %+.3f\n", NICE[[v]], cor(D$mv, D[[v]], use = "complete.obs")))

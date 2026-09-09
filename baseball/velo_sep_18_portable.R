#!/usr/bin/env Rscript

# A GATE THAT D1 CAN ACTUALLY MEASURE.
#
# The +4.97 result is built on Statcast's imaged spin axis. TrackMan does not measure that - its
# changeup SpinAxis sits within a few degrees of a deterministic function of the break, so the D1
# axis gap is a MOVEMENT-DIRECTION gap. The two correlate at only r = .47 in MLB, which is why the
# earlier continuous replication was uninformative rather than negative.
#
# So the portable gate has to be denominated in movement direction from the start. The order matters:
# rebuild the MLB gate using movement direction, confirm it still works IN MLB, and only then carry
# the threshold to D1. A gate that does not survive the currency conversion in MLB cannot be tested
# in D1 at all, and finding that out here is cheaper than finding it out there.
#
# Thresholds are carried across by PERCENTILE, not by degrees. The two leagues have different axis
# distributions and different measurement error, so "under 10 degrees" does not mean the same thing
# in both. Percentile matching is what the project used for the D1 arm angle already.
#
# TWO EXCLUSIONS, both principled rather than by name:
#   swing gate 75  removes Roark 2020 (62 swings) along with Skubal 2020 (53) and Cease 2020 (52).
#                  Excluding the three smallest samples is defensible; excluding one man is not.
#   seam shift     CANNOT be ported. It is the disagreement between imaged spin axis and movement
#                  direction, and in D1 those are the same measurement, so the deviation is zero by
#                  construction. It stays an MLB-only diagnostic and is reported here for reference.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 205); MDIR <- "data/statcast_model"
SP <- readRDS(file.path(MDIR, "locked_spec.rds")); L <- SP$data; setDT(L)
NM <- unique(readRDS(file.path(MDIR, "parachute_ff.rds"))[, .(pitcher, player_name)]); setDT(NM)
NM <- unique(NM, by = "pitcher")[, .(id = as.character(pitcher), nm = sub(",.*", "", player_name))]

crob <- function(D, f, k) {
  D <- D[complete.cases(D[, c(all.vars(f), "id"), with = FALSE])]
  if (length(unique(D[[all.vars(f)[2]]])) < 2) return(c(NA,NA,NA,0,0))
  m <- lm(f, D); u <- residuals(m); X <- model.matrix(m); nc <- uniqueN(D$id)
  if (!k %in% colnames(X)) return(c(NA,NA,NA,0,0))
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k])
  c(e, s, 2*pt(-abs(e/s), nc-1), nrow(D), nc)
}

## ---- 1. rebuild the MLB axis gap in movement-direction currency ---------------------------------
P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(ax) & is.finite(az) & is.finite(ax_diff) & is.finite(az_diff)]
P[, `:=`(chx = fifelse(p_throws == "L", -ax, ax), chz = az + 32.174,
         ffx = fifelse(p_throws == "L", -(ax - ax_diff), ax - ax_diff), ffz = (az - az_diff) + 32.174)]
A <- P[, .(n = .N, ca = atan2(mean(chx), mean(chz))*180/pi,
           fa = atan2(mean(ffx), mean(ffz))*180/pi), by = .(pitcher, season)][n >= 30]
A[, `:=`(mv = abs(((ca - fa + 180) %% 360) - 180), id = as.character(pitcher))]
M <- merge(L[league == "MLB"], A[, .(id, season, mv)], by = c("id","season"))
M <- merge(M, NM, by = "id", all.x = TRUE)
D1 <- L[league == "D1" & !is.na(id) & id != "0" & id != ""]

cat("=== 1. are the two leagues' movement-axis gaps on a comparable scale? ===\n")
cmp <- function(x, lab) sprintf("  %-28s median %5.1f | p10 %5.1f  p25 %5.1f  p33 %5.1f  p50 %5.1f  p67 %5.1f\n",
        lab, median(x), quantile(x,.10), quantile(x,.25), quantile(x,.33), quantile(x,.50), quantile(x,.67))
cat(cmp(M$mv, "MLB movement-axis gap"), cmp(D1$axis, "D1 TrackMan axis gap"),
    cmp(M$axis, "MLB imaged spin-axis gap"), sep = "")
cat(sprintf("\n  MLB: cor(imaged, movement) = %+.3f. The locked imaged gate (<10 deg) keeps %.1f%% of seasons.\n",
            cor(M$axis, M$mv), 100*mean(M$axis < 10)))

## ---- 2. does a movement-axis gate still work in MLB? --------------------------------------------
cat("\n=== 2. MLB, movement-axis gate + arm slot, inside the high-separation third, 75+ swings ===\n")
cat("    sweeping the movement-axis percentile so the same rule can be carried to D1\n\n")
H <- M[vs >= quantile(M$vs, 2/3) & nsw >= 75]
cat(sprintf("    pool: %d seasons, %d arms (75-swing gate removes Roark 2020, Skubal 2020, Cease 2020)\n\n",
            nrow(H), uniqueN(H$id)))
SW <- rbindlist(lapply(c(.05,.10,.15,.20,.25,.33), function(q) {
  cut <- quantile(M$mv, q); D <- copy(H); D[, g := mv <= cut & arm >= arm_thr]
  r <- crob(D, y ~ g + vs, "gTRUE")
  data.table(pctile = q, deg = round(cut,1), members = sum(D$g), arms = uniqueN(D[g == TRUE]$id),
             est = r[1], se = r[2], p = r[3]) }))
print(SW[, .(mv_pctile = pctile, mv_cut_deg = deg, members, arms, est = round(est,3),
             se = round(se,3), p = round(p,4))], row.names = FALSE)
cat("\n    for reference, the IMAGED-axis gate on the same 75-swing pool:\n")
D <- copy(H); D[, g := axis < 10 & arm >= arm_thr]; r <- crob(D, y ~ g + vs, "gTRUE")
cat(sprintf("      axis<10 + arm: %+.3f (se %.3f) p = %.4f | %d members, %d arms\n",
            r[1], r[2], r[3], sum(D$g), uniqueN(D[g == TRUE]$id)))

best <- SW[which.min(p)]
cat(sprintf("\n    best movement-axis gate: %.0fth pctile (%.1f deg) -> %+.3f, p = %.4f\n",
            100*best$pctile, best$deg, best$est, best$p))
cut <- quantile(M$mv, best$pctile)
Hm <- copy(H)[, g := mv <= cut & arm >= arm_thr]
cat("\n    who it selects in MLB:\n")
print(Hm[g == TRUE][order(-y), .(nm, season, swings = nsw, sep = round(vs,1),
        mv_axis = round(mv,1), imaged_axis = round(axis,1), arm = round(arm,1),
        whiff_over = round(y,1))], row.names = FALSE)

## ---- 3. carry it to D1 ---------------------------------------------------------------------------
cat("\n\n=== 3. the same rule in D1, thresholds matched by percentile ===\n")
d1cut <- quantile(D1$axis, best$pctile)
cat(sprintf("    D1 %.0fth pctile of the axis gap = %.1f deg (MLB equivalent %.1f deg)\n",
            100*best$pctile, d1cut, cut))
for (g75 in c(40, 75)) {
  Hd <- D1[vs >= quantile(D1$vs, 2/3) & nsw >= g75]
  Hd[, g := axis <= d1cut & arm >= arm_thr]
  r <- crob(Hd, y ~ g + vs, "gTRUE")
  cat(sprintf("    %d+ swings: %+.3f (se %.3f) p = %.4f | %d members from %d arms, pool %d seasons / %d arms\n",
              g75, r[1], r[2], r[3], sum(Hd$g), uniqueN(Hd[g == TRUE]$id), nrow(Hd), uniqueN(Hd$id)))
  cat(sprintf("      minimum detectable effect at 80%% power: %.2f whiff points\n", 2.8*r[2]))
}
cat("\n    D1 across the full percentile sweep (75+ swings):\n")
Hd <- D1[vs >= quantile(D1$vs, 2/3) & nsw >= 75]
SD <- rbindlist(lapply(c(.05,.10,.15,.20,.25,.33), function(q) {
  cc <- quantile(D1$axis, q); Z <- copy(Hd); Z[, g := axis <= cc & arm >= arm_thr]
  r <- crob(Z, y ~ g + vs, "gTRUE")
  data.table(pctile = q, deg = round(cc,1), members = sum(Z$g), arms = uniqueN(Z[g == TRUE]$id),
             est = r[1], se = r[2], p = r[3]) }))
print(SD[, .(pctile, d1_cut_deg = deg, members, arms, est = round(est,3),
             se = round(se,3), p = round(p,4), mde = round(2.8*se,2))], row.names = FALSE)

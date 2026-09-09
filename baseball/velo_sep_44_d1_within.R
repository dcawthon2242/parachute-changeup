#!/usr/bin/env Rscript

# THE TEACHING TEST, RUN WHERE IT HAS POWER
#
# The MLB within-arm test could not settle anything. Season run value repeats at r = .20 and the
# smallest effect the design could detect was about 1.0 runs per 100 per SD, three to five times the
# effect being looked for. Whiff repeats at r = .67 and is the outcome that survives.
#
# D1 is the better laboratory for a "can it be taught" question for three reasons: college arms
# change their pitches far more than major leaguers, the 2023-2025 window covers the kick-change
# adoption wave, and there are several times as many of them. It carries no run value, so this is a
# whiff test only - which is the outcome that has enough signal to test anyway.
#
# One measurement caveat governs the reading. TrackMan's changeup spin axis is close to a
# deterministic function of the break, so the D1 axis gap is a movement-direction gap and not imaged
# rotation. Separation and arm slot survive the translation; the axis component does not. So D1 can
# test the velocity-separation half of the profile cleanly and the axis half only weakly.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 215); MDIR <- "data/statcast_model"

R <- readRDS(file.path(MDIR, "ncaa_d1_pitches.rds")); setDT(R)
R <- R[!is.na(PitcherId) & PitcherId != 0 & is.finite(RelSpeed)]
R[, `:=`(id = as.character(PitcherId), lh = PitcherThrows == "Left")]
R[, pt2 := fifelse(AutoPitchType %in% c("Four-Seam","Fastball"), "FF",
           fifelse(AutoPitchType == "Sinker", "SI",
           fifelse(AutoPitchType %in% c("Changeup","Splitter"), "CH", "other")))]

# fastball anchor: four-seam if he throws enough of them, otherwise sinker, matching the MLB build
FB <- R[pt2 %in% c("FF","SI"), .(n = .N, velo = mean(RelSpeed, na.rm = TRUE),
        ivb = mean(InducedVertBreak, na.rm = TRUE), hb = mean(HorzBreak, na.rm = TRUE),
        relh = mean(RelHeight, na.rm = TRUE), rels = mean(RelSide, na.rm = TRUE),
        ax = mean(ax, na.rm = TRUE), az = mean(az, na.rm = TRUE)), by = .(id, season, pt2)]
setorder(FB, id, season, -n)
FB <- FB[, .SD[fifelse(any(pt2 == "FF" & n >= 100), which(pt2 == "FF" & n >= 100)[1], 1L)],
         by = .(id, season)]
setnames(FB, c("velo","ivb","hb","relh","rels","ax","az","n","pt2"),
         c("f_velo","f_ivb","f_hb","f_relh","f_rels","f_ax","f_az","f_n","f_type"))

CH <- R[pt2 == "CH"]
CS <- CH[, .(np = .N, nm = Pitcher[1], lh = lh[1],
             velo = mean(RelSpeed, na.rm = TRUE), spin = mean(SpinRate, na.rm = TRUE),
             ivb = mean(InducedVertBreak, na.rm = TRUE), hb = mean(HorzBreak, na.rm = TRUE),
             relh = mean(RelHeight, na.rm = TRUE), rels = mean(RelSide, na.rm = TRUE),
             ax = mean(ax, na.rm = TRUE), az = mean(az, na.rm = TRUE),
             nsw = sum(swing == 1, na.rm = TRUE),
             whiff = 100*mean(whiff[swing == 1] == 1, na.rm = TRUE)), by = .(id, season)]
S <- merge(CS, FB[, .(id, season, f_velo, f_ivb, f_hb, f_relh, f_rels, f_ax, f_az, f_n, f_type)],
           by = c("id","season"))
S <- S[f_n >= 100 & is.finite(whiff)]

# mirror lefties so both hands share a frame, then build the same four traits as the MLB panel
S[, `:=`(mx = fifelse(lh, -ax, ax), fmx = fifelse(lh, -f_ax, f_ax))]
S[, `:=`(sep = f_velo - velo,
         slot = atan2(relh - 5.0, abs(rels))*180/pi,
         f_slot = atan2(f_relh - 5.0, abs(f_rels))*180/pi)]
S[, axis := abs(((atan2(mx, az + 32.174)*180/pi - atan2(fmx, f_az + 32.174)*180/pi + 180) %% 360) - 180)]
S[, armgap := abs(slot - f_slot)]
S <- S[is.finite(sep) & is.finite(axis) & is.finite(slot) & np >= 40 & nsw >= 25]

cat(sprintf("D1 panel: %d pitcher-seasons, %d arms, %d changeups.\n", nrow(S), uniqueN(S$id), sum(S$np)))
k <- S[, .N, by = id]
cat(sprintf("  arms with 2+ seasons: %d   with 3: %d\n", k[N >= 2, .N], k[N >= 3, .N]))
cat(sprintf("  median changeups per season: %.0f, median swings: %.0f\n\n",
            median(S$np), median(S$nsw)))

VARS <- c("sep","axis","slot"); GOODDIR <- c(sep = 1, axis = -1, slot = 1)
NICE <- c(sep = "Separation (mph)", axis = "Movement-axis gap (deg)", slot = "Arm slot (deg)")

## =============================================================================================
cat("=== 1. DO D1 ARMS MOVE THESE TRAITS MORE THAN BIG LEAGUERS? ===\n\n")
setorder(S, id, season)
S[, gap_yr := season - shift(season), by = id]
for (v in c(VARS,"whiff","velo","spin","f_velo")) S[, (paste0("d_",v)) := get(v) - shift(get(v)), by = id]
for (v in c("whiff","sep","axis","slot")) S[, (paste0("l_",v)) := shift(get(v)), by = id]
DD <- S[gap_yr == 1 & is.finite(d_sep)]
cat(sprintf("  %d consecutive-season pairs, %d arms.\n\n", nrow(DD), uniqueN(DD$id)))
cat(sprintf("  %-26s %11s %12s %8s %14s\n", "trait", "between SD", "yr-to-yr SD", "ratio",
            "% moving 1+ SD"))
for (v in VARS) {
  b <- S[, sd(get(v))]; w <- DD[, sd(get(paste0("d_",v)))]
  cat(sprintf("  %-26s %11.2f %12.2f %8.2f %13.0f%%\n", NICE[[v]], b, w, w/b,
              100*DD[, mean(abs(get(paste0("d_",v))) > b)]))
}
cat("\n  MLB for reference: separation 0.36, axis 0.32, slot 0.28.\n")

## =============================================================================================
cat("\n=== 2. WITHIN-ARM: DOES MOVING TOWARD THE PROFILE BUY WHIFFS? ===\n\n")
ci <- function(dat, f, k, wt = "nsw") {
  environment(f) <- environment()
  dat <- dat[complete.cases(dat[, c(all.vars(f), "id", wt), with = FALSE])]
  w <- as.numeric(dat[[wt]]); m <- lm(f, dat, weights = w)
  u <- residuals(m)*sqrt(w); X <- model.matrix(m)*sqrt(w); nc <- uniqueN(dat$id)
  keep <- !is.na(coef(m)); X <- X[, keep, drop = FALSE]
  b <- tryCatch(solve(crossprod(X)), error = function(e) MASS::ginv(crossprod(X)))
  V <- b %*% crossprod(rowsum(X*u, dat$id)) %*% b * (nc/(nc-1))
  j <- which(names(coef(m))[keep] == k)
  e <- unname(coef(m)[keep][j]); s <- unname(sqrt(diag(V))[j]); t <- qt(.975, nc-1)
  c(est = e, se = s, lo = e-t*s, hi = e+t*s, p = 2*pt(-abs(e/s), nc-1), n = nrow(dat), arms = nc)
}
sr <- function(r) sprintf("%+7.3f %7.3f  [%+6.2f,%+6.2f] %5d%s", r["est"], r["se"], r["lo"], r["hi"],
                          r["arms"], ifelse(r["p"]<.01," **",ifelse(r["p"]<.05," * ","   ")))
S2 <- copy(S)
for (v in VARS) S2[, (paste0("z_",v)) := scale(get(v))[,1]*GOODDIR[[v]]]
S2[, toward := rowMeans(cbind(z_sep, z_axis))]
cat(sprintf("  %-40s %7s %7s %16s %5s\n", "specification", "est", "se", "95% interval", "arms"))
for (v in c(VARS, "toward")) {
  k <- if (v == "toward") "toward" else paste0("z_", v)
  lab <- if (v == "toward") "COMPOSITE (sep + axis)" else NICE[[v]]
  cat(sprintf("  %-40s %s\n", paste("between arms:", lab),
              sr(ci(S2, as.formula(paste0("whiff ~ ", k)), k))))
}
cat("\n")
for (v in c(VARS, "toward")) {
  k <- if (v == "toward") "toward" else paste0("z_", v)
  lab <- if (v == "toward") "COMPOSITE (sep + axis)" else NICE[[v]]
  f <- as.formula(paste0("whiff ~ ", k, " + factor(id) + factor(season)"))
  cat(sprintf("  %-40s %s\n", paste("within arms:", lab), sr(ci(S2, f, k))))
}

cat("\n  first differences, with last year's whiff controlled for reversion:\n\n")
Z <- copy(DD)
for (v in VARS) Z[, (paste0("z_",v)) := scale(get(paste0("d_",v)))[,1]*GOODDIR[[v]]]
Z[, toward := rowMeans(cbind(z_sep, z_axis))]
for (v in c(VARS, "toward")) {
  k <- if (v == "toward") "toward" else paste0("z_", v)
  lab <- if (v == "toward") "COMPOSITE (sep + axis)" else NICE[[v]]
  cat(sprintf("  %-40s %s\n", lab, sr(ci(Z, as.formula(paste0("d_whiff ~ ", k, " + l_whiff + factor(season)")), k))))
}

cat("\n  dose-response by size of move:\n\n")
Z[, tq := cut(toward, quantile(toward, 0:5/5), include.lowest = TRUE,
              labels = c("away, hard","away","flat","toward","toward, hard"))]
print(Z[, .(pairs = .N, `mean move` = round(mean(toward),2), `d sep` = round(mean(d_sep),2),
            `d axis` = round(mean(d_axis),1), `d whiff` = round(weighted.mean(d_whiff, nsw),2),
            `whiff before` = round(weighted.mean(l_whiff, nsw),1)), by = tq][order(tq)], row.names = FALSE)

## =============================================================================================
cat("\n=== 3. THE CLEANEST VERSION: ARMS WHO ADDED REAL SEPARATION ===\n\n")
# separation is the trait D1 measures without ambiguity and the one a pitching coach can actually
# target. a half-mph change is noise; a 1.5 mph change in one offseason is an intervention.
DD[, grp := fifelse(d_sep >= 1.5, "added 1.5+ mph",
            fifelse(d_sep <= -1.5, "lost 1.5+ mph", "held steady"))]
G <- DD[, .(pairs = .N, arms = uniqueN(id), `d sep` = round(mean(d_sep),2),
            `whiff before` = round(weighted.mean(l_whiff, nsw),2),
            `whiff after` = round(weighted.mean(whiff, nsw),2),
            change = round(weighted.mean(d_whiff, nsw),2)), by = grp]
print(G[order(-pairs)], row.names = FALSE)
for (g in c("added 1.5+ mph","lost 1.5+ mph")) {
  X <- DD[grp %in% c(g, "held steady")][, hit := grp == g]
  r <- ci(X, d_whiff ~ hit + l_whiff, "hitTRUE")
  cat(sprintf("\n  %-16s vs held steady: %+.2f whiff points  [%+.2f, %+.2f]  p = %.4f\n",
              g, r["est"], r["lo"], r["hi"], r["p"]))
}
cat(sprintf("\n  minimum detectable effect at 80%% power, within-arm composite: %.2f whiff points per SD\n",
            2.8*ci(S2, whiff ~ toward + factor(id) + factor(season), "toward")["se"]))

## =============================================================================================
cat("\n=== 4. HOW DID THEY GET THE SEPARATION? ===\n\n")
# this decides whether the finding is coachable. separation can widen because the changeup slowed
# down, which a coach can instruct, or because the fastball sped up, which is a different project
# entirely and not a changeup intervention at all.
setnames(DD, "d_f_velo", "d_fvelo")
MECH <- DD[is.finite(d_fvelo), .(pairs = .N,
   `d separation` = round(mean(d_sep),2),
   `from slower CH` = round(mean(-d_velo),2),
   `from faster FB` = round(mean(d_fvelo),2),
   `CH share of move` = paste0(round(100*mean(-d_velo)/mean(d_sep)), "%"),
   `d CH spin` = round(mean(d_spin, na.rm = TRUE)),
   `d whiff` = round(weighted.mean(d_whiff, nsw),2)), by = grp]
print(MECH[order(-pairs)], row.names = FALSE)
cat("\n  among the arms who added 1.5+ mph, how each one did it:\n\n")
A <- DD[grp == "added 1.5+ mph" & is.finite(d_fvelo)]
cat(sprintf("    changeup slowed by 1+ mph, fastball flat   %2d of %d\n",
            A[-d_velo >= 1 & abs(d_fvelo) < 1, .N], nrow(A)))
cat(sprintf("    fastball gained 1+ mph, changeup flat      %2d of %d\n",
            A[d_fvelo >= 1 & abs(d_velo) < 1, .N], nrow(A)))
cat(sprintf("    both moved                                 %2d of %d\n",
            A[-d_velo >= 1 & d_fvelo >= 1, .N], nrow(A)))
cat(sprintf("    neither cleanly (small mixed moves)        %2d of %d\n",
            A[!((-d_velo >= 1 & abs(d_fvelo) < 1) | (d_fvelo >= 1 & abs(d_velo) < 1) |
                (-d_velo >= 1 & d_fvelo >= 1)), .N], nrow(A)))
cat("\n  splitting the effect by route:\n\n")
A2 <- DD[grp %in% c("added 1.5+ mph","held steady") & is.finite(d_fvelo)]
A2[, route := fifelse(grp == "held steady", "held steady",
             fifelse(-d_velo >= d_fvelo, "added via slower changeup", "added via faster fastball"))]
for (r in c("added via slower changeup","added via faster fastball")) {
  X <- A2[route %in% c(r, "held steady")][, hit := route == r]
  z <- ci(X, d_whiff ~ hit + l_whiff, "hitTRUE")
  cat(sprintf("  %-30s %+.2f whiff points  [%+.2f, %+.2f]  p = %.4f  (n = %d)\n",
              r, z["est"], z["lo"], z["hi"], z["p"], X[hit == TRUE, .N]))
}

## =============================================================================================
cat("\n=== 5. THE SAME TEST IN MLB, AND WHY IT COMES BACK EMPTY ===\n\n")
TT <- readRDS(file.path(MDIR, "teachable.rds"))$pairs; setDT(TT)
TT[, grp := fifelse(d_sep >= 1.5, "added 1.5+ mph",
           fifelse(d_sep <= -1.5, "lost 1.5+ mph", "held steady"))]
cat(sprintf("  D1 arms making a 1.5+ mph separation change: %d of %d pairs (%.0f%%)\n",
            DD[grp != "held steady", .N], nrow(DD), 100*DD[, mean(grp != "held steady")]))
cat(sprintf("  MLB arms making the same change:             %d of %d pairs (%.0f%%)\n\n",
            TT[grp != "held steady", .N], nrow(TT), 100*TT[, mean(grp != "held steady")]))
print(TT[, .(pairs = .N, `d sep` = round(mean(d_sep),2),
             `whiff before` = round(weighted.mean(l_whiff, nsw2),2),
             `whiff after` = round(weighted.mean(whiff, nsw2),2),
             change = round(weighted.mean(d_whiff, nsw2),2)), by = grp][order(-pairs)], row.names = FALSE)
X <- TT[grp %in% c("added 1.5+ mph","held steady")][, hit := grp == "added 1.5+ mph"]
if (X[hit == TRUE, .N] >= 5) {
  z <- ci(X, d_whiff ~ hit + l_whiff, "hitTRUE", "nsw2")
  cat(sprintf("\n  MLB, added 1.5+ mph vs held steady: %+.2f whiff points  [%+.2f, %+.2f]  p = %.4f  (n = %d)\n",
              z["est"], z["lo"], z["hi"], z["p"], X[hit == TRUE, .N]))
  cat(sprintf("  smallest effect that MLB sample could detect: %.2f whiff points\n", 2.8*z["se"]))
}
saveRDS(list(panel = S, pairs = DD), file.path(MDIR, "d1_within.rds"))

#!/usr/bin/env Rscript

# THE SWING MINIMUM, APPLIED HONESTLY.
#
# Raising the swing floor is a legitimate pre-registration choice: a season-level whiff rate built
# on 50 swings is a noisy estimate of anything, and dropping those should tighten every standard
# error in the study. But it only counts as that if it is applied to the whole population. Raising
# the floor for bin members while leaving the comparison group at 40 is not a precision gain, it is
# a selection rule wearing one.
#
# So: the floor is re-applied to the full MLB universe, the separation tercile is recomputed inside
# the surviving population, and every bin is rebuilt from scratch at each floor. The sweep runs
# 40 / 50 / 60 / 75 / 100 / 125 so a reader can see whether 75 sits on a smooth trend or on a cliff.
#
# One thing to watch, stated up front rather than discovered later: at a floor of 75 the seasons
# that drop out of the original bin are Roark 2020 (62 swings), Skubal 2020 (53) and Cease 2020
# (52) - and Cease 2020 is the worst-performing member in the set at -9.6. Any improvement is
# therefore part precision and part removal of the biggest negative. Section 3 separates the two.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
set.seed(53); options(width = 205); MDIR <- "data/statcast_model"

L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L); M <- L[league == "MLB"]
NM <- unique(readRDS(file.path(MDIR, "parachute_ff.rds"))[, .(pitcher, player_name)]); setDT(NM)
NM <- unique(NM, by = "pitcher")[, .(id = as.character(pitcher), nm = sub(",.*", "", player_name))]
M <- merge(M, NM, by = "id", all.x = TRUE)
P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(ax) & is.finite(az) & is.finite(ax_diff) & is.finite(az_diff)]
P[, `:=`(mx = fifelse(p_throws == "L", -ax, ax), mz = az + 32.174,
         fx = fifelse(p_throws == "L", -(ax - ax_diff), ax - ax_diff), fz = (az - az_diff) + 32.174)]
A <- P[, .(n = .N, ca = atan2(mean(mx), mean(mz))*180/pi,
           fa = atan2(mean(fx), mean(fz))*180/pi), by = .(pitcher, season)][n >= 40]
A[, `:=`(id = as.character(pitcher), mv = abs(((ca - fa + 180) %% 360) - 180))]
M <- merge(M, A[, .(id, season, mv)], by = c("id","season")); M[, disc := mv - axis]
CH <- readRDS(file.path(MDIR, "mlb_chase_resid.rds")); setDT(CH)
CH <- CH[, .(ych = 100*mean(r_all, na.rm = TRUE)), by = .(pitcher, season)][, id := as.character(pitcher)]
VR <- readRDS(file.path(MDIR, "velo_sep_resid.rds"))
RV <- as.data.table(VR$R)[, .(yrv = 100*mean(rv - q_aware)), by = .(pitcher, season)][, id := as.character(pitcher)]
WB <- as.data.table(VR$W)[, .(blind = 100*mean(whiff - p_blind)), by = .(pitcher, season)][, id := as.character(pitcher)]
M <- merge(M, CH[, .(id, season, ych)], by = c("id","season"), all.x = TRUE)
M <- merge(M, RV[, .(id, season, yrv)], by = c("id","season"), all.x = TRUE)
M <- merge(M, WB[, .(id, season, blind)], by = c("id","season"), all.x = TRUE)

crob <- function(D, f, k, w = FALSE) {
  environment(f) <- environment()
  D <- D[complete.cases(D[, c(all.vars(f), "id", "nsw"), with = FALSE])]
  if (uniqueN(D$id) < 4) return(c(NA, NA, NA))
  wt <- if (w) as.numeric(D$nsw) else rep(1, nrow(D))
  m <- lm(f, D, weights = wt); u <- residuals(m)*sqrt(wt); X <- model.matrix(m)*sqrt(wt)
  nc <- uniqueN(D$id); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k]); c(e, s, 2*pt(-abs(e/s), nc-1))
}
build <- function(floor) {
  U <- M[nsw >= floor]
  vc <- quantile(U$vs, 2/3); H <- U[vs >= vc]
  H[, `:=`(orig = axis < 10 & arm >= arm_thr,
           lock = axis < 10 & arm >= arm_thr & disc <= 19.71,
           wide = axis < 11 & arm >= quantile(U$arm, .55, na.rm = TRUE) & disc <= 27.5)]
  list(U = U, H = H, vc = vc)
}

## =============================================================================================
cat("=== 1. SWEEP OF THE SWING FLOOR, APPLIED TO THE WHOLE POPULATION ===\n\n")
S <- rbindlist(lapply(c(40, 50, 60, 75, 100, 125), function(fl) {
  B <- build(fl); H <- B$H
  rbindlist(lapply(c("orig","lock","wide"), function(b) {
    H[, gg := H[[b]]]
    r <- crob(H, y ~ gg + vs, "ggTRUE"); rw <- crob(H, y ~ gg + vs, "ggTRUE", TRUE)
    rc <- crob(H, y ~ gg + vs + axis + arm, "ggTRUE", TRUE)
    data.table(floor = fl, spec = b, pop = nrow(B$U), hi = nrow(H), n = sum(H$gg),
               arms = uniqueN(H[gg == TRUE]$id), est = r[1], se = r[2], p = r[3],
               estw = rw[1], pw = rw[3], ctrl = rc[1], pctrl = rc[3])
  }))
}))
for (b in c("orig","lock","wide")) {
  lab <- c(orig = "ORIGINAL  axis<10 + high slot", lock = "LOCKED    + discrepancy <= 19.71",
           wide = "WIDE      axis<11, 55th pct, disc <= 27.5")[b]
  cat(sprintf("--- %s ---\n", lab))
  print(S[spec == b, .(floor, population = pop, high_sep = hi, bin_n = n, arms,
                       unweighted = round(est,2), se = round(se,2), p = round(p,4),
                       weighted = round(estw,2), p_w = round(pw,4),
                       with_controls = round(ctrl,2), p_ctrl = round(pctrl,4))], row.names = FALSE)
  cat("\n")
}

## =============================================================================================
cat("=== 2. THE ORIGINAL SET AT A 75-SWING FLOOR ===\n\n")
B <- build(75); H <- B$H
cat(sprintf("population %d seasons (from %d at floor 40) | separation cut %.2f mph | high-sep third %d\n",
            nrow(B$U), M[, .N], B$vc, nrow(H)))
for (b in c("orig","lock","wide")) {
  H[, gg := H[[b]]]
  cat(sprintf("\n%s: %d seasons, %d arms, %d swings\n", toupper(b), sum(H$gg),
              uniqueN(H[gg == TRUE]$id), H[gg == TRUE, sum(nsw)]))
  for (o in list(c("y","whiff"), c("ych","chase"), c("yrv","run value"), c("blind","whiff vs blind"))) {
    r <- crob(H[is.finite(get(o[1]))], as.formula(paste(o[1], "~ gg + vs")), "ggTRUE", TRUE)
    cat(sprintf("   %-16s %+7.3f (se %5.3f) p = %.4f\n", o[2], r[1], r[2], r[3]))
  }
  r <- crob(H, y ~ gg + vs + axis + arm, "ggTRUE", TRUE)
  cat(sprintf("   %-16s %+7.3f (se %5.3f) p = %.4f\n", "whiff + controls", r[1], r[2], r[3]))
}
H[, gg := orig]
cat("\nthe original bin at floor 75:\n")
print(H[orig == TRUE][order(-y)][, .(pitcher = nm, season, swings = nsw, sep = round(vs,1),
      axis = round(axis,1), arm = round(arm,1), disc = round(disc,1), whiff = round(y,1),
      chase = round(ych,1), rv = round(yrv,2), blind = round(blind,1))], row.names = FALSE)
D40 <- build(40)$H
cat("\ndropped by the floor (were in the original bin at 40 swings):\n")
print(D40[orig == TRUE & nsw < 75][order(-y)][, .(pitcher = nm, season, swings = nsw,
      whiff = round(y,1), rv = round(yrv,2))], row.names = FALSE)

## =============================================================================================
cat("\n=== 3. PRECISION OR SELECTION? ===\n")
# A floor that only buys precision shrinks the standard error and leaves the point estimate alone.
# A floor that is really a selection rule moves the point estimate. Decompose by re-estimating the
# floor-40 bin on the floor-40 data but weighting each season the way the floor-75 sample would.
o40 <- crob(D40, y ~ orig + vs, "origTRUE"); o75 <- crob(H, y ~ orig + vs, "origTRUE")
cat(sprintf("  floor 40: %+.3f (se %.3f) p %.4f on %d seasons\n", o40[1], o40[2], o40[3], sum(D40$orig)))
cat(sprintf("  floor 75: %+.3f (se %.3f) p %.4f on %d seasons\n", o75[1], o75[2], o75[3], sum(H$orig)))
cat(sprintf("  point estimate moved %+.3f, standard error moved %+.3f\n", o75[1]-o40[1], o75[2]-o40[2]))
# how much of the estimate move is just deleting Cease 2020?
noc <- crob(D40[!(nm == "Cease" & season == 2020)], y ~ orig + vs, "origTRUE")
cat(sprintf("  floor 40 with only Cease 2020 deleted: %+.3f (se %.3f) p %.4f\n", noc[1], noc[2], noc[3]))
cat(sprintf("     -> of the %+.3f estimate gain, %+.3f comes from that one season alone (%.0f%%)\n",
            o75[1]-o40[1], noc[1]-o40[1], 100*(noc[1]-o40[1])/(o75[1]-o40[1])))
# does raising the floor help random bins too? if yes it is a general precision effect
cat("\n  does the floor help RANDOM bins of the same size? 1500 draws at each floor:\n")
for (fl in c(40, 75)) {
  BB <- build(fl); HH <- BB$H; k <- sum(HH$orig)
  d <- replicate(1500, { HH[, gp := FALSE]; HH[sample(.N, k), gp := TRUE]
    crob(HH, y ~ gp + vs, "gpTRUE")[c(1,2)] })
  cat(sprintf("    floor %3d: random-bin |estimate| mean %.3f, mean se %.3f\n",
              fl, mean(abs(d[1,])), mean(d[2,])))
}

## =============================================================================================
cat("\n=== 4. ROBUSTNESS OF THE ORIGINAL BIN AT FLOOR 75 ===\n")
H[, gg := orig]; obs <- crob(H, y ~ gg + vs, "ggTRUE")[1]
k <- sum(H$gg); ka <- uniqueN(H[gg == TRUE]$id); arms <- unique(H$id)
p1 <- mean(replicate(4000, { H[, gp := FALSE]; H[sample(.N, k), gp := TRUE]
  crob(H, y ~ gp + vs, "gpTRUE")[1] }) >= obs)
p2 <- mean(replicate(2000, { s <- sample(arms, ka); H[, gp := id %in% s]
  crob(H, y ~ gp + vs, "gpTRUE")[1] }) >= obs)
cat(sprintf("  permutation: season-shuffle p = %.4f | arm-shuffle p = %.4f\n", p1, p2))
ar <- unique(H[gg == TRUE]$id)
L1 <- rbindlist(lapply(ar, function(a) { x <- crob(H[id != a], y ~ gg + vs, "ggTRUE")
  data.table(drop = H[id == a, unique(nm)], n = sum(H[id != a]$gg), est = x[1], p = x[3]) }))
setorder(L1, -p); print(L1[, .(drop, n, est = round(est,2), p = round(p,4))], row.names = FALSE)
L2 <- rbindlist(apply(combn(length(ar), 2), 2, function(ix) { a <- ar[ix]
  x <- crob(H[!id %in% a], y ~ gg + vs, "ggTRUE")
  data.table(d = paste(sort(H[id %in% a, unique(nm)]), collapse = " + "), est = x[1], p = x[3]) }))
cat(sprintf("  LOO %d/%d hold p<.05 | L2O %d/%d hold, worst %+.2f (p %.4f, %s)\n",
            sum(L1$p < .05), nrow(L1), sum(L2$p < .05), nrow(L2),
            L2[which.max(p)]$est, max(L2$p), L2[which.max(p)]$d))
d <- crob(H[!nm %in% c("Skubal","Cease")], y ~ gg + vs, "ggTRUE")
cat(sprintf("  Skubal and Cease both removed: %+.3f (se %.3f) p = %.4f on %d seasons\n",
            d[1], d[2], d[3], H[!nm %in% c("Skubal","Cease") & orig == TRUE, .N]))
LO <- M[nsw >= 75][vs < B$vc]
LO[, gl := axis < 10 & arm >= arm_thr]
r <- crob(LO, y ~ gl + vs, "glTRUE")
cat(sprintf("  low-separation placebo: %+.3f (se %.3f) p = %.4f on %d seasons\n",
            r[1], r[2], r[3], sum(LO$gl)))
saveRDS(list(sweep = S, H75 = H), file.path(MDIR, "swing75.rds"))

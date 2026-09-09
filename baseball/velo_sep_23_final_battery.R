#!/usr/bin/env Rscript

# EVERYTHING, RUN AGAINST THE 10-MEMBER SET.
#
# The definition is now three conditions, all fixed before this script:
#   (1) velocity separation in the top third of changeups
#   (2) imaged spin-axis gap to the four-seam under 10 degrees
#   (3) arm slot at or above the population threshold
#   (4) movement-axis gap minus imaged-axis gap at or under 19.71 degrees
# The fourth is the new one, and its cut came from a blind largest-gap split rather than from
# naming who had to fall out.
#
# The point of this pass is to find out which claims survive contact with everything, not to
# assemble a list of the ones that look good. Sections are ordered so the ones most likely to kill
# the finding run first: multiplicity, then permutation, then the outcomes that were never the
# target. A result that only appears in the outcome the definition was tuned on is not a result.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
set.seed(41); options(width = 200); MDIR <- "data/statcast_model"

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
M <- merge(M, A[, .(id, season, mv)], by = c("id","season"))
M[, disc := mv - axis]

VSCUT <- quantile(M$vs, 2/3); DCUT <- 19.71
H <- M[vs >= VSCUT]
H[, `:=`(g13 = axis < 10 & arm >= arm_thr, g = axis < 10 & arm >= arm_thr & disc <= DCUT)]

crob <- function(D, f, k, w = FALSE) {
  environment(f) <- environment()
  D <- D[complete.cases(D[, c(all.vars(f), "id", "nsw"), with = FALSE])]
  wt <- if (w) as.numeric(D$nsw) else rep(1, nrow(D))
  m <- lm(f, D, weights = wt); u <- residuals(m)*sqrt(wt); X <- model.matrix(m)*sqrt(wt)
  nc <- uniqueN(D$id); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k]); c(e, s, 2*pt(-abs(e/s), nc-1))
}
line <- function(lab, r, extra = "") cat(sprintf("  %-42s %+7.3f (se %5.3f)  p = %.4f  %s\n",
                                                 lab, r[1], r[2], r[3], extra))

cat(sprintf("MLB universe %d pitcher-seasons | high-separation third: %d (vs >= %.2f mph)\n",
            nrow(M), nrow(H), VSCUT))
cat(sprintf("gated set: %d seasons, %d arms | prior 13-season set: %d seasons, %d arms\n\n",
            sum(H$g), uniqueN(H[g == TRUE]$id), sum(H$g13), uniqueN(H[g13 == TRUE]$id)))

## =============================================================================================
cat("=== 1. MULTIPLICITY: where does this cell rank in the space it was found in? ===\n")
# The three thresholds were each chosen with the data in view. The honest question is how many
# other threshold triples would have produced an effect this large. Grid over plausible values of
# all three, score every cell that retains at least 6 seasons, and rank the chosen one.
grid <- CJ(ax = seq(6, 20, by = 1), ap = seq(0.50, 0.90, by = 0.05), dc = seq(10, 45, by = 2.5))
sc <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  th <- quantile(M$arm, grid$ap[i], na.rm = TRUE)
  H[, gg := axis < grid$ax[i] & arm >= th & disc <= grid$dc[i]]
  if (sum(H$gg) < 6 || uniqueN(H[gg == TRUE]$id) < 4) return(NULL)
  r <- try(crob(H, y ~ gg + vs, "ggTRUE"), silent = TRUE); if (inherits(r, "try-error")) return(NULL)
  data.table(ax = grid$ax[i], ap = grid$ap[i], dc = grid$dc[i], n = sum(H$gg), est = r[1], p = r[3])
}))
setorder(sc, -est)
sc[, rank := .I]
own <- sc[ax == 10 & abs(ap - 0.70) < 1e-9 & dc == 20]
cat(sprintf("  %d valid threshold triples scored. %d reach p<.05 (%.0f%%), %d reach est > +4.0\n",
            nrow(sc), sum(sc$p < .05), 100*mean(sc$p < .05), sum(sc$est > 4)))
cat(sprintf("  median effect across the whole grid: %+.3f | 90th pct %+.3f | max %+.3f\n",
            median(sc$est), quantile(sc$est, .9), max(sc$est)))
if (nrow(own)) cat(sprintf("  nearest grid neighbour of the chosen cell ranks %d of %d (est %+.3f, p %.4f)\n",
                           own$rank[1], nrow(sc), own$est[1], own$p[1]))
cat("  top 6 cells by effect:\n")
print(head(sc[, .(axis_lt = ax, arm_pct = ap, disc_le = dc, n, est = round(est,2), p = round(p,4))], 6),
      row.names = FALSE)

## =============================================================================================
cat("\n=== 2. PERMUTATION: how often does a random 10-season set beat it? ===\n")
# Membership is shuffled within the high-separation third, holding set size fixed, so the null
# keeps the outcome distribution and the separation control intact and only breaks the link to
# who is actually in the bin.
obs <- crob(H, y ~ g + vs, "gTRUE")[1]
nperm <- 4000; idx <- which(H$g); k <- length(idx)
pm <- replicate(nperm, {
  H[, gp := FALSE]; H[sample(.N, k), gp := TRUE]
  suppressWarnings(try(crob(H, y ~ gp + vs, "gpTRUE")[1], silent = TRUE)) })
pm <- as.numeric(pm[!is.na(suppressWarnings(as.numeric(pm)))])
cat(sprintf("  observed %+.3f | permutation mean %+.3f sd %.3f | p_perm = %.4f (%d of %d exceed)\n",
            obs, mean(pm), sd(pm), mean(pm >= obs), sum(pm >= obs), length(pm)))
# and a second null that respects arm clustering: shuffle at the pitcher level
arms <- unique(H$id); ka <- uniqueN(H[g == TRUE]$id)
pa <- replicate(2000, {
  s <- sample(arms, ka); H[, gp := id %in% s]
  suppressWarnings(try(crob(H, y ~ gp + vs, "gpTRUE")[1], silent = TRUE)) })
pa <- as.numeric(pa[!is.na(suppressWarnings(as.numeric(pa)))])
cat(sprintf("  arm-level null (shuffle %d whole arms):  p_perm = %.4f\n", ka, mean(pa >= obs)))

## =============================================================================================
cat("\n=== 3. THE EFFECT, under every control I can think of ===\n")
for (w in c(FALSE, TRUE)) {
  tag <- if (w) "[swing-weighted] " else "[unweighted]    "
  line(paste0(tag, "separation only"),          crob(H, y ~ g + vs, "gTRUE", w))
  line(paste0(tag, "+ axis, arm continuous"),   crob(H, y ~ g + vs + axis + arm, "gTRUE", w))
  line(paste0(tag, "+ both efficiencies"),      crob(H, y ~ g + vs + axis + arm + ec + ef, "gTRUE", w))
  line(paste0(tag, "+ discrepancy continuous"), crob(H, y ~ g + vs + axis + arm + ec + ef + disc, "gTRUE", w))
  line(paste0(tag, "+ season fixed effects"),   crob(H, y ~ g + vs + axis + arm + factor(season), "gTRUE", w))
  line(paste0(tag, "+ sample size"),            crob(H, y ~ g + vs + axis + arm + log(nsw), "gTRUE", w))
  cat("\n")
}
cat("  reference points:\n")
line("  prior 13-season set, no discrepancy gate", crob(H, y ~ g13 + vs, "g13TRUE"))
line("  the 3 excluded seasons vs the other 10",   crob(H[g13 == TRUE][, x := !g], y ~ x, "xTRUE"))
line("  gate applied to the LOW two-thirds of sep",
     crob(M[vs < VSCUT][, gl := axis < 10 & arm >= arm_thr & disc <= DCUT], y ~ gl + vs, "glTRUE"),
     sprintf("n=%d", M[vs < VSCUT & axis < 10 & arm >= arm_thr & disc <= DCUT, .N]))

## =============================================================================================
cat("\n=== 4. OUTCOMES THE DEFINITION WAS NEVER TUNED ON ===\n")
CH <- readRDS(file.path(MDIR, "mlb_chase_resid.rds")); setDT(CH)
CH <- CH[, .(nch = .N, ych = 100*mean(r_all, na.rm = TRUE)), by = .(pitcher, season)][nch >= 40]
CH[, id := as.character(pitcher)]
VR <- readRDS(file.path(MDIR, "velo_sep_resid.rds"))
RV <- as.data.table(VR$R)[, .(nrv = .N, yrv = 100*mean(rv - q_aware)), by = .(pitcher, season)]
RV[, id := as.character(pitcher)]
WB <- as.data.table(VR$W)[, .(nwb = .N, blind = 100*mean(whiff - p_blind),
                              aware = 100*mean(whiff - p_aware)), by = .(pitcher, season)]
WB[, id := as.character(pitcher)]
H2 <- merge(H, CH[, .(id, season, nch, ych)], by = c("id","season"), all.x = TRUE)
H2 <- merge(H2, RV[, .(id, season, nrv, yrv)], by = c("id","season"), all.x = TRUE)
H2 <- merge(H2, WB[, .(id, season, nwb, blind, aware)], by = c("id","season"), all.x = TRUE)
line("chase residual (pts)",              crob(H2[is.finite(ych)], ych ~ g + vs, "gTRUE"))
line("run value residual (runs/100)",     crob(H2[is.finite(yrv)], yrv ~ g + vs, "gTRUE"))
line("whiff vs fastball-BLIND model",     crob(H2[is.finite(blind)], blind ~ g + vs, "gTRUE"))
line("whiff vs separation-AWARE model",   crob(H2[is.finite(aware)], aware ~ g + vs, "gTRUE"))
cat("\n  (run value is signed so negative = better for the pitcher)\n")

## =============================================================================================
cat("\n=== 5. ROBUSTNESS TO ANY ONE ARM ===\n")
ga <- unique(H[g == TRUE]$id)
L1 <- rbindlist(lapply(ga, function(a) { x <- crob(H[id != a], y ~ g + vs, "gTRUE")
  data.table(drop = H[id == a, unique(nm)], n = sum(H[id != a]$g), est = x[1], p = x[3]) }))
setorder(L1, -p); print(L1[, .(drop, n, est = round(est,2), p = round(p,4))], row.names = FALSE)
L2 <- rbindlist(apply(combn(length(ga), 2), 2, function(ix) { a <- ga[ix]
  x <- crob(H[!id %in% a], y ~ g + vs, "gTRUE")
  data.table(d = paste(sort(H[id %in% a, unique(nm)]), collapse = " + "), est = x[1], p = x[3]) }))
cat(sprintf("\n  leave-one-arm-out:  %d of %d keep p<.05, weakest %+.3f (p %.4f)\n",
            sum(L1$p < .05), nrow(L1), L1$est[1], L1$p[1]))
cat(sprintf("  leave-two-arms-out: %d of %d keep p<.05, weakest %+.3f (p %.4f, %s)\n",
            sum(L2$p < .05), nrow(L2), L2[which.max(p)]$est, max(L2$p), L2[which.max(p)]$d))
d <- crob(H[!nm %in% c("Skubal","Cease")], y ~ g + vs, "gTRUE")
cat(sprintf("  Skubal AND Cease both removed (%d seasons left in bin): %+.3f p = %.4f\n",
            H[!nm %in% c("Skubal","Cease") & g == TRUE, .N], d[1], d[3]))

## =============================================================================================
cat("\n=== 6. WHAT IT IS WORTH ===\n")
sw <- H[g == TRUE, sum(nsw)]; est <- crob(H, y ~ g + vs, "gTRUE", TRUE)[1]
cat(sprintf("  effect %+.2f whiff points per swing, over %d swings in %d seasons\n", est, sw, sum(H$g)))
cat(sprintf("  a full season at the bin median (%d swings) is %+.1f extra whiffs above model\n",
            as.integer(median(H[g == TRUE]$nsw)), est/100*median(H[g == TRUE]$nsw)))
q <- M[, .(n = .N, pct = 100*mean(vs >= VSCUT & axis < 10 & arm >= arm_thr & disc <= DCUT))]
cat(sprintf("  prevalence: %d of %d qualifying MLB pitcher-seasons, %.1f%%\n",
            round(q$n*q$pct/100), q$n, q$pct))

## =============================================================================================
cat("\n=== 7. THE SET ===\n")
S <- H[g13 == TRUE][order(-y)]
S <- merge(S, H2[, .(id, season, ych, yrv, blind)], by = c("id","season"), all.x = TRUE)[order(-y)]
print(S[, .(pitcher = nm, season, swings = nsw, sep = round(vs,1), axis = round(axis,1),
            arm = round(arm,0), disc = round(disc,1), whiff_over = round(y,1),
            chase_over = round(ych,1), rv_over = round(yrv,2), blind_over = round(blind,1),
            IN = g)], row.names = FALSE)
C <- S[g == TRUE, .(seasons = .N, swings = sum(nsw), sep = round(weighted.mean(vs, nsw),1),
                    whiff_over = round(weighted.mean(y, nsw),2)), by = nm][order(-whiff_over)]
cat("\n  by arm:\n"); print(C, row.names = FALSE)
saveRDS(list(H = H2, set = S, grid = sc, perm = pm), file.path(MDIR, "final_gate_battery.rds"))

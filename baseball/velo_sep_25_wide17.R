#!/usr/bin/env Rscript

# The deduped grid surfaced a wider specification that looks better than the locked one on the
# test the locked one failed. Before recommending it, it gets the same battery the 10-season set
# got, plus the comparison run side by side so the choice is on visible evidence.
#
#   locked : axis < 10 | arm >= 70th pct | discrepancy <= 19.71   -> 10 seasons,  7 arms
#   wide17 : axis < 11 | arm >= 55th pct | discrepancy <= 27.5    -> 17 seasons, 10 arms
#
# The wide set is not a loosening of the same idea in one direction - it trades a stricter
# discrepancy cut for a much lower arm-slot bar. Whether that is a real improvement or just a
# larger sample buying significance is exactly what the permutation and leave-out tests decide.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
set.seed(47); options(width = 200); MDIR <- "data/statcast_model"
B <- readRDS(file.path(MDIR, "final_gate_battery.rds")); H <- B$H; setDT(H)
M <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(M); M <- M[league == "MLB"]
VSCUT <- quantile(M$vs, 2/3)

crob <- function(D, f, k, w = FALSE) {
  environment(f) <- environment()
  D <- D[complete.cases(D[, c(all.vars(f), "id", "nsw"), with = FALSE])]
  wt <- if (w) as.numeric(D$nsw) else rep(1, nrow(D))
  m <- lm(f, D, weights = wt); u <- residuals(m)*sqrt(wt); X <- model.matrix(m)*sqrt(wt)
  nc <- uniqueN(D$id); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k]); c(e, s, 2*pt(-abs(e/s), nc-1))
}
th55 <- quantile(M$arm, .55, na.rm = TRUE)
H[, w17 := axis < 11 & arm >= th55 & disc <= 27.5]
cat(sprintf("locked: %d seasons %d arms %d swings | wide17: %d seasons %d arms %d swings\n\n",
            sum(H$g), uniqueN(H[g == TRUE]$id), H[g == TRUE, sum(nsw)],
            sum(H$w17), uniqueN(H[w17 == TRUE]$id), H[w17 == TRUE, sum(nsw)]))

sp <- function(r) sprintf("%+7.3f (se %5.3f) p%.4f", r[1], r[2], r[3])
cat(sprintf("%-42s %28s   %28s\n", "", "LOCKED (10 seasons)", "WIDE17 (17 seasons)"))
row <- function(lab, f, kl, kw, D = H, w = FALSE)
  cat(sprintf("%-42s %s   %s\n", lab, sp(crob(D, as.formula(sub("BIN","g",f)), kl, w)),
              sp(crob(D, as.formula(sub("BIN","w17",f)), kw, w))))

cat("\n--- whiff, adding controls cumulatively (swing-weighted) ---\n")
row("separation only",                "y ~ BIN + vs", "gTRUE", "w17TRUE", H, TRUE)
row("+ axis, arm continuous",         "y ~ BIN + vs + axis + arm", "gTRUE", "w17TRUE", H, TRUE)
row("+ discrepancy continuous",       "y ~ BIN + vs + axis + arm + disc", "gTRUE", "w17TRUE", H, TRUE)
row("+ efficiencies, season FE",      "y ~ BIN + vs + axis + arm + disc + ec + ef + factor(season)",
    "gTRUE", "w17TRUE", H, TRUE)
cat("\n--- and unweighted, the harder test ---\n")
row("separation only",                "y ~ BIN + vs", "gTRUE", "w17TRUE")
row("+ axis, arm continuous",         "y ~ BIN + vs + axis + arm", "gTRUE", "w17TRUE")
row("+ discrepancy continuous",       "y ~ BIN + vs + axis + arm + disc", "gTRUE", "w17TRUE")

cat("\n--- outcomes never tuned on ---\n")
row("chase above model",              "ych ~ BIN + vs", "gTRUE", "w17TRUE", H[is.finite(ych)])
row("run value (neg = better)",       "yrv ~ BIN + vs", "gTRUE", "w17TRUE", H[is.finite(yrv)])
row("whiff vs fastball-blind model",  "blind ~ BIN + vs", "gTRUE", "w17TRUE", H[is.finite(blind)])

cat("\n--- does it need separation? same gate on the low two-thirds ---\n")
A <- readRDS(file.path(MDIR, "final_gate_battery.rds"))$H
LO <- M[vs < VSCUT]
P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(ax) & is.finite(az) & is.finite(ax_diff) & is.finite(az_diff)]
P[, `:=`(mx = fifelse(p_throws == "L", -ax, ax), mz = az + 32.174,
         fx = fifelse(p_throws == "L", -(ax - ax_diff), ax - ax_diff), fz = (az - az_diff) + 32.174)]
AA <- P[, .(n = .N, ca = atan2(mean(mx), mean(mz))*180/pi,
            fa = atan2(mean(fx), mean(fz))*180/pi), by = .(pitcher, season)][n >= 40]
AA[, `:=`(id = as.character(pitcher), mv = abs(((ca - fa + 180) %% 360) - 180))]
LO <- merge(LO, AA[, .(id, season, mv)], by = c("id","season")); LO[, disc := mv - axis]
LO[, `:=`(gl = axis < 10 & arm >= arm_thr & disc <= 19.71, wl = axis < 11 & arm >= th55 & disc <= 27.5)]
cat(sprintf("%-42s %s   %s\n", "low-separation placebo",
            sp(crob(LO, y ~ gl + vs, "glTRUE")), sp(crob(LO, y ~ wl + vs, "wlTRUE"))))
cat(sprintf("%-42s %28d   %28d\n", "   (seasons in the placebo bin)", sum(LO$gl), sum(LO$wl)))

cat("\n--- permutation, membership shuffled within the high-separation third ---\n")
for (spec in c("g","w17")) {
  H[, gg := H[[spec]]]; obs <- crob(H, y ~ gg + vs, "ggTRUE")[1]; k <- sum(H$gg)
  ka <- uniqueN(H[gg == TRUE]$id); arms <- unique(H$id)
  p1 <- mean(replicate(4000, { H[, gp := FALSE]; H[sample(.N, k), gp := TRUE]
    crob(H, y ~ gp + vs, "gpTRUE")[1] }) >= obs)
  p2 <- mean(replicate(2000, { s <- sample(arms, ka); H[, gp := id %in% s]
    crob(H, y ~ gp + vs, "gpTRUE")[1] }) >= obs)
  cat(sprintf("  %-8s observed %+.3f | season-shuffle p = %.4f | arm-shuffle p = %.4f\n",
              spec, obs, p1, p2))
}

cat("\n--- leave one arm out ---\n")
for (spec in c("g","w17")) {
  H[, gg := H[[spec]]]; ar <- unique(H[gg == TRUE]$id)
  L <- rbindlist(lapply(ar, function(a) { x <- crob(H[id != a], y ~ gg + vs, "ggTRUE")
    data.table(drop = H[id == a, unique(nm)], est = x[1], p = x[3]) }))
  L2 <- rbindlist(apply(combn(length(ar), 2), 2, function(ix) { a <- ar[ix]
    x <- crob(H[!id %in% a], y ~ gg + vs, "ggTRUE"); data.table(p = x[3]) }))
  cat(sprintf("  %-8s LOO %d/%d hold p<.05 (worst %+.2f p %.4f, dropping %s) | L2O %d/%d hold\n",
              spec, sum(L$p < .05), nrow(L), L[which.max(p)]$est, max(L$p), L[which.max(p)]$drop,
              sum(L2$p < .05), nrow(L2)))
}

cat("\n--- who is in wide17 but not locked, and vice versa ---\n")
H[, tag := fifelse(g & w17, "both", fifelse(w17, "wide17 only", fifelse(g, "locked only", "")))]
print(H[tag != "", .(pitcher = nm, season, swings = nsw, sep = round(vs,1), axis = round(axis,1),
                     arm = round(arm,0), disc = round(disc,1), whiff = round(y,1),
                     rv = round(yrv,2), tag)][order(tag, -whiff)], row.names = FALSE)
saveRDS(H, file.path(MDIR, "wide17_compare.rds"))

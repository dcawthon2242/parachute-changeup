#!/usr/bin/env Rscript

# TWO QUESTIONS.
#
# (1) Which threshold combinations are actually strongest, stated honestly. The raw grid ranking is
#     misleading because many triples select the identical set of pitcher-seasons - once the
#     discrepancy cut is above the largest member value it stops binding, so a dozen "different"
#     cells are one cell wearing different labels. Dedupe on membership, then rank. And rank on the
#     t-statistic as well as the point estimate, because a six-season cell with a huge coefficient
#     is mostly telling you it has six seasons.
#
# (2) Whether "gradient not archetype" is a real claim or a hedge. The testable version: if these
#     pitches are a distinct type, membership should buy something beyond where the pitcher sits on
#     the underlying continuous scales. A step function should beat a smooth one. If instead the
#     smooth function absorbs the step, then the bin is a way of naming the tail, not a category.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
set.seed(43); options(width = 200); MDIR <- "data/statcast_model"
B <- readRDS(file.path(MDIR, "final_gate_battery.rds")); H <- B$H; setDT(H); G <- B$grid; setDT(G)
M <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(M); M <- M[league == "MLB"]

crob <- function(D, f, k, w = FALSE) {
  environment(f) <- environment()
  D <- D[complete.cases(D[, c(all.vars(f), "id", "nsw"), with = FALSE])]
  wt <- if (w) as.numeric(D$nsw) else rep(1, nrow(D))
  m <- lm(f, D, weights = wt); u <- residuals(m)*sqrt(wt); X <- model.matrix(m)*sqrt(wt)
  nc <- uniqueN(D$id); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k]); c(e, s, 2*pt(-abs(e/s), nc-1), e/s)
}

## =============================================================================================
cat("=== 1. THE GRID, DEDUPED ON WHO IT ACTUALLY SELECTS ===\n\n")
G[, memb := ""]
for (i in seq_len(nrow(G))) {
  th <- quantile(M$arm, G$ap[i], na.rm = TRUE)
  s <- H[axis < G$ax[i] & arm >= th & disc <= G$dc[i]]
  G[i, memb := paste(sort(paste0(s$nm, s$season)), collapse = "|")]
}
U <- G[, .SD[1], by = memb]  # one representative per distinct member set
cat(sprintf("  %d scored triples collapse to %d distinct pitcher-season sets.\n", nrow(G), nrow(U)))
cat(sprintf("  the top 6 by raw effect were %d labels for %d real set(s).\n\n",
            6, uniqueN(G[order(-est)][1:6]$memb)))

# refit each distinct set to get a t-statistic, and record the widest thresholds that produce it
info <- rbindlist(lapply(seq_len(nrow(U)), function(i) {
  sub <- G[memb == U$memb[i]]
  th <- quantile(M$arm, U$ap[i], na.rm = TRUE)
  H[, gg := axis < U$ax[i] & arm >= th & disc <= U$dc[i]]
  r <- crob(H, y ~ gg + vs, "ggTRUE"); rw <- crob(H, y ~ gg + vs, "ggTRUE", TRUE)
  rc <- crob(H, y ~ gg + vs + axis + arm, "ggTRUE", TRUE)
  data.table(n = sum(H$gg), arms = uniqueN(H[gg == TRUE]$id),
             ax_lo = min(sub$ax), ax_hi = max(sub$ax), ap = U$ap[i],
             dc_lo = min(sub$dc), dc_hi = max(sub$dc), ncell = nrow(sub),
             est = r[1], t = r[4], p = r[3], estw = rw[1], pw = rw[3], est_ctrl = rc[1], p_ctrl = rc[3],
             memb = U$memb[i])
}))
show <- function(D, lab, k = 3) {
  cat(sprintf("--- %s ---\n", lab))
  for (i in 1:k) {
    r <- D[i]
    cat(sprintf("  #%d  axis < %g%s | arm >= %.0fth pct (%.1f deg) | discrepancy <= %g%s\n",
                i, r$ax_lo, if (r$ax_hi != r$ax_lo) sprintf("-%g", r$ax_hi) else "",
                100*r$ap, quantile(M$arm, r$ap, na.rm = TRUE), r$dc_lo,
                if (r$dc_hi != r$dc_lo) sprintf(" (any cut up to %g)", r$dc_hi) else ""))
    cat(sprintf("      %d seasons / %d arms | unweighted %+.2f (t %.2f, p %.4f) | weighted %+.2f (p %.4f)\n",
                r$n, r$arms, r$est, r$t, r$p, r$estw, r$pw))
    cat(sprintf("      with axis+arm as continuous controls: %+.2f (p %.4f)\n", r$est_ctrl, r$p_ctrl))
    cat(sprintf("      members: %s\n\n", gsub("\\|", ", ", r$memb)))
  }
}
show(info[order(-est)], "STRONGEST BY POINT ESTIMATE (biased toward tiny sets)")
show(info[order(-t)],   "STRONGEST BY t-STATISTIC (what you should actually rank on)")
show(info[n >= 10][order(-t)], "STRONGEST AMONG SETS OF 10+ SEASONS")

cat("--- where the chosen spec sits among the 4 rankings ---\n")
own <- info[grepl("Skubal2020", memb) & n == 10 & abs(ap - .70) < 1e-9]
if (nrow(own)) {
  cat(sprintf("  by estimate: %d of %d | by t: %d of %d\n",
              which(info[order(-est)]$memb == own$memb[1]), nrow(info),
              which(info[order(-t)]$memb == own$memb[1]), nrow(info)))
}
cat(sprintf("\n  across the %d distinct sets: median t %.2f, %d%% reach p<.05, median size %d seasons\n",
            nrow(info), median(info$t), round(100*mean(info$p < .05)), as.integer(median(info$n))))
cat(sprintf("  how often does a given pitcher appear across all %d sets?\n", nrow(info)))
who <- sort(table(unlist(strsplit(paste(info$memb, collapse = "|"), "\\|"))), decreasing = TRUE)
print(data.table(season = names(who), sets = as.integer(who),
                 pct = sprintf("%.0f%%", 100*as.integer(who)/nrow(info)))[1:12], row.names = FALSE)

## =============================================================================================
cat("\n\n=== 2. GRADIENT OR ARCHETYPE? ===\n")
# The archetype claim says something changes at the boundary. The gradient claim says the bin is
# just the far end of a smooth slope. Four tests.

cat("\n-- (a) is there a smooth slope at all, across the whole high-separation third? --\n")
for (v in c("axis", "arm", "disc")) {
  r <- crob(H, as.formula(paste("y ~", v, "+ vs")), v, TRUE)
  cat(sprintf("   %-5s slope %+7.3f whiff pts per unit (se %.3f) p = %.4f\n", v, r[1], r[2], r[3]))
}
r <- crob(H, y ~ axis + arm + vs, c("axis"), TRUE)
r2 <- crob(H, y ~ axis + arm + vs, c("arm"), TRUE)
cat(sprintf("   jointly: axis %+.3f (p %.4f), arm %+.3f (p %.4f)\n", r[1], r[3], r2[1], r2[3]))

cat("\n-- (b) does the step beat the smooth function? nested model comparison --\n")
D <- H[complete.cases(H[, .(y, axis, arm, disc, vs, nsw, id)])]
wt <- as.numeric(D$nsw)
m_sm <- lm(y ~ axis + arm + disc + vs, D, weights = wt)
m_st <- lm(y ~ g + vs, D, weights = wt)
m_bo <- lm(y ~ axis + arm + disc + vs + g, D, weights = wt)
cat(sprintf("   smooth only (axis+arm+disc+sep)   R2 = %.4f\n", summary(m_sm)$r.squared))
cat(sprintf("   step only   (bin+sep)             R2 = %.4f\n", summary(m_st)$r.squared))
cat(sprintf("   both                              R2 = %.4f\n", summary(m_bo)$r.squared))
a <- anova(m_sm, m_bo)
cat(sprintf("   does the step add anything to the smooth?  F = %.2f, p = %.4f\n", a$F[2], a$`Pr(>F)`[2]))
a2 <- anova(m_st, m_bo)
cat(sprintf("   does the smooth add anything to the step?  F = %.2f, p = %.4f\n", a2$F[2], a2$`Pr(>F)`[2]))

cat("\n-- (c) is the boundary a cliff or a ramp? deciles of a continuous 'parachute score' --\n")
# One continuous index built from the same three quantities, standardised, higher = more parachute.
z <- function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
H[, score := -z(axis) + z(arm) - z(disc)]
H[, dec := cut(score, quantile(score, 0:10/10, na.rm = TRUE), include.lowest = TRUE, labels = 1:10)]
print(H[, .(seasons = .N, swings = sum(nsw), score = round(mean(score),2),
            whiff_over = round(weighted.mean(y, nsw),2), in_bin = sum(g)), by = dec][order(dec)],
      row.names = FALSE)
rs <- crob(H, y ~ score + vs, "score", TRUE)
cat(sprintf("\n   continuous score slope: %+.3f whiff pts per SD (se %.3f) p = %.4f\n", rs[1], rs[2], rs[3]))
rb <- crob(H, y ~ score + g + vs, "gTRUE", TRUE)
cat(sprintf("   bin membership ON TOP of the continuous score: %+.3f (se %.3f) p = %.4f\n", rb[1], rb[2], rb[3]))

cat("\n-- (d) bin members vs near-misses: is there a discontinuity? --\n")
# Near-miss = fails exactly one of the three conditions, and fails it narrowly.
tha <- quantile(M$arm, .70, na.rm = TRUE)
H[, nmiss := !g & (axis < 13 & arm >= tha - 4 & disc <= 25) ]
cat(sprintf("   %d bin members (score %.2f), %d near-misses (score %.2f), %d rest (score %.2f)\n",
            sum(H$g), H[g == TRUE, mean(score)], sum(H$nmiss), H[nmiss == TRUE, mean(score)],
            H[!g & !nmiss, .N], H[!g & !nmiss, mean(score)]))
cat(sprintf("   whiff above model: bin %+.2f | near-miss %+.2f | rest %+.2f\n",
            H[g == TRUE, weighted.mean(y, nsw)], H[nmiss == TRUE, weighted.mean(y, nsw)],
            H[!g & !nmiss, weighted.mean(y, nsw)]))
rn <- crob(H[!g == TRUE], y ~ nmiss + vs, "nmissTRUE", TRUE)
cat(sprintf("   near-misses vs the rest (bin excluded): %+.3f (se %.3f) p = %.4f\n", rn[1], rn[2], rn[3]))
rg <- crob(H[!nmiss == TRUE], y ~ g + vs, "gTRUE", TRUE)
cat(sprintf("   bin vs the rest (near-misses excluded): %+.3f (se %.3f) p = %.4f\n", rg[1], rg[2], rg[3]))

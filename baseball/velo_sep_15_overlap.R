#!/usr/bin/env Rscript

# THE TOP HIGH-SEPARATION OVERPERFORMERS LOOK LIKE PARACHUTE CHANGEUPS. IS THAT SOMETHING?
#
# Skubal, Cease and Ragans sit at the top of the high-separation overperformer list and at least two
# of them are locked-bin members. That is a striking overlap and it is also the exact shape of a
# trap: the names were selected BY the outcome, so noticing that they share a trait is not evidence
# the trait caused the outcome. The bin effect inside high separation was already measured at +5.175
# (p = .0017) and already known to rest on seven seasons from five arms.
#
# So this script does not re-measure the overlap. It asks the two questions that can actually
# distinguish a real archetype from two great pitchers who happen to match:
#
#   1. GATE SENSITIVITY. A real archetype degrades gracefully. Loosen the axis and efficiency gates
#      and membership should grow while the effect decays slowly, because the pitchers being
#      recruited are genuinely similar. Two lucky arms behave differently: the effect collapses the
#      moment anyone else is let in, because the estimate was never about the trait.
#
#   2. LEAVE-ONE-ARM-OUT. Drop each arm in turn and refit. If removing one or two names takes the
#      effect to nothing, the coefficient is a description of those names.
#
# Ragans is the useful new case. If he is a near-miss on one gate rather than a clear non-member,
# the gate is drawn too tight and the bin is undercounting the archetype. If he misses badly, the
# user's grouping of the three is the pattern-match rather than the data's.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 200); MDIR <- "data/statcast_model"
SPEC <- readRDS(file.path(MDIR, "locked_spec.rds"))
L <- SPEC$data; setDT(L); M <- L[league == "MLB"]
NM <- unique(readRDS(file.path(MDIR, "parachute_ff.rds"))[, .(pitcher, player_name)]); setDT(NM)
NM <- unique(NM, by = "pitcher")[, .(id = as.character(pitcher), nm = sub(",.*", "", player_name))]
M <- merge(M, NM, by = "id", all.x = TRUE)
thr <- quantile(M$vs, 2/3); M[, hi := vs >= thr]
cat(sprintf("MLB: %d seasons, %d arms | high-separation third = %.1f mph+ | locked gates: eff>=%.2f, axis<%d, arm>=p%d\n",
            nrow(M), uniqueN(M$id), thr, SPEC$spec$eff_min, SPEC$spec$axis_max,
            round(100*SPEC$spec$arm_pctile)))

## ---- 1. where do the three names actually sit? ------------------------------------------------------
cat("\n=== 1. gate-by-gate for the arms in question, plus every other high-sep overperformer ===\n")
cat("    pass/fail per gate; 'miss' is how far from qualifying on the failed gates\n\n")
show <- function(D) {
  D[, `:=`(g_ec = ec >= SPEC$spec$eff_min, g_ef = ef >= SPEC$spec$eff_min,
           g_ax = axis < SPEC$spec$axis_max, g_arm = arm >= arm_thr)]
  D[, gates := paste0(fifelse(g_ec,"E","."), fifelse(g_ef,"F","."),
                      fifelse(g_ax,"A","."), fifelse(g_arm,"S","."))]
  D[, miss := round(pmax(0, SPEC$spec$eff_min-ec) + pmax(0, SPEC$spec$eff_min-ef), 3)]
  D[, axis_over := round(pmax(0, axis - SPEC$spec$axis_max), 1)]
  D[order(-y), .(nm, season, sep = round(vs,1), whiff_over = round(y,1), ec = round(ec,3),
                 ef = round(ef,3), axis = round(axis,1), arm = round(arm,1),
                 gates, eff_short = miss, axis_over, bin)]
}
print(show(copy(M[nm %in% c("Skubal","Cease","Ragans")]))[order(nm, season)], row.names = FALSE)
cat("\n   legend: E = changeup efficiency gate, F = four-seam efficiency, A = axis gap, S = arm slot\n")
cat("\n   the 12 best high-separation seasons by whiff over the model, whether or not they are in the bin:\n")
print(show(copy(M[hi == TRUE]))[1:12], row.names = FALSE)

## ---- 2. gate sensitivity ---------------------------------------------------------------------------
crob <- function(D, f, k) {
  D <- D[complete.cases(D[, c(all.vars(f), "id"), with = FALSE])]
  if (uniqueN(D[[all.vars(f)[2]]]) < 2) return(c(NA, NA, NA))
  m <- lm(f, D); u <- residuals(m); X <- model.matrix(m); nc <- uniqueN(D$id)
  if (!k %in% colnames(X)) return(c(NA, NA, NA))
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k])
  c(est = e, se = s, p = 2*pt(-abs(e/s), nc-1))
}
cat("\n\n=== 2. gate sensitivity: does the effect decay gracefully as the definition widens? ===\n")
cat("    each row loosens one gate; the bin is rebuilt and refit inside the high-separation third\n\n")
G <- CJ(ax = c(8, 10, 14, 18, 25, 400), eff = c(0.92, 0.88, 0.85, 0.80, 0.00))
R <- rbindlist(lapply(seq_len(nrow(G)), function(i) {
  D <- copy(M[hi == TRUE])
  D[, b2 := ec >= G$eff[i] & ef >= G$eff[i] & axis < G$ax[i] & arm >= arm_thr]
  r <- crob(D, y ~ b2 + vs, "b2TRUE")
  data.table(axis_max = G$ax[i], eff_min = G$eff[i], seasons = sum(D$b2),
             arms = uniqueN(D[b2 == TRUE]$id), est = r[1], se = r[2], p = r[3])
}))
R[, `:=`(est = round(est,3), se = round(se,3), p = round(p,4))]
print(dcast(R, eff_min ~ axis_max, value.var = "est"), row.names = FALSE)
cat("\n   membership (seasons / arms) at each cell:\n")
print(dcast(R, eff_min ~ axis_max, value.var = "seasons"), row.names = FALSE)
cat("\n   p-values:\n")
print(dcast(R, eff_min ~ axis_max, value.var = "p"), row.names = FALSE)

## ---- 3. leave one arm out ---------------------------------------------------------------------------
cat("\n=== 3. leave-one-arm-out on the locked bin inside high separation ===\n")
H <- copy(M[hi == TRUE]); full <- crob(H, y ~ bin + vs, "binTRUE")
cat(sprintf("    full estimate: %+.3f (se %.3f, p = %.4f) on %d bin seasons from %d arms\n\n",
            full[1], full[2], full[3], sum(H$bin), uniqueN(H[bin == TRUE]$id)))
LO <- rbindlist(lapply(unique(H[bin == TRUE]$id), function(a) {
  r <- crob(H[id != a], y ~ bin + vs, "binTRUE")
  data.table(dropped = H[id == a, nm[1]], seasons_lost = H[id == a & bin == TRUE, .N],
             est = r[1], se = r[2], p = r[3]) }))
print(LO[order(est), .(dropped, seasons_lost, est = round(est,3), se = round(se,3),
                       p = round(p,4), still_sig = p < .05)], row.names = FALSE)
cat(sprintf("\n    effect survives dropping any single arm: %s\n",
            if (all(LO$p < .05)) "YES" else paste0("NO - fails when dropping ",
            paste(LO[p >= .05]$dropped, collapse = ", "))))

## ---- 4. is the overlap itself surprising? ------------------------------------------------------------
cat("\n=== 4. how surprising is the roster overlap, given the bin's size? ===\n")
K <- 10; top <- H[order(-y)][1:K]
q <- sum(top$bin); N <- nrow(H); Kb <- sum(H$bin)
cat(sprintf("    %d of the top %d high-sep seasons are bin members (bin is %d of %d seasons, %.1f%%)\n",
            q, K, Kb, N, 100*Kb/N))
cat(sprintf("    hypergeometric p(>= %d of top %d) = %.4f\n", q, K, phyper(q-1, Kb, N-Kb, K, lower.tail = FALSE)))
cat("    NOTE: this is a descriptive check, not an independent test - the same seasons drive the\n")
cat("    coefficient in section 3, so it cannot corroborate it.\n")

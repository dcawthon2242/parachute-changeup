#!/usr/bin/env Rscript

# TWO FOLLOW-UPS THE GATE SWEEP DEMANDS.
#
# 1. DILUTION BENCHMARK. The bin effect shrinking as the axis gate widens is not by itself evidence
#    of anything - widening always shrinks an estimate toward the population mean, purely
#    arithmetically, if the newly admitted members carry no signal. The informative comparison is
#    against that null. If the original 7 seasons carried everything and the 28 added at axis<18
#    carried nothing, the estimate would fall to 7/35 of its value. Anything above that line means
#    the marginal members are genuinely similar pitchers, which is what "real archetype" means.
#
# 2. LEAVE-TWO-OUT. Single-arm deletion already survived, but Skubal and Cease supply four of the
#    seven seasons between them and single deletion cannot rule out that the pair carries it.
#
# Also decomposed here: which gate is actually doing the work. The sweep hinted that loosening the
# efficiency floor STRENGTHENS the effect, which would mean the active spin requirements - half the
# original parachute definition - are not just inert but harmful.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 200); MDIR <- "data/statcast_model"
SPEC <- readRDS(file.path(MDIR, "locked_spec.rds")); L <- SPEC$data; setDT(L)
NM <- unique(readRDS(file.path(MDIR, "parachute_ff.rds"))[, .(pitcher, player_name)]); setDT(NM)
NM <- unique(NM, by = "pitcher")[, .(id = as.character(pitcher), nm = sub(",.*", "", player_name))]
M <- merge(L[league == "MLB"], NM, by = "id", all.x = TRUE)
H <- M[vs >= quantile(M$vs, 2/3)]

crob <- function(D, f, k) {
  D <- D[complete.cases(D[, c(all.vars(f), "id"), with = FALSE])]
  m <- lm(f, D); u <- residuals(m); X <- model.matrix(m); nc <- uniqueN(D$id)
  if (!k %in% colnames(X)) return(c(NA, NA, NA))
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k])
  c(est = e, se = s, p = 2*pt(-abs(e/s), nc-1))
}

## ---- 1. leave-two-out --------------------------------------------------------------------------
cat("=== 1. leave-two-arms-out, every pair of the five bin arms ===\n")
full <- crob(H, y ~ bin + vs, "binTRUE")
cat(sprintf("    full: %+.3f (se %.3f, p = %.4f)\n\n", full[1], full[2], full[3]))
arms <- unique(H[bin == TRUE]$id)
P <- rbindlist(apply(combn(length(arms), 2), 2, function(ix) {
  a <- arms[ix]; r <- crob(H[!id %in% a], y ~ bin + vs, "binTRUE")
  data.table(dropped = paste(sort(H[id %in% a, unique(nm)]), collapse = " + "),
             seasons_left = H[!id %in% a & bin == TRUE, .N],
             est = r[1], se = r[2], p = r[3]) }))
print(P[order(est), .(dropped, seasons_left, est = round(est,3), se = round(se,3),
                      p = round(p,4), sig = p < .05)], row.names = FALSE)
cat(sprintf("\n    survives every pair deletion: %s\n",
            if (all(P$p < .05, na.rm = TRUE)) "YES" else
            paste0("NO - fails on ", paste(P[p >= .05]$dropped, collapse = "; "))))

## ---- 2. dilution benchmark ---------------------------------------------------------------------
cat("\n=== 2. does widening the axis gate decay slower than pure dilution? ===\n")
cat("    pure-dilution prediction = (locked seasons / wider seasons) x locked estimate.\n")
cat("    that is what you would see if every newly admitted season carried ZERO effect.\n\n")
base <- crob(H, y ~ bin + vs, "binTRUE")[1]; n0 <- sum(H$bin)
DB <- rbindlist(lapply(c(8, 10, 14, 18, 25, 400), function(ax) {
  D <- copy(H); D[, b2 := ec >= .85 & ef >= .85 & axis < ax & arm >= arm_thr]
  r <- crob(D, y ~ b2 + vs, "b2TRUE")
  data.table(axis_max = ax, seasons = sum(D$b2), observed = r[1], p = r[3],
             dilution = base * n0 / sum(D$b2)) }))
DB[, ratio := observed / dilution]
print(DB[, .(axis_max, seasons, observed = round(observed,3), if_new_members_null = round(dilution,3),
             ratio = round(ratio,2), p = round(p,4))], row.names = FALSE)
cat("    ratio above 1 means the seasons admitted by the wider gate carry real effect of their own\n")

## ---- 3. which gate carries it -------------------------------------------------------------------
cat("\n=== 3. gate decomposition inside high separation: what is each requirement worth? ===\n")
combos <- list(
  "arm slot only"                 = quote(arm >= arm_thr),
  "axis < 10 only"                = quote(axis < 10),
  "efficiency only"               = quote(ec >= .85 & ef >= .85),
  "axis + arm"                    = quote(axis < 10 & arm >= arm_thr),
  "axis + efficiency"             = quote(axis < 10 & ec >= .85 & ef >= .85),
  "arm + efficiency"              = quote(arm >= arm_thr & ec >= .85 & ef >= .85),
  "LOCKED: all three"             = quote(axis < 10 & arm >= arm_thr & ec >= .85 & ef >= .85))
G <- rbindlist(lapply(names(combos), function(k) {
  D <- copy(H); D[, b2 := eval(combos[[k]])]
  r <- crob(D, y ~ b2 + vs, "b2TRUE")
  data.table(definition = k, seasons = sum(D$b2), arms = uniqueN(D[b2 == TRUE]$id),
             est = r[1], se = r[2], p = r[3]) }))
print(G[order(-est), .(definition, seasons, arms, est = round(est,3), se = round(se,3),
                       p = round(p,4))], row.names = FALSE)

## ---- 4. and does the axis gate work OUTSIDE high separation? --------------------------------------
cat("\n=== 4. the same definitions in the bottom two thirds of separation ===\n")
Lo <- M[vs < quantile(M$vs, 2/3)]
G2 <- rbindlist(lapply(names(combos), function(k) {
  D <- copy(Lo); D[, b2 := eval(combos[[k]])]
  r <- crob(D, y ~ b2 + vs, "b2TRUE")
  data.table(definition = k, seasons = sum(D$b2), est = r[1], se = r[2], p = r[3]) }))
print(G2[, .(definition, seasons, est = round(est,3), se = round(se,3), p = round(p,4))],
      row.names = FALSE)

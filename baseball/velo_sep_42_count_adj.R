#!/usr/bin/env Rscript

# COUNT-ADJUSTED RUN VALUE FOR THE CHANGEUP ARCHETYPES.
#
# Raw run value per 100 mixes two things: how well the pitch performs inside a given count, and how
# favourable the counts are that it gets thrown in. A pitch used mostly with two strikes will look
# good on raw run value even if it is ordinary, because two-strike pitches are worth more to the
# pitcher on average. The archetypes do not all share a count profile, so the raw comparison is not
# clean.
#
# Three adjustments, in increasing strictness:
#   COUNT        - standardise to the league changeup count distribution (12 cells)
#   COUNT + HAND - add batter handedness, since these archetypes skew opposite-handed (24 cells)
#   + STUFF      - on top of that, subtract the stuff model's prediction, so what is left is
#                  performance beyond shape, beyond count, beyond platoon
#
# Cell means are computed leave-one-arm-out so a pitcher is never benchmarked against himself.
# Then an Oaxaca split says how much of each archetype's raw edge is the counts and how much is
# the pitch.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 240); MDIR <- "data/statcast_model"

P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(rv)][, id := as.character(pitcher)]
P[, cnt := paste0(balls, "-", strikes)]
P[, cell2 := paste0(cnt, "|", same_hand)]

# the stuff model was fit on a slightly narrower row set, so its predictions cannot be attached
# pitch by pitch. that does not matter here: the model carries no count features, so its prediction
# is constant with respect to count and can be applied as a season-level offset instead.
OOF <- readRDS(file.path(MDIR, "absorb_preds.rds"))$STUFF$r
PRED <- OOF[, .(pred = 100*mean(p)), by = .(id, season)]
LGPRED <- OOF[, 100*mean(p)]

M <- readRDS(file.path(MDIR, "archetype_roster.rds")); setDT(M)
TAGS <- c("A2_wide","A4_broad","A5_seam","A6_extreme","A7_trad","A8_mismatch")
LAB <- c(A2_wide = "Matched axis, wide", A4_broad = "Broad (axis + separation)",
         A5_seam = "Seam-shifted matched", A6_extreme = "Extreme seam shift",
         A7_trad = "Traditional, low seam dev", A8_mismatch = "Mismatch (underperforms)")
P <- merge(P, M[, c("id","season","nm","nsw", TAGS), with = FALSE], by = c("id","season"))
P <- P[nsw >= 75]
cat(sprintf("%d changeups, %d pitcher-seasons, %d arms.\n\n", nrow(P),
            uniqueN(P[, .(id, season)]), uniqueN(P$id)))

## ---- leave-one-arm-out cell means ---------------------------------------------------------
loo <- function(D, cellcol, val) {
  tot <- D[, .(S = sum(get(val)), N = .N), by = c(cellcol)]
  arm <- D[, .(s = sum(get(val)), n = .N), by = c(cellcol, "id")]
  z <- merge(arm, tot, by = cellcol)
  z[, m := (S - s)/pmax(N - n, 1)]
  merge(D, z[, c(cellcol, "id", "m"), with = FALSE], by = c(cellcol, "id"), sort = FALSE)$m
}
P[, m1 := loo(P, "cnt", "rv")]
P[, m2 := loo(P, "cell2", "rv")]
P[, `:=`(a_cnt = rv - m1, a_hand = rv - m2)]
# stuff enters as a season-level offset: how much better the shape was predicted to be than average
P <- merge(P, PRED, by = c("id","season"), all.x = TRUE)
P[, a_stuff := a_hand - (pred - LGPRED)/100]
cat(sprintf("stuff prediction attached for %.1f%% of pitches; the rest fall outside the model's\n",
            100*mean(is.finite(P$pred))))
cat("feature filter and are dropped from the stuff-adjusted column only.\n\n")

## ---- how different are the count profiles at all? ------------------------------------------
cat("=== 1. DO THE ARCHETYPES SIT IN DIFFERENT COUNTS? ===\n\n")
CNTS <- c("0-0","1-0","2-0","3-0","0-1","1-1","2-1","3-1","0-2","1-2","2-2","3-2")
CNTS <- CNTS[CNTS %in% unique(P$cnt)]
P[, cnt := factor(cnt, levels = CNTS)]
LG <- P[, .(w = .N/nrow(P), rv = 100*mean(rv)), by = cnt][order(cnt)]
SHR <- rbindlist(c(
  list(data.table(grp = "league", dcast(P[, .N, by = cnt][, s := 100*N/sum(N)], . ~ cnt,
                                        value.var = "s")[, -1])),
  lapply(TAGS, function(t) data.table(grp = LAB[[t]],
    dcast(P[get(t) == TRUE, .N, by = cnt][, s := 100*N/sum(N)], . ~ cnt, value.var = "s")[, -1]))),
  fill = TRUE)
print(SHR[, lapply(.SD, function(z) if (is.numeric(z)) round(z,1) else z)], row.names = FALSE)
cat("\n  league run value per 100 in each count, which is what the mix is worth:\n")
print(dcast(LG, . ~ cnt, value.var = "rv")[, -1][, lapply(.SD, round, 2)], row.names = FALSE)
cat(sprintf("\n  two-strike share: league %.1f%%", 100*P[, mean(strikes == 2)]))
for (t in TAGS) cat(sprintf(" | %s %.1f%%", sub(" .*", "", LAB[[t]]), 100*P[get(t) == TRUE, mean(strikes == 2)]))
cat("\n")

## ---- the adjusted numbers -------------------------------------------------------------------
crob <- function(D, yv, gv) {
  D <- D[is.finite(get(yv))]
  A <- D[, .(n = .N, y = mean(get(yv)), g = get(gv)[1]), by = .(id, season)]
  w <- as.numeric(A$n); m <- lm(y ~ g, A, weights = w)
  u <- residuals(m)*sqrt(w); X <- model.matrix(m)*sqrt(w); nc <- uniqueN(A$id)
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, A$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[2]); s <- unname(sqrt(diag(V))[2]); c(100*e, 2*pt(-abs(e/s), nc-1))
}
st <- function(r) sprintf("%+7.3f%s", r[1], ifelse(r[2] < .01, "**", ifelse(r[2] < .05, "* ", "  ")))

cat("\n=== 2. RAW AGAINST COUNT-ADJUSTED, RUN VALUE PER 100 ===\n\n")
cat("  each column is the archetype minus everyone else, arm-clustered, weighted by pitches.\n\n")
cat(sprintf("  %-28s %8s %11s %11s %11s %11s\n", "archetype", "seasons", "raw",
            "count adj", "+ hand adj", "+ stuff"))
RES <- rbindlist(lapply(TAGS, function(t) {
  r0 <- crob(P, "rv", t); r1 <- crob(P, "a_cnt", t); r2 <- crob(P, "a_hand", t)
  r3 <- crob(P, "a_stuff", t)
  cat(sprintf("  %-28s %8d %11s %11s %11s %11s\n", LAB[[t]], uniqueN(P[get(t) == TRUE, .(id, season)]),
              st(r0), st(r1), st(r2), st(r3)))
  data.table(tag = t, raw = r0[1], p_raw = r0[2], cnt = r1[1], p_cnt = r1[2],
             hand = r2[1], p_hand = r2[2], stuff = r3[1], p_stuff = r3[2])
}))

cat("\n  and the archetypes' own levels, not differences:\n\n")
cat(sprintf("  %-28s %10s %12s %12s\n", "archetype", "pitches", "raw rv/100", "count+hand adj"))
for (t in c(TAGS, "ALL")) {
  X <- if (t == "ALL") P else P[get(t) == TRUE]
  cat(sprintf("  %-28s %10d %12.3f %12.3f\n", if (t == "ALL") "everyone" else LAB[[t]],
              nrow(X), 100*mean(X$rv), 100*mean(X$a_hand)))
}

## ---- Oaxaca split ---------------------------------------------------------------------------
cat("\n=== 3. HOW MUCH OF THE RAW EDGE IS THE COUNTS? ===\n\n")
cat("  raw edge split into the part from performing better inside counts and the part from\n")
cat("  sitting in better counts. the two add to the raw edge.\n\n")
cat(sprintf("  %-28s %11s %14s %14s %10s\n", "archetype", "raw edge", "within count",
            "count mix", "mix share"))
for (t in TAGS) {
  X <- P[get(t) == TRUE]; Y <- P[get(t) == FALSE]
  a <- X[, .(w = .N/nrow(X), m = mean(rv)), by = cnt]
  b <- Y[, .(w = .N/nrow(Y), m = mean(rv)), by = cnt]
  z <- merge(a, b, by = "cnt", suffixes = c("_a","_b"), all = TRUE)
  z[is.na(z)] <- 0
  raw <- 100*(sum(z$w_a*z$m_a) - sum(z$w_b*z$m_b))
  within <- 100*sum(z$w_b*(z$m_a - z$m_b))
  mix <- 100*sum((z$w_a - z$w_b)*z$m_a)
  cat(sprintf("  %-28s %11.3f %14.3f %14.3f %9.0f%%\n", LAB[[t]], raw, within, mix,
              100*mix/(abs(within) + abs(mix))))
}

## ---- count-adjusted season table --------------------------------------------------------------
cat("\n=== 4. THE ADJUSTED LEADERBOARD ===\n\n")
S <- P[, .(np = .N, nm = nm[1], raw = 100*mean(rv), adj = 100*mean(a_hand),
           adjs = 100*mean(a_stuff, na.rm = TRUE), two = 100*mean(strikes == 2),
           opp = 100*mean(same_hand == 0)), by = .(id, season)]
S <- merge(S, M[, c("id","season", TAGS), with = FALSE], by = c("id","season"))
S[, tag := apply(as.matrix(S[, ..TAGS]), 1, function(r) {
  k <- sub("^A[0-9]_", "", TAGS[which(r == TRUE)]); if (!length(k)) "" else paste(k, collapse = "/") })]
S <- S[np >= 200]
cat(sprintf("  %d pitcher-seasons with 200+ changeups. shift is adjusted minus raw.\n", nrow(S)))
cat("\n  fifteen best by count-and-hand adjusted run value:\n\n")
print(S[order(-adj)][1:15, .(pitcher = nm, season, changeups = np, raw = round(raw,2),
        adjusted = round(adj,2), shift = round(adj - raw,2), `+stuff` = round(adjs,2),
        `2str %` = round(two,1), `opp %` = round(opp,1), archetype = tag)], row.names = FALSE)
cat("\n  the fifteen whose ranking the adjustment changes most (adjusted well below raw):\n\n")
print(S[order(adj - raw)][1:15, .(pitcher = nm, season, changeups = np, raw = round(raw,2),
        adjusted = round(adj,2), shift = round(adj - raw,2), `2str %` = round(two,1),
        archetype = tag)], row.names = FALSE)
cat(sprintf("\n  correlation of raw and adjusted across seasons: %+.3f. mean absolute shift %.3f runs.\n",
            S[, cor(raw, adj)], S[, mean(abs(adj - raw))]))

cat("\n=== 5. DOES THE ADJUSTMENT MAKE THE EDGE MORE REPEATABLE? ===\n\n")
for (v in c("raw","adj","adjs")) {
  X <- S[, .(id, season, r = get(v), n = np)]
  Y <- copy(X)[, season := season - 1][, .(id, season, r2 = r)]
  J <- merge(X, Y, by = c("id","season"))
  ct <- cor.test(J$r, J$r2)
  cat(sprintf("  %-24s %3d pairs   sd %5.2f   year to year %+.3f   p %.4f\n",
              c(raw = "raw", adj = "count + hand adj", adjs = "+ stuff")[[v]],
              nrow(J), sd(X$r), ct$estimate, ct$p.value))
}
saveRDS(list(season = S, arch = RES), file.path(MDIR, "count_adjusted_rv.rds"))

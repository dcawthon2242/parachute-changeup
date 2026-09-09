#!/usr/bin/env Rscript

# IS THIS PITCH BEING MISUSED?
#
# Two specific claims to test: that it should be thrown more often below the zone for chase, and
# that it is being thrown in the wrong counts.
#
# "Wrong" needs a counterfactual. A changeup thrown is a hard pitch not thrown, so the quantity that
# matters is not the changeup's run value but the GAP between the changeup and the same pitcher's
# hard pitches in the same count. If that gap is large in a count where the changeup is rarely
# thrown, there are runs sitting on the table.
#
# The arsenal file carries four-seam, sinker, cutter and changeup but no breaking balls, so the
# alternative here is "his hard stuff", not "everything else". That understates the alternative for
# arms with a good slider and the usage shares are shares of hard pitches, not of all pitches.
#
# The allocation test at the end is the one that decides it: with pitcher fixed effects, does the
# gap shrink as a pitcher leans on the changeup harder? If it stays positive at high usage the
# pitch is underthrown; if it decays to zero the pitcher has already found the optimum.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 235); MDIR <- "data/statcast_model"

A <- readRDS(file.path(MDIR, "parachute_rv.rds")); setDT(A)
A <- A[is.finite(rv) & pitch_type %in% c("FF","SI","FC","CH")]
A[, `:=`(id = as.character(pitcher), ch = pitch_type == "CH", cnt = paste0(balls, "-", strikes))]
M <- readRDS(file.path(MDIR, "archetype_roster.rds")); setDT(M)
M[, para := A4_broad | A5_seam]
A <- merge(A, M[, .(id, season, nm, para, A4_broad, A5_seam, nsw, vs)], by = c("id","season"))
A <- A[nsw >= 75]
cat(sprintf("%d pitches, %d pitcher-seasons, %d arms. %d seasons carry a parachute flag.\n",
            nrow(A), uniqueN(A[, .(id, season)]), uniqueN(A$id), uniqueN(A[para == TRUE, .(id, season)])))

CNTS <- c("0-0","1-0","2-0","3-0","0-1","1-1","2-1","3-1","0-2","1-2","2-2","3-2")
A[, cnt := factor(cnt, levels = CNTS)]

## =============================================================================================
cat("\n=== 1. THE GAP BY COUNT: CHANGEUP MINUS HIS OWN HARD PITCHES ===\n\n")
bycnt <- function(X, lab) {
  G <- X[, .(n = .N, rv = 100*mean(rv)), by = .(id, season, cnt, ch)]
  G <- dcast(G, id + season + cnt ~ ch, value.var = c("n","rv"))
  setnames(G, c("n_FALSE","n_TRUE","rv_FALSE","rv_TRUE"), c("nh","nc","rvh","rvc"))
  G <- G[!is.na(nc) & !is.na(nh) & nc >= 5 & nh >= 15]
  G[, `:=`(gap = rvc - rvh, use = 100*nc/(nc + nh))]
  R <- G[, .(seasons = .N, ch = sum(nc), use = 100*sum(nc)/sum(nc + nh),
             rvc = weighted.mean(rvc, nc), rvh = weighted.mean(rvh, nh),
             gap = weighted.mean(gap, nc),
             se = sqrt(sum(nc^2 * (gap - weighted.mean(gap, nc))^2)/sum(nc)^2)), by = cnt][order(cnt)]
  R[, `:=`(t = gap/se, grp = lab)]
  R
}
ALL <- bycnt(A, "all arms")
print(ALL[, .(count = cnt, seasons, changeups = ch, `CH use %` = round(use,1),
              `CH rv/100` = round(rvc,2), `hard rv/100` = round(rvh,2),
              gap = round(gap,2), t = round(t,1))], row.names = FALSE)
cat(sprintf("\n  correlation across counts between the gap and how often the changeup is thrown: %+.3f\n",
            ALL[, cor(gap, use)]))
cat("  (positive means usage already tracks value; negative means it is pointed the wrong way)\n")

cat("\n  parachute seasons only:\n\n")
PAR <- bycnt(A[para == TRUE], "parachute")
print(PAR[, .(count = cnt, seasons, changeups = ch, `CH use %` = round(use,1),
              `CH rv/100` = round(rvc,2), `hard rv/100` = round(rvh,2),
              gap = round(gap,2), t = round(t,1))], row.names = FALSE)
cat(sprintf("\n  same correlation, parachute seasons: %+.3f\n", PAR[, cor(gap, use)]))

cat("\n  collapsed to strike count, which is where the usage decision actually lives:\n\n")
A[, sc := paste0(strikes, " strikes")]
byS <- function(X, lab) {
  G <- X[, .(n = .N, rv = 100*mean(rv)), by = .(id, season, sc, ch)]
  G <- dcast(G, id + season + sc ~ ch, value.var = c("n","rv"))
  setnames(G, c("n_FALSE","n_TRUE","rv_FALSE","rv_TRUE"), c("nh","nc","rvh","rvc"))
  G <- G[!is.na(nc) & !is.na(nh) & nc >= 10 & nh >= 25][, gap := rvc - rvh]
  m <- G[, .(seasons = .N, ch = sum(nc), use = 100*sum(nc)/sum(nc + nh),
             rvc = weighted.mean(rvc, nc), rvh = weighted.mean(rvh, nh),
             gap = weighted.mean(gap, nc)), by = sc][order(sc)]
  # arm-clustered test of the gap against zero
  m[, p := sapply(sc, function(s) {
    D <- G[sc == s]; w <- as.numeric(D$nc); f <- lm(gap ~ 1, D, weights = w)
    u <- residuals(f)*sqrt(w); X2 <- model.matrix(f)*sqrt(w); nc2 <- uniqueN(D$id)
    b <- solve(crossprod(X2)); V <- b %*% crossprod(rowsum(X2*u, D$id)) %*% b * (nc2/(nc2-1))
    2*pt(-abs(coef(f)[1]/sqrt(diag(V))[1]), nc2-1) })]
  m[, grp := lab][]
}
print(rbind(byS(A, "all arms"), byS(A[para == TRUE], "parachute"))[
  , .(grp, strikes = sc, seasons, changeups = ch, `CH use %` = round(use,1),
      `CH rv/100` = round(rvc,2), `hard rv/100` = round(rvh,2), gap = round(gap,2),
      p = signif(p,3))], row.names = FALSE)

## =============================================================================================
cat("\n=== 2. LOCATION: IS IT THROWN BELOW THE ZONE OFTEN ENOUGH? ===\n\n")
P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(rv)][, id := as.character(pitcher)]
P <- merge(P, M[, .(id, season, nm, para, nsw)], by = c("id","season"))[nsw >= 75]
P[, wide := abs(plate_x) > 0.83]
P[, loc := fifelse(!wide & z_rel_bot >= 0 & z_rel_top <= 0, "in zone",
           fifelse(z_rel_bot >= -0.40 & z_rel_bot < 0 & !wide, "just below",
           fifelse(z_rel_bot >= -1.00 & z_rel_bot < -0.40 & !wide, "chase, below",
           fifelse(z_rel_bot < -1.00 & !wide, "buried",
           fifelse(z_rel_top > 0, "elevated", "wide")))))]
LV <- c("in zone","just below","chase, below","buried","elevated","wide")
P[, loc := factor(loc, levels = LV)]
LOC <- P[, .(pitches = .N, share = 100*.N/nrow(P), swing = 100*mean(is_swing == 1),
             whiff = 100*sum(whiff == 1)/pmax(sum(is_swing == 1),1),
             rv = 100*mean(rv)), by = loc][order(loc)]
print(LOC[, .(location = loc, pitches, `share %` = round(share,1), `swing %` = round(swing,1),
              `whiff/swing %` = round(whiff,1), `rv/100` = round(rv,2))], row.names = FALSE)

cat("\n  by strike count, run value per 100 in each location (all arms):\n\n")
W <- dcast(P[, .(rv = 100*mean(rv), n = .N), by = .(strikes, loc)], loc ~ strikes,
           value.var = c("rv","n"))
print(W[, .(location = loc, `0str rv` = round(rv_0,2), `1str rv` = round(rv_1,2),
            `2str rv` = round(rv_2,2), `0str n` = n_0, `1str n` = n_1, `2str n` = n_2)],
      row.names = FALSE)

cat("\n  share of changeups in each location, by strike count and by group:\n\n")
SH <- P[, .(n = .N), by = .(grp = fifelse(para, "parachute", "rest"), strikes, loc)]
SH[, share := 100*n/sum(n), by = .(grp, strikes)]
print(dcast(SH, grp + strikes ~ loc, value.var = "share")[
  , lapply(.SD, function(z) if (is.numeric(z)) round(z,1) else z)], row.names = FALSE)

cat("\n  parachute minus rest, share points in each location (2 strikes only):\n\n")
D2 <- dcast(SH[strikes == 2], loc ~ grp, value.var = "share")
D2[, diff := parachute - rest]
print(D2[, .(location = loc, parachute = round(parachute,1), rest = round(rest,1),
             diff = round(diff,1))], row.names = FALSE)

## =============================================================================================
cat("\n=== 3. DOES THE EDGE SURVIVE MORE USAGE? ===\n\n")
cat("  pitcher fixed effects: within an arm's own career, do the seasons where he leaned on the\n")
cat("  changeup harder show a smaller gap over his hard pitches? a flat line means the marginal\n")
cat("  changeup is worth as much as the average one, which is what underuse looks like.\n\n")
U <- A[, .(nc = sum(ch), nh = sum(!ch), rvc = 100*mean(rv[ch]), rvh = 100*mean(rv[!ch])),
       by = .(id, season, nm, para)]
U <- U[nc >= 150 & nh >= 400][, `:=`(use = 100*nc/(nc + nh), gap = rvc - rvh)]
U[, nseas := .N, by = id]
cat(sprintf("  %d pitcher-seasons with 150+ changeups and 400+ hard pitches, %d arms;\n",
            nrow(U), uniqueN(U$id)))
cat(sprintf("  %d arms have 3+ such seasons and carry the within-arm identification.\n\n",
            uniqueN(U[nseas >= 3]$id)))
fe <- function(D, lab) {
  D <- copy(D)[nseas >= 3]
  D[, `:=`(use_c = use - mean(use), gap_c = gap - mean(gap)), by = id]
  w <- as.numeric(D$nc); m <- lm(gap_c ~ use_c, D, weights = w)
  u <- residuals(m)*sqrt(w); X <- model.matrix(m)*sqrt(w); nc <- uniqueN(D$id)
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- coef(m)[2]; s <- sqrt(diag(V))[2]
  cat(sprintf("  %-22s %4d seasons %3d arms   slope %+6.3f runs/100 per usage point  p %.3f\n",
              lab, nrow(D), nc, e, 2*pt(-abs(e/s), nc-1)))
  invisible(e)
}
fe(U, "all arms"); fe(U[para == TRUE], "parachute arms")
cat("\n  raw picture: gap by usage quintile (between arms, so confounded by who throws a lot)\n\n")
U[, q := cut(use, quantile(use, 0:5/5), include.lowest = TRUE, labels = paste0("Q",1:5))]
print(U[, .(seasons = .N, `CH use %` = round(mean(use),1), `CH rv/100` = round(weighted.mean(rvc, nc),2),
            `hard rv/100` = round(weighted.mean(rvh, nh),2),
            gap = round(weighted.mean(gap, nc),2)), by = q][order(q)], row.names = FALSE)

## =============================================================================================
cat("\n=== 4. THE TARGET BAND, IN FINE SLICES ===\n\n")
cat("  vertical position relative to the bottom of the strike zone, in three-inch bands.\n")
cat("  negative is below the knees. horizontal misses excluded so this is the vertical story only.\n\n")
V <- P[abs(plate_x) <= 0.83 & z_rel_bot > -1.75 & z_rel_bot < 1.75]
V[, band := cut(z_rel_bot, seq(-1.75, 1.75, 0.25), include.lowest = TRUE)]
VB <- V[, .(n = .N, rv = 100*mean(rv)), by = .(band, strikes)]
VB[, share := 100*n/sum(n), by = strikes]
BW <- dcast(VB, band ~ strikes, value.var = c("rv","share","n"))
BW[, lo := 12*(-1.75 + 0.25*(as.numeric(band) - 1))]
print(BW[, .(`inches vs knees` = sprintf("%+.0f to %+.0f", lo, lo + 3),
             `0str rv` = round(rv_0,2), `1str rv` = round(rv_1,2), `2str rv` = round(rv_2,2),
             `0str share` = round(share_0,1), `1str share` = round(share_1,1),
             `2str share` = round(share_2,1), `2str n` = n_2)], row.names = FALSE)
cat("\n  share of changeups thrown more than 9 inches below the knees, where every count is negative:\n")
print(V[, .(`0 strikes` = round(100*mean(z_rel_bot < -0.75 & strikes == 0)/mean(strikes == 0),1),
            `1 strike` = round(100*mean(z_rel_bot < -0.75 & strikes == 1)/mean(strikes == 1),1),
            `2 strikes` = round(100*mean(z_rel_bot < -0.75 & strikes == 2)/mean(strikes == 2),1)),
        by = .(group = fifelse(para, "parachute", "rest"))], row.names = FALSE)

## =============================================================================================
cat("\n=== 5. WHAT WOULD REALLOCATION BE WORTH? ===\n\n")
# hold each season's total changeup count fixed and shift changeups out of no-strike counts into
# two-strike counts, capped at the 90th percentile of two-strike usage already observed. the gap
# used is the league profile by strike count, which is what a coach could actually act on.
SG <- byS(A, "all")[, .(sc, gap)]
CAP2 <- U[, quantile(100*0, .5)]                       # placeholder, real cap computed below
SU <- A[, .(nc = sum(ch), n = .N), by = .(id, season, sc)]
SU[, use := 100*nc/n]
cap2 <- SU[sc == "2 strikes", quantile(use, .90)]
cat(sprintf("  two-strike changeup usage: median %.1f%%, 90th percentile %.1f%%.\n",
            SU[sc == "2 strikes", median(use)], cap2))
W2 <- dcast(SU, id + season ~ sc, value.var = c("nc","n","use"))
setnames(W2, names(W2), gsub("[ _]", "", names(W2)))
W2 <- W2[complete.cases(W2)]
SGP <- byS(A[para == TRUE], "para")[, .(sc, gap)]
g0 <- SG[sc == "0 strikes", gap]; g2 <- SG[sc == "2 strikes", gap]
p0 <- SGP[sc == "0 strikes", gap]; p2 <- SGP[sc == "2 strikes", gap]
W2[, room := pmax(cap2/100*n2strikes - nc2strikes, 0)]
W2[, move := pmin(room, nc0strikes)]                   # cannot move more than he throws early
W2[, tot := n0strikes + n1strikes + n2strikes]
W2 <- merge(W2, unique(A[, .(id, season, nm, para)]), by = c("id","season"))
# each group priced at its own gap profile, since the parachute edge is far more count-dependent
W2[, gain := move*fifelse(para, p2 - p0, g2 - g0)/100]
cat(sprintf("  moving one changeup from a no-strike to a two-strike count is worth %+.3f runs for a\n",
            (g2 - g0)/100))
cat(sprintf("  typical arm and %+.3f runs for a parachute arm, whose edge is far more count-dependent.\n",
            (p2 - p0)/100))
print(W2[, .(seasons = .N, `changeups moved` = round(mean(move)),
             `runs gained per season` = round(mean(gain),2),
             `runs per 100 pitches` = round(100*mean(gain)/mean(tot),3)),
         by = .(group = fifelse(para, "parachute", "rest"))], row.names = FALSE)
cat("\n  the ten seasons with the most room, by this rule:\n\n")
print(W2[order(-gain)][1:10, .(pitcher = nm, season, changeups = nc0strikes + nc1strikes + nc2strikes,
        `0str CH` = nc0strikes, `2str CH` = nc2strikes, `2str use %` = round(use2strikes,1),
        move = round(move), `runs gained` = round(gain,2))], row.names = FALSE)
saveRDS(list(bycount = ALL, para = PAR, loc = LOC, band = BW, usage = U, realloc = W2),
        file.path(MDIR, "usage_study.rds"))

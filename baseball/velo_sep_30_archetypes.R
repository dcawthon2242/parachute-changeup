#!/usr/bin/env Rscript

# EVERY ARCHETYPE THE STUDY HAS DEFINED, WITH ITS ROSTER.
#
# Nine definitions have accumulated across this work, some superseded, some still standing, one
# that predicts underperformance. Assembling them in one place makes two things visible that the
# separate runs could not: which pitcher-seasons satisfy several definitions at once, and which
# definitions are really the same population under different names.
#
# Everything is rebuilt from source at a 40-swing floor so the superseded definitions can be shown
# as they were originally stated, with the 75-swing versions marked separately.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 215); MDIR <- "data/statcast_model"

L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L); M <- L[league == "MLB"]
F <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F)
NM <- unique(F[, .(pitcher, player_name)])[, .(id = as.character(pitcher), nm = player_name)]
M <- merge(M, unique(NM, by = "id"), by = "id", all.x = TRUE)
AD <- F[is.finite(arm_diff), .(armdiff = mean(arm_diff)), by = .(pitcher, season)][, id := as.character(pitcher)]
M <- merge(M, AD[, .(id, season, armdiff)], by = c("id","season"), all.x = TRUE)
SD <- readRDS(file.path(MDIR, "ch_seam_deviation.rds")); setDT(SD)
M <- merge(M, SD[, .(id, season, sswsign = dev, ssw = absdev)], by = c("id","season"), all.x = TRUE)
P <- F[is.finite(ax) & is.finite(az) & is.finite(ax_diff) & is.finite(az_diff)]
P[, `:=`(mx = fifelse(p_throws == "L", -ax, ax), mz = az + 32.174,
         fx = fifelse(p_throws == "L", -(ax - ax_diff), ax - ax_diff), fz = (az - az_diff) + 32.174)]
DS <- P[, .(n = .N, ca = atan2(mean(mx), mean(mz))*180/pi,
            fa = atan2(mean(fx), mean(fz))*180/pi), by = .(pitcher, season)][n >= 40]
DS[, `:=`(id = as.character(pitcher), mv = abs(((ca - fa + 180) %% 360) - 180))]
M <- merge(M, DS[, .(id, season, mv)], by = c("id","season"), all.x = TRUE); M[, disc := mv - axis]
CH <- readRDS(file.path(MDIR, "mlb_chase_resid.rds")); setDT(CH)
CH <- CH[, .(ych = 100*mean(r_all, na.rm = TRUE)), by = .(pitcher, season)][, id := as.character(pitcher)]
VR <- readRDS(file.path(MDIR, "velo_sep_resid.rds"))
RV <- as.data.table(VR$R)[, .(yrv = 100*mean(rv - q_aware)), by = .(pitcher, season)][, id := as.character(pitcher)]
WB <- as.data.table(VR$W)[, .(blind = 100*mean(whiff - p_blind)), by = .(pitcher, season)][, id := as.character(pitcher)]
M <- merge(M, CH[, .(id, season, ych)], by = c("id","season"), all.x = TRUE)
M <- merge(M, RV[, .(id, season, yrv)], by = c("id","season"), all.x = TRUE)
M <- merge(M, WB[, .(id, season, blind)], by = c("id","season"), all.x = TRUE)
M[, `:=`(effgap = ef - ec, slot = arm)]

# thresholds, each stated in the units the definition was written in
SEP3 <- quantile(M$vs, 2/3)                       # top third of separation, 40-swing population
A75  <- M[nsw >= 75]
TH55 <- quantile(A75$arm, .55, na.rm = TRUE)
QAX  <- quantile(A75$axis, .25); QVS <- quantile(A75$vs, .60)
QSSW25 <- quantile(A75$ssw, .25, na.rm = TRUE); QSSW75 <- quantile(A75$ssw, .75, na.rm = TRUE)
QSGN75 <- quantile(A75$sswsign, .75, na.rm = TRUE); MEDSSW <- median(A75$ssw, na.rm = TRUE)
QAXhi <- quantile(A75$axis, .25); QEFF <- quantile(A75$effgap, .60); QARM <- quantile(A75$armdiff, .75)

M[, `:=`(
  A1_locked   = axis < 10 & arm >= arm_thr & is.finite(disc) & disc <= 19.71 & vs >= SEP3,
  A2_wide     = axis < 11 & arm >= TH55 & is.finite(disc) & disc <= 27.5 & vs >= SEP3,
  A3_wide75   = axis < 11 & arm >= TH55 & is.finite(disc) & disc <= 27.5 & vs >= SEP3 & nsw >= 75,
  A4_broad    = nsw >= 75 & axis <= QAX & vs >= QVS,
  A5_seam     = nsw >= 75 & axis <= QAX & vs >= QVS & ssw >= QSSW25,
  A6_extreme  = nsw >= 75 & axis <= QAX & ssw >= QSSW75 & sswsign >= QSGN75,
  A7_trad     = nsw >= 75 & axis <= QAX & vs >= QVS & ssw <= MEDSSW,
  A8_mismatch = nsw >= 75 & season >= 2023 & axis >= QAXhi & effgap <= QEFF & armdiff >= QARM
)]
SPEC <- c(
  A1_locked   = sprintf("axis<10, slot>=pop thr, disc<=19.7, sep>=%.1f  [superseded]", SEP3),
  A2_wide     = sprintf("axis<11, slot>=%.1f, disc<=27.5, sep>=%.1f", TH55, SEP3),
  A3_wide75   = "same as A2 plus a 75-swing floor",
  A4_broad    = sprintf("axis<=%.1f, sep>=%.1f  [simplest surviving rule]", QAX, QVS),
  A5_seam     = sprintf("axis<=%.1f, sep>=%.1f, seam dev>=%.1f  [current best]", QAX, QVS, QSSW25),
  A6_extreme  = sprintf("axis<=%.1f, seam dev>=%.1f, signed>=%.1f", QAX, QSSW75, QSGN75),
  A7_trad     = sprintf("axis<=%.1f, sep>=%.1f, seam dev<=%.1f  [the half that fails]", QAX, QVS, MEDSSW),
  A8_mismatch = sprintf("axis>=%.1f, effgap<=%.2f, armdiff>=%.1f, 2023+  [underperforms]", QAXhi, QEFF, QARM))
AN <- names(SPEC)

crob <- function(D, f, k, w = TRUE) {
  environment(f) <- environment()
  D <- D[complete.cases(D[, c(all.vars(f), "id", "nsw"), with = FALSE])]
  wt <- if (w) as.numeric(D$nsw) else rep(1, nrow(D))
  m <- lm(f, D, weights = wt); u <- residuals(m)*sqrt(wt); X <- model.matrix(m)*sqrt(wt)
  nc <- uniqueN(D$id); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k]); c(e, s, 2*pt(-abs(e/s), nc-1))
}

## =============================================================================================
cat("=== 1. THE EIGHT DEFINITIONS SIDE BY SIDE ===\n\n")
SUM <- rbindlist(lapply(AN, function(a) {
  D <- if (a == "A8_mismatch") M[nsw >= 75 & season >= 2023] else if (grepl("A[4-7]", a)) M[nsw >= 75] else M
  D <- copy(D)[, gg := get(a)]
  r <- crob(D, y ~ gg + vs, "ggTRUE"); rc <- crob(D[is.finite(ych)], ych ~ gg + vs, "ggTRUE")
  rv <- crob(D[is.finite(yrv)], yrv ~ gg + vs, "ggTRUE")
  data.table(archetype = a, definition = SPEC[[a]], seasons = sum(D$gg), arms = uniqueN(D[gg == TRUE]$id),
             swings = D[gg == TRUE, sum(nsw)], whiff = r[1], p_wh = r[3], chase = rc[1], rv = rv[1])
}))
print(SUM[, .(archetype, seasons, arms, swings, whiff = round(whiff,2), p = round(p_wh,4),
              chase = round(chase,2), rv = round(rv,2), definition)], row.names = FALSE)

## =============================================================================================
cat("\n=== 2. ROSTERS ===\n")
for (a in AN) {
  D <- M[get(a) == TRUE][order(-y)]
  cat(sprintf("\n--- %s : %s ---\n  %d seasons, %d arms\n", a, SPEC[[a]], nrow(D), uniqueN(D$id)))
  print(D[, .(pitcher = nm, season, swings = nsw, axis = round(axis,1), sep = round(vs,1),
              slot = round(slot,0), seam_dev = round(ssw,1), signed = round(sswsign,1),
              armdiff = round(armdiff,1), disc = round(disc,1), whiff = round(y,1),
              chase = round(ych,1), rv = round(yrv,2))], row.names = FALSE)
}

## =============================================================================================
cat("\n\n=== 3. SEASONS THAT SATISFY SEVERAL DEFINITIONS AT ONCE ===\n")
GOOD <- c("A1_locked","A2_wide","A3_wide75","A4_broad","A5_seam","A6_extreme")
M[, nhits := rowSums(as.matrix(M[, ..GOOD]))]
X <- M[nhits >= 2][order(-nhits, -y)]
X[, which := apply(as.matrix(X[, ..GOOD]), 1, function(r) paste(sub("^A[0-9]_", "", GOOD[r == 1]), collapse = ", "))]
print(X[, .(pitcher = nm, season, swings = nsw, hits = nhits, whiff = round(y,1), chase = round(ych,1),
            rv = round(yrv,2), definitions = which)], row.names = FALSE)

cat("\n=== 4. ONE FLAGSHIP PER ARCHETYPE ===\n")
# the season with the most swings among those beating the model, so the example is not a small-sample fluke
for (a in AN) {
  s <- M[get(a) == TRUE & ((a == "A8_mismatch" & y < 0) | (a != "A8_mismatch" & y > 0))][order(-nsw)][1]
  if (!nrow(s) || is.na(s$nm)) next
  cat(sprintf("  %-12s %-22s %d  %3d swings | axis %4.1f  sep %4.1f  slot %2.0f  seam %4.1f | whiff %+5.1f  chase %+5.1f  rv %+5.2f\n",
              a, s$nm, s$season, s$nsw, s$axis, s$vs, s$slot, s$ssw, s$y, s$ych, s$yrv))
}
cat("\n=== 5. HOW MUCH DO THE DEFINITIONS OVERLAP? (shared seasons) ===\n")
O <- outer(AN, AN, Vectorize(function(a, b) M[get(a) == TRUE & get(b) == TRUE, .N]))
dimnames(O) <- list(AN, AN); print(O)
saveRDS(M, file.path(MDIR, "archetype_roster.rds"))

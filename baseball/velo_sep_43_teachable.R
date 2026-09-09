#!/usr/bin/env Rscript

# CAN THIS BE TAUGHT, AND IS IT WORTH TEACHING?
#
# Everything so far is cross-sectional: arms WITH the profile outperform arms WITHOUT it. That is
# consistent with the profile causing the gain and equally consistent with good pitchers happening
# to have it. Neither supports "teach this" on its own.
#
# The design that does is within-pitcher. Take one arm across consecutive seasons, watch his
# changeup's shape move, and see whether his run value moves with it. Every fixed thing about him -
# fastball, command, arm, organisation - differences out. Four questions in order:
#
#   1. Is the profile rare, and is it spreading?
#   2. Are the traits actually movable, or is a changeup's shape fixed at birth?
#   3. When an arm moves toward the profile, does his run value follow? (dose-response)
#   4. What happens to arms who cross into the archetype, against arms who never do? (event study)
#
# Outcome is the count-and-handedness-adjusted run value built in the previous step, so none of this
# can be an artifact of when the pitch gets thrown.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 235); MDIR <- "data/statcast_model"

P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(rv) & is.finite(ax) & is.finite(az) & is.finite(sax) & is.finite(cax) &
       is.finite(speed_diff) & is.finite(axis_diff) & is.finite(arm_angle) & is.finite(arm_diff)]
P[, `:=`(id = as.character(pitcher), lh = p_throws == "L")]
P[, `:=`(mx = fifelse(lh, -ax, ax), mz = az + 32.174, sep = -speed_diff)]
P[, spin_axis := (atan2(sax, cax) * 180/pi) %% 360]
P[, spin_axis := fifelse(lh, (360 - spin_axis) %% 360, spin_axis)]
wrap <- function(d) ((d + 180) %% 360) - 180
P[, dev := wrap(wrap(atan2(mx, mz)*180/pi + (spin_axis - 180)) -
                median(wrap(atan2(mx, mz)*180/pi + (spin_axis - 180)), na.rm = TRUE))]
P[, absdev := abs(dev)]

TR <- P[, .(np = .N, nm = nm <- player_name[1],
            sep = mean(sep), axis = mean(axis_diff), seam = mean(absdev), signed = mean(dev),
            slot = mean(arm_angle), velo = mean(release_speed),
            spin = mean(release_spin_rate)), by = .(id, season)]
CA <- readRDS(file.path(MDIR, "count_adjusted_rv.rds"))$season
TR <- merge(TR, CA[, .(id, season, adj, adjs, raw)], by = c("id","season"))
M <- readRDS(file.path(MDIR, "archetype_roster.rds")); setDT(M)
TAGS <- c("A4_broad","A5_seam")
TR <- merge(TR, M[, c("id","season","nsw", TAGS), with = FALSE], by = c("id","season"))
TR[, para := A4_broad | A5_seam]
W <- P[is_swing == 1, .(nsw2 = .N, whiff = 100*mean(whiff == 1)), by = .(id, season)]
TR <- merge(TR, W, by = c("id","season"))
TR <- TR[np >= 120]
cat(sprintf("%d pitcher-seasons with 120+ changeups, %d arms.\n", nrow(TR), uniqueN(TR$id)))

VARS <- c("sep","axis","seam","slot")
NICE <- c(sep = "Separation (mph)", axis = "Spin-axis gap (deg)", seam = "Seam deviation (deg)",
          slot = "Arm slot (deg)")
GOODDIR <- c(sep = 1, axis = -1, seam = 1, slot = 1)   # direction that moves toward the profile

## =============================================================================================
cat("\n=== 1. HOW RARE IS IT, AND IS IT SPREADING? ===\n\n")
PR <- TR[, .(seasons = .N, para = sum(para), share = 100*mean(para),
             sep = mean(sep), seam = mean(seam), axis = mean(axis),
             adj = weighted.mean(adj, np)), by = season][order(season)]
print(PR[, .(season, seasons, `parachute seasons` = para, `share %` = round(share,1),
             `mean separation` = round(sep,2), `mean seam dev` = round(seam,2),
             `mean axis gap` = round(axis,1), `adj rv/100` = round(adj,3))], row.names = FALSE)
for (v in VARS) {
  f <- lm(get(v) ~ season, TR, weights = np); s <- summary(f)$coefficients[2,]
  cat(sprintf("  %-24s league trend %+.3f per season  (p %.4f)\n", NICE[[v]], s[1], s[4]))
}

## =============================================================================================
cat("\n=== 2. ARE THE TRAITS MOVABLE? ===\n\n")
setorder(TR, id, season)
D <- TR[, .SD[order(season)], by = id]
D[, gap_yr := season - shift(season), by = id]
for (v in c(VARS, "adj","adjs","whiff","raw")) D[, (paste0("d_", v)) := get(v) - shift(get(v)), by = id]
# lags must come off the full panel, not the differenced subset, or the first pair of every arm
# loses its predecessor and gets dropped
for (v in c("para","adj","whiff","sep","axis","seam")) D[, (paste0("l_", v)) := shift(get(v)), by = id]
DD <- D[gap_yr == 1 & !is.na(d_sep)]
cat(sprintf("  %d consecutive-season pairs, %d arms.\n\n", nrow(DD), uniqueN(DD$id)))
cat(sprintf("  %-24s %10s %12s %12s %10s %12s\n", "trait", "between SD", "yr-to-yr SD",
            "ratio", "median |move|", "% moving 1+ SD"))
for (v in VARS) {
  bsd <- TR[, sd(get(v))]; wsd <- DD[, sd(get(paste0("d_", v)))]
  cat(sprintf("  %-24s %10.2f %12.2f %12.2f %10.2f %11.0f%%\n", NICE[[v]], bsd, wsd, wsd/bsd,
              DD[, median(abs(get(paste0("d_", v))))],
              100*DD[, mean(abs(get(paste0("d_", v))) > bsd)]))
}
cat("\n  a ratio near 1 means an arm can move the trait about as far in one offseason as the whole\n")
cat("  league spans. a ratio near 0 means the trait is a fixed property of the pitcher.\n")

## =============================================================================================
cat("\n=== 3. WHEN AN ARM MOVES TOWARD THE PROFILE, DOES RUN VALUE FOLLOW? ===\n\n")
clus <- function(dat, f, k, wt = "np") {
  environment(f) <- environment()
  dat <- dat[complete.cases(dat[, c(all.vars(f), "id", wt), with = FALSE])]
  w <- as.numeric(dat[[wt]]); m <- lm(f, dat, weights = w)
  u <- residuals(m)*sqrt(w); X <- model.matrix(m)*sqrt(w); nc <- uniqueN(dat$id)
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, dat$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k]); c(e, s, 2*pt(-abs(e/s), nc-1))
}
st <- function(r, d = 3) if (is.na(r[1])) "      -   " else
  sprintf(paste0("%+", d+5, ".", d, "f%s"), r[1], ifelse(r[3] < .01, "**", ifelse(r[3] < .05, "* ", "  ")))
Z <- copy(DD); for (v in VARS) Z[, (paste0("z_", v)) := scale(get(paste0("d_", v)))[,1]]
cat("  change in outcome regressed on change in trait, one SD of change, season effects included.\n")
cat("  every fixed attribute of the pitcher differences out.\n\n")
cat(sprintf("  %-24s %14s %14s %14s\n", "change in trait", "adj rv/100", "+ stuff", "whiff pts"))
for (v in VARS) {
  cells <- sapply(c("d_adj","d_adjs","d_whiff"), function(y)
    st(clus(Z, as.formula(paste0(y, " ~ z_", v, " + factor(season)")), paste0("z_", v)),
       if (y == "d_whiff") 2 else 3))
  cat(sprintf("  %-24s %14s %14s %14s\n", NICE[[v]], cells[1], cells[2], cells[3]))
}
cat("\n  all four at once, so each is net of the others:\n\n")
f <- as.formula(paste("d_adj ~", paste(paste0("z_", VARS), collapse = " + "), "+ factor(season)"))
for (v in VARS) cat(sprintf("  %-24s %14s\n", NICE[[v]], st(clus(Z, f, paste0("z_", v)))))

cat("\n  a composite move: standardised distance travelled toward the profile\n\n")
Z[, toward := rowMeans(cbind(z_sep*GOODDIR[["sep"]], z_axis*GOODDIR[["axis"]],
                             z_seam*GOODDIR[["seam"]]))]
Z[, tq := cut(toward, quantile(toward, 0:5/5), include.lowest = TRUE,
              labels = c("away, hard","away","flat","toward","toward, hard"))]
print(Z[, .(pairs = .N, `mean move (SD)` = round(mean(toward),2),
            `d sep` = round(mean(d_sep),2), `d axis` = round(mean(d_axis),1),
            `d seam` = round(mean(d_seam),2),
            `d adj rv` = round(weighted.mean(d_adj, np),3),
            `d whiff` = round(weighted.mean(d_whiff, nsw2),2)), by = tq][order(tq)], row.names = FALSE)
r <- clus(Z, d_adj ~ toward + factor(season), "toward")
cat(sprintf("\n  slope on the composite: %s runs per 100 per SD of movement toward the profile\n", st(r)))
rw <- clus(Z, d_whiff ~ toward + factor(season), "toward", "nsw2")
cat(sprintf("  and on whiff:            %s points per SD\n", st(rw, 2)))

cat("\n  --- is this just mean reversion? ---\n\n")
cat("  the worry: arms who got hit hard last year both tinker more AND regress upward anyway,\n")
cat("  which would manufacture a slope with no causal content. two checks.\n\n")
rp <- clus(Z, toward ~ l_adj + factor(season), "l_adj")
cat(sprintf("  (a) does last year's run value predict how much an arm moves?  %s SD per run/100\n", st(rp)))
cat(sprintf("      if that is flat, selection into moving is not outcome-driven.\n\n"))
cat("  (b) the slope with last year's run value controlled, which absorbs reversion directly:\n\n")
for (lab in c("no control","+ prior rv","+ prior rv and prior shape")) {
  f <- switch(lab,
    "no control" = d_adj ~ toward + factor(season),
    "+ prior rv" = d_adj ~ toward + l_adj + factor(season),
    d_adj ~ toward + l_adj + l_sep + l_axis + l_seam + factor(season))
  fw <- switch(lab,
    "no control" = d_whiff ~ toward + factor(season),
    "+ prior rv" = d_whiff ~ toward + l_whiff + factor(season),
    d_whiff ~ toward + l_whiff + l_sep + l_axis + l_seam + factor(season))
  cat(sprintf("  %-28s adj rv %s   whiff %s\n", lab, st(clus(Z, f, "toward")),
              st(clus(Z, fw, "toward", "nsw2"), 2)))
}

## =============================================================================================
cat("\n=== 3b. BETWEEN ARMS VS WITHIN ARMS, SAME SCALE ===\n\n")
# first differences throw away every arm-season that has no neighbour and double the measurement
# noise. the two-way fixed effect estimator uses all 577 seasons and every within-arm contrast, so
# it is the specification with the most power. reported next to the cross-sectional estimate the
# whole argument has rested on.
ci <- function(dat, f, k, wt) {
  environment(f) <- environment()
  dat <- dat[complete.cases(dat[, c(all.vars(f), "id", wt), with = FALSE])]
  w <- as.numeric(dat[[wt]]); m <- lm(f, dat, weights = w)
  u <- residuals(m)*sqrt(w); X <- model.matrix(m)*sqrt(w); nc <- uniqueN(dat$id)
  keep <- !is.na(coef(m)); X <- X[, keep, drop = FALSE]
  b <- tryCatch(solve(crossprod(X)), error = function(e) MASS::ginv(crossprod(X)))
  V <- b %*% crossprod(rowsum(X*u, dat$id)) %*% b * (nc/(nc-1))
  j <- which(names(coef(m))[keep] == k)
  e <- unname(coef(m)[keep][j]); s <- unname(sqrt(diag(V))[j]); t <- qt(.975, nc-1)
  c(est = e, se = s, lo = e - t*s, hi = e + t*s, p = 2*pt(-abs(e/s), nc-1))
}
TR2 <- copy(TR)
for (v in VARS) TR2[, (paste0("z_", v)) := scale(get(v))[,1] * GOODDIR[[v]]]
TR2[, toward := rowMeans(cbind(z_sep, z_axis, z_seam))]
sdw <- TR2[, sd(toward)]
cat("  effect of a one-SD move toward the profile, 95% intervals, arm-clustered.\n")
cat("  'toward' combines separation, axis gap and seam deviation.\n\n")
cat(sprintf("  %-34s %9s %9s %20s %8s\n", "outcome / specification", "estimate", "std err",
            "95% interval", "n"))
SPECS <- list(
  list("adj rv/100", "between arms (cross-section)", TR2, adj ~ toward, "np"),
  list("adj rv/100", "within arms (arm + season FE)", TR2, adj ~ toward + factor(id) + factor(season), "np"),
  list("whiff pts",  "between arms (cross-section)", TR2, whiff ~ toward, "nsw2"),
  list("whiff pts",  "within arms (arm + season FE)", TR2, whiff ~ toward + factor(id) + factor(season), "nsw2"))
for (s in SPECS) {
  r <- ci(s[[3]], s[[4]], "toward", s[[5]])
  cat(sprintf("  %-12s %-21s %+9.3f %9.3f   [%+7.3f, %+7.3f] %6d%s\n", s[[1]], s[[2]],
              r["est"], r["se"], r["lo"], r["hi"], nrow(s[[3]]),
              ifelse(r["p"] < .01, " **", ifelse(r["p"] < .05, " * ", "   "))))
}
cat("\n  the same for first differences, for comparison:\n\n")
for (o in c("d_adj","d_whiff")) {
  r <- ci(Z, as.formula(paste(o, "~ toward + l_adj + factor(season)")), "toward",
          if (o == "d_whiff") "nsw2" else "np")
  cat(sprintf("  %-12s %-21s %+9.3f %9.3f   [%+7.3f, %+7.3f] %6d\n",
              ifelse(o == "d_adj", "adj rv/100", "whiff pts"), "first differences",
              r["est"], r["se"], r["lo"], r["hi"], nrow(Z)))
}

cat("\n  --- how much power is there, really? ---\n\n")
# split-half reliability of each season-level outcome tells us how much of the variance the design
# is even allowed to explain. an outcome that does not repeat cannot be predicted by anything.
rel <- function(col, nwt) {
  X <- TR2[get(nwt) >= 60]
  s <- X[, .(a = get(col), n = get(nwt))]
  1 - (mean(s$a^2, na.rm = TRUE) - mean(s$a, na.rm = TRUE)^2)^0 * 0  # placeholder, see YoY below
}
YY <- merge(TR[, .(id, season, adj, whiff, np, nsw2)],
            TR[, .(id, season = season - 1, adj1 = adj, whiff1 = whiff)], by = c("id","season"))
cat(sprintf("  year-over-year correlation, adj rv   r = %+.3f  (n = %d)\n",
            YY[, cor(adj, adj1)], nrow(YY)))
cat(sprintf("  year-over-year correlation, whiff    r = %+.3f  (n = %d)\n\n",
            YY[, cor(whiff, whiff1)], nrow(YY)))
mde <- function(spec) { r <- ci(spec[[3]], spec[[4]], "toward", spec[[5]]); 2.8*r["se"] }
cat(sprintf("  smallest effect the within-arm design could detect at 80%% power:\n"))
cat(sprintf("    adj rv/100  %.3f runs per SD\n", mde(SPECS[[2]])))
cat(sprintf("    whiff       %.3f points per SD\n", mde(SPECS[[4]])))

## =============================================================================================
cat("\n=== 4. EVENT STUDY: ARMS WHO CROSS INTO THE ARCHETYPE ===\n\n")
DD[, ev := fifelse(para & !l_para, "entered",
           fifelse(!para & l_para, "left",
           fifelse(para & l_para, "stayed in", "stayed out")))]
E <- DD[, .(pairs = .N, arms = uniqueN(id),
            d_sep = mean(d_sep), d_axis = mean(d_axis), d_seam = mean(d_seam),
            before = weighted.mean(adj - d_adj, np), after = weighted.mean(adj, np),
            d_adj = weighted.mean(d_adj, np), d_whiff = weighted.mean(d_whiff, nsw2)), by = ev]
print(E[order(-pairs), .(event = ev, pairs, arms, `d sep` = round(d_sep,2),
                         `d axis` = round(d_axis,1), `d seam` = round(d_seam,2),
                         `rv before` = round(before,3), `rv after` = round(after,3),
                         `change` = round(d_adj,3), `d whiff` = round(d_whiff,2))], row.names = FALSE)
cat("\n  difference in differences against the arms who stayed out:\n\n")
for (e in c("entered","left","stayed in")) {
  X <- DD[ev %in% c(e, "stayed out")][, g := ev == e]
  a <- clus(X, d_adj ~ g, "gTRUE"); b <- clus(X, d_whiff ~ g, "gTRUE", "nsw2")
  cat(sprintf("  %-12s vs stayed out   adj rv %s   whiff %s   (n = %d)\n",
              e, st(a), st(b, 2), nrow(X[g == TRUE])))
}

## =============================================================================================
cat("\n=== 5. HOW MANY ARMS ARE WITHIN REACH? ===\n\n")
QAX <- quantile(M[nsw >= 75]$axis, .25); QVS <- quantile(M[nsw >= 75]$vs, .60)
cat(sprintf("  the broad gate is axis gap <= %.1f deg and separation >= %.1f mph.\n\n", QAX, QVS))
L <- TR[season == max(season) - 1 | season == max(season)]
L <- L[, .SD[which.max(season)], by = id]
L[, `:=`(need_ax = pmax(axis - QAX, 0), need_sep = pmax(QVS - sep, 0))]
L[, dist := sqrt((need_ax/DD[, sd(d_axis)])^2 + (need_sep/DD[, sd(d_sep)])^2)]
cat(sprintf("  %d arms in the most recent season.\n", nrow(L)))
cat(sprintf("    already inside the gate                       %3d (%.0f%%)\n",
            L[need_ax == 0 & need_sep == 0, .N], 100*L[, mean(need_ax == 0 & need_sep == 0)]))
cat(sprintf("    within one typical offseason move of it       %3d (%.0f%%)\n",
            L[dist > 0 & dist <= 1, .N], 100*L[, mean(dist > 0 & dist <= 1)]))
cat(sprintf("    within two                                    %3d (%.0f%%)\n",
            L[dist > 1 & dist <= 2, .N], 100*L[, mean(dist > 1 & dist <= 2)]))
cat(sprintf("    further than that                             %3d (%.0f%%)\n",
            L[dist > 2, .N], 100*L[, mean(dist > 2)]))
cat("\n  the twenty closest arms currently outside the gate:\n\n")
print(L[dist > 0][order(dist)][1:20, .(pitcher = nm, season, changeups = np,
        sep = round(sep,1), `sep needed` = round(need_sep,1), axis = round(axis,1),
        `axis needed` = round(need_ax,1), seam = round(seam,1),
        `adj rv/100` = round(adj,2), moves = round(dist,2))], row.names = FALSE)
saveRDS(list(pairs = DD, latest = L, trend = PR), file.path(MDIR, "teachable.rds"))

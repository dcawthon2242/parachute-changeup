#!/usr/bin/env Rscript

# DO THE WHIFF GAINS TURN INTO RUNS, AND WHAT DO THEY COST?
#
# The archetypes beat the whiff model by 4 to 13 points. A whiff is worth about 0.15 runs more than
# a ball in play and 0.08 more than a foul, so those gains imply a specific run value surplus. The
# question is whether the surplus shows up, and if it falls short, which channel is leaking.
#
# Every pitch lands in exactly one of five channels and the five run value contributions add back to
# the season's run value per 100 pitches:
#     ball (taken outside the zone), called strike (taken inside), whiff, foul, ball in play.
# Regressing each channel on the archetype flag while holding the stuff model's predicted run value
# fixed says where an archetype's runs come from and where they go.
#
# Drawbacks tested: zone avoidance and the walk bill, damage on the contact they do allow, home runs,
# the platoon split, how much of the effect needs two strikes, and whether the run value edge is as
# repeatable as the whiff edge.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 235); MDIR <- "data/statcast_model"

P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(rv)][, id := as.character(pitcher)]
P[, inzone := (z_rel_bot >= 0 & z_rel_top <= 0 & abs(plate_x) <= 0.83)]
P[, ch := fifelse(is_swing == 1 & whiff == 1, "whiff",
          fifelse(is_swing == 1 & is_bip == 1, "bip",
          fifelse(is_swing == 1, "foul",
          fifelse(inzone, "cstrike", "ball"))))]
CH <- c("ball","cstrike","whiff","foul","bip")
CN <- c(ball = "Ball", cstrike = "Called strike", whiff = "Whiff", foul = "Foul", bip = "Ball in play")

cat("=== THE EXCHANGE RATE ===\n\n")
EX <- P[, .(n = .N, share = 100*.N/nrow(P), rv = mean(rv)), by = ch][order(-n)]
print(EX[, .(channel = CN[ch], pitches = n, `share %` = round(share,1), `mean run value` = round(rv,4))],
      row.names = FALSE)
sw <- P[is_swing == 1]; nsw <- nrow(sw)
pB <- sw[ch == "bip", .N]/(nsw - sw[ch == "whiff", .N]); pF <- 1 - pB
gain <- pB*(EX[ch=="whiff", rv] - EX[ch=="bip", rv]) + pF*(EX[ch=="whiff", rv] - EX[ch=="foul", rv])
swrate <- nsw/nrow(P)
PERPT <- gain * swrate            # runs per 100 pitches for one extra whiff point
cat(sprintf("\n  swings are %.1f%% of pitches. a swing that is not a whiff is %.0f%% a ball in play,\n",
            100*swrate, 100*pB))
cat(sprintf("  %.0f%% a foul. converting one of them to a whiff is worth %+.4f runs.\n", 100*pF, gain))
cat(sprintf("  SO: one extra whiff point should be worth %+.3f runs per 100 pitches.\n", 100*PERPT/100))

## ---- season-level panel -----------------------------------------------------------------------
S <- P[, .(np = .N, rv100 = 100*mean(rv), nsw = sum(is_swing == 1),
           whiff = 100*sum(whiff == 1)/sum(is_swing == 1),
           zone = 100*mean(inzone), swing = 100*mean(is_swing == 1),
           chase = 100*sum(is_swing == 1 & !inzone)/sum(!inzone),
           zsw = 100*sum(is_swing == 1 & inzone)/sum(inzone),
           hr100bip = 100*sum(launch_speed >= 100 & !is.na(bb_type) & bb_type == "fly_ball", na.rm = TRUE)/
                      sum(is_bip == 1),
           xwcon = 1000*mean(estimated_woba_using_speedangle, na.rm = TRUE),
           sep = mean(-speed_diff)), by = .(id, season)]
for (c in CH) {
  A <- P[ch == c, .(r = 100*sum(rv), k = .N), by = .(id, season)]
  S <- merge(S, A, by = c("id","season"), all.x = TRUE)
  S[is.na(r), r := 0][is.na(k), k := 0]
  setnames(S, c("r","k"), paste0(c("rv_","n_"), c))
}
for (c in CH) { S[, (paste0("rv_", c)) := get(paste0("rv_", c))/np]
                S[, (paste0("rt_", c)) := 100*get(paste0("n_", c))/np] }

OOF <- readRDS(file.path(MDIR, "absorb_preds.rds"))$STUFF
PR <- OOF$r[, .(pred_rv = 100*mean(p), act_rv = 100*mean(y)), by = .(id, season)]
PW <- OOF$w[, .(pred_wh = 100*mean(p), act_wh = 100*mean(y), nw = .N), by = .(id, season)]
S <- merge(merge(S, PR, by = c("id","season")), PW, by = c("id","season"))
S[, `:=`(r_rv = act_rv - pred_rv, r_wh = act_wh - pred_wh)]

M <- readRDS(file.path(MDIR, "archetype_roster.rds")); setDT(M)
TAGS <- c("A2_wide","A4_broad","A5_seam","A6_extreme","A7_trad","A8_mismatch")
LAB <- c(A2_wide = "Matched axis, wide", A4_broad = "Broad (axis + separation)",
         A5_seam = "Seam-shifted matched", A6_extreme = "Extreme seam shift",
         A7_trad = "Traditional, low seam dev", A8_mismatch = "Mismatch (underperforms)")
S <- merge(S, M[, c("id","season","nm", TAGS), with = FALSE], by = c("id","season"))
D <- S[nw >= 75]
cat(sprintf("\n%d pitcher-seasons with 75+ swings, %d arms.\n", nrow(D), uniqueN(D$id)))

crob <- function(dat, f, k, wt = "np") {
  environment(f) <- environment()
  dat <- dat[complete.cases(dat[, c(all.vars(f), "id", wt), with = FALSE])]
  w <- as.numeric(dat[[wt]]); m <- lm(f, dat, weights = w)
  u <- residuals(m)*sqrt(w); X <- model.matrix(m)*sqrt(w)
  nc <- uniqueN(dat$id); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X*u, dat$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k]); c(e, s, 2*pt(-abs(e/s), nc-1))
}
st <- function(r, d = 2) if (is.na(r[1])) "    -  " else
  sprintf(paste0("%+", d+5, ".", d, "f%s"), r[1], ifelse(r[3] < .01, "**", ifelse(r[3] < .05, "* ", "  ")))

## =============================================================================================
cat("\n=== 1. DOES THE WHIFF SURPLUS ARRIVE AS RUNS? ===\n\n")
cat("  whiff surplus in points, the run value it implies at the exchange rate above, and the run\n")
cat("  value surplus actually observed. all against the rest of the population, arm-clustered.\n\n")
cat(sprintf("  %-28s %8s %10s %12s %12s %10s\n", "archetype", "seasons", "whiff", "implied rv",
            "actual rv", "delivered"))
CONV <- rbindlist(lapply(TAGS, function(t) {
  X <- copy(D)[, g := get(t)]
  w <- crob(X, r_wh ~ g, "gTRUE", "nw"); v <- crob(X, r_rv ~ g, "gTRUE", "np")
  imp <- w[1] * PERPT
  data.table(tag = t, seasons = sum(X$g), wh = w[1], pw = w[3], imp = imp, rv = v[1], pv = v[3],
             del = 100*v[1]/imp)
}))
for (i in seq_len(nrow(CONV))) with(CONV[i], cat(sprintf(
  "  %-28s %8d %10s %12s %12s %9.0f%%\n", LAB[[tag]], seasons,
  st(c(wh, NA, pw)), sprintf("%+.3f", imp), st(c(rv, NA, pv), 3), del)))
cat("\n  the same on the continuous gradient, per 1 SD:\n\n")
for (v in c("sep","whiff")) {
  X <- copy(D)[, z := scale(get(v))[,1]]
  w <- crob(X, r_wh ~ z, "z", "nw"); r <- crob(X, r_rv ~ z, "z", "np")
  cat(sprintf("    %-22s whiff %s   implied %+.3f   actual rv %s   delivered %3.0f%%\n",
              v, st(w), w[1]*PERPT, st(r, 3), 100*r[1]/(w[1]*PERPT)))
}

## =============================================================================================
cat("\n=== 2. WHICH CHANNEL DOES THE RUN VALUE COME FROM, AND WHICH ONE LEAKS? ===\n\n")
cat("  run value per 100 pitches booked in each channel, archetype minus the rest, holding the\n")
cat("  stuff model's predicted run value fixed. the five columns add to the total.\n\n")
cat(sprintf("  %-28s %s %10s\n", "archetype",
            paste(sprintf("%12s", CN[CH]), collapse = ""), "total"))
for (t in TAGS) {
  X <- copy(D)[, g := get(t)]
  cells <- sapply(CH, function(c) st(crob(X, as.formula(paste0("rv_", c, " ~ g + pred_rv")), "gTRUE"), 3))
  tot <- st(crob(X, rv100 ~ g + pred_rv, "gTRUE"), 3)
  cat(sprintf("  %-28s %s %10s\n", LAB[[t]], paste(sprintf("%12s", cells), collapse = ""), tot))
}
cat("\n  and the rate of each channel, per 100 pitches (same controls):\n\n")
cat(sprintf("  %-28s %s\n", "archetype", paste(sprintf("%12s", CN[CH]), collapse = "")))
for (t in TAGS) {
  X <- copy(D)[, g := get(t)]
  cells <- sapply(CH, function(c) st(crob(X, as.formula(paste0("rt_", c, " ~ g + pred_rv")), "gTRUE")))
  cat(sprintf("  %-28s %s\n", LAB[[t]], paste(sprintf("%12s", cells), collapse = "")))
}

## =============================================================================================
cat("\n=== 3. THE DRAWBACKS ===\n\n")
cat("  (a) command and usage, weighted by pitches\n\n")
COST <- c(zone = "Zone rate %", chase = "Chase rate %", zsw = "In-zone swing %",
          swing = "Swing rate %", np = "Changeups thrown")
cat(sprintf("  %-28s %s\n", "archetype", paste(sprintf("%16s", COST), collapse = "")))
for (t in TAGS) {
  X <- copy(D)[, g := get(t)]
  cells <- sapply(names(COST), function(c)
    st(crob(X, as.formula(paste0(c, " ~ g + pred_rv")), "gTRUE"), if (c == "np") 0 else 2))
  cat(sprintf("  %-28s %s\n", LAB[[t]], paste(sprintf("%16s", cells), collapse = "")))
}

# contact has to be weighted by balls in play, not by pitches. weighting a per-ball-in-play average
# by the season's total changeup count overweights arms who threw many and put few in play, which
# was inflating both the effect sizes and their significance.
cat("\n  (b) the contact they do allow, weighted by BALLS IN PLAY\n\n")
BP <- P[is_bip == 1 & !is.na(bb_type) & bb_type != "",
        .(nbip = .N, xw = 1000*mean(estimated_woba_using_speedangle, na.rm = TRUE),
          hr = 100*mean(rv < -1.0), gbr = 100*mean(bb_type == "ground_ball"),
          rvbip = 100*mean(rv)), by = .(id, season)]
DB <- merge(D, BP, by = c("id","season"))[nbip >= 30]
CQ <- c(xw = "xwOBA on contact", hr = "Home runs per 100 BIP", gbr = "Ground ball %",
        rvbip = "Run value per 100 BIP")
cat(sprintf("  %-28s %8s %s\n", "archetype", "seasons", paste(sprintf("%22s", CQ), collapse = "")))
for (t in TAGS) {
  X <- copy(DB)[, g := get(t)]
  cells <- sapply(names(CQ), function(c)
    st(crob(X, as.formula(paste0(c, " ~ g + pred_rv")), "gTRUE", "nbip"), 1))
  cat(sprintf("  %-28s %8d %s\n", LAB[[t]], sum(X$g), paste(sprintf("%22s", cells), collapse = "")))
}

## =============================================================================================
cat("\n=== 4. PLATOON: DOES IT WORK AGAINST SAME-HANDED HITTERS? ===\n\n")
PL <- rbindlist(lapply(c(0,1), function(h) {
  Q <- P[same_hand == h]
  A <- Q[, .(np = .N, nw = sum(is_swing == 1), rv100 = 100*mean(rv),
             whiff = 100*sum(whiff == 1)/sum(is_swing == 1)), by = .(id, season)]
  A[, same := h]
}))
OW <- OOF$w; OW[, k := .I]
PWH <- merge(P[is_swing == 1 & is.finite(whiff), .(id, season, same_hand)][, k := .I],
             OW[, .(k, p)], by = "k")
PB <- PWH[, .(pw = 100*mean(p)), by = .(id, season, same = same_hand)]
PR2 <- merge(P[, .(id, season, same_hand)][, k := .I], OOF$r[, .(k = .I, p)], by = "k")
PB2 <- PR2[, .(pr = 100*mean(p)), by = .(id, season, same = same_hand)]
PL <- merge(merge(PL, PB, by = c("id","season","same")), PB2, by = c("id","season","same"))
PL <- merge(PL, M[, c("id","season", TAGS), with = FALSE], by = c("id","season"))
PL[, `:=`(r_wh = whiff - pw, r_rv = rv100 - pr)]
cat(sprintf("  %-28s %26s %26s\n", "", "OPPOSITE-handed hitters", "SAME-handed hitters"))
cat(sprintf("  %-28s %8s %8s %8s %8s %8s %8s\n", "archetype", "seasons", "whiff+", "rv+",
            "seasons", "whiff+", "rv+"))
for (t in c("A2_wide","A4_broad","A5_seam","A6_extreme")) {
  out <- lapply(c(0,1), function(h) {
    X <- PL[same == h & nw >= 30][, g := get(t)]
    c(nrow(X[g == TRUE]), st(crob(X, r_wh ~ g, "gTRUE", "nw")), st(crob(X, r_rv ~ g, "gTRUE", "np"), 3))
  })
  cat(sprintf("  %-28s %8s %8s %8s %8s %8s %8s\n", LAB[[t]],
              out[[1]][1], out[[1]][2], out[[1]][3], out[[2]][1], out[[2]][2], out[[2]][3]))
}

## =============================================================================================
cat("\n=== 5. IS THE RUN VALUE EDGE AS REPEATABLE AS THE WHIFF EDGE? ===\n\n")
for (v in c("r_wh","r_rv")) {
  X <- D[, .(id, season, r = get(v), w = if (v == "r_wh") nw else np)]
  Y <- copy(X)[, season := season - 1][, .(id, season, r2 = r)]
  J <- merge(X, Y, by = c("id","season")); ct <- cor.test(J$r, J$r2)
  cat(sprintf("  %-22s %d pairs   sd %6.2f   year to year %+.3f   p %.4f\n",
              if (v == "r_wh") "whiff residual" else "run value residual",
              nrow(J), sd(X$r), ct$estimate, ct$p.value))
}
cat("\n  split-half within season is not available here, so the year-over-year figure is the test.\n")

## =============================================================================================
cat("\n=== 6. HOW BIG IS THIS IN RUNS PER SEASON? ===\n\n")
cat(sprintf("  %-28s %8s %10s %12s %14s\n", "archetype", "seasons", "CH/season", "rv/100", "runs/season"))
for (t in TAGS) {
  X <- D[get(t) == TRUE]; v <- crob(copy(D)[, g := get(t)], r_rv ~ g, "gTRUE", "np")
  cat(sprintf("  %-28s %8d %10.0f %12s %14s\n", LAB[[t]], nrow(X), mean(X$np),
              st(v, 3), sprintf("%+.2f", v[1]*mean(X$np)/100)))
}
saveRDS(list(season = D, platoon = PL, conv = CONV, perpt = PERPT),
        file.path(MDIR, "rv_conversion.rds"))

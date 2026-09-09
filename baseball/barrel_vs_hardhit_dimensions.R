#!/usr/bin/env Rscript

# Which swing dimension drives barrel% vs hard-hit%, and can pitchers be split into
# traditional contact-manager types?
#
# The two rates are nested: hard-hit is an EXIT VELOCITY condition (>=95 mph); barrel is
# hard-hit PLUS a launch-angle condition. So barrel% factors exactly:
#     barrel% = hardhit% x conversion,   conversion = P(barrel | hard-hit)
# Section 2 decomposes exit velocity and launch angle separately against the swing
# dimensions, which is what actually answers "which dimension contributes to which".
# Section 4 splits pitchers on the two factors, and Section 5 asks the question that
# decides whether those types are real: are they repeatable?

suppressPackageStartupMessages({ library(data.table); library(splines) })
set.seed(5)

cols <- c("game_year","game_type","player_name","pitcher","batter","stand","pitch_type",
          "description","bb_type","balls","strikes","plate_x","plate_z","sz_top","sz_bot",
          "release_speed","launch_speed","launch_angle","launch_speed_angle","delta_run_exp",
          "estimated_woba_using_speedangle","bat_speed","swing_length","attack_angle",
          "attack_direction","swing_path_tilt",
          "intercept_ball_minus_batter_pos_x_inches",
          "intercept_ball_minus_batter_pos_y_inches")

dt <- rbindlist(lapply(c(2025, 2026), function(yr)
  fread(file.path("data", sprintf("statcast_%d", yr), sprintf("statcast_%d_all.csv", yr)),
        showProgress = FALSE, select = cols)))
setnames(dt, c("intercept_ball_minus_batter_pos_y_inches",
               "intercept_ball_minus_batter_pos_x_inches"), c("depth","latx"))
dt <- dt[game_type == "R" & pitch_type != "" & balls <= 3 & strikes <= 2]
dt[, pit_rv := -delta_run_exp]

FB<-c("FF","SI","FC"); BR<-c("SL","ST","CU","KC","SV","CS"); OS<-c("CH","FS","FO")
dt[, pgrp := fifelse(pitch_type %in% FB,"FB", fifelse(pitch_type %in% BR,"BR",
             fifelse(pitch_type %in% OS,"OS",NA_character_)))]
dt[, px_bat := fifelse(stand=="R", -plate_x, plate_x)]
dt[, pz_rel := (plate_z - sz_bot)/pmax(sz_top - sz_bot, 0.1)]

# ---- balls in play with the full dimension set ------------------------------
bip <- dt[description == "hit_into_play" & bb_type != "" & !is.na(pgrp) &
          is.finite(launch_speed) & is.finite(launch_angle) & is.finite(depth) &
          is.finite(bat_speed) & is.finite(attack_angle) & is.finite(swing_path_tilt) &
          is.finite(attack_direction) & is.finite(swing_length) & is.finite(latx) &
          is.finite(px_bat) & is.finite(pz_rel) & is.finite(release_speed)]
bip <- bip[is.finite(launch_speed_angle)]
bip[, hard   := launch_speed >= 95]
bip[, barrel := launch_speed_angle == 6]
bip[, sweet  := launch_angle >= 8 & launch_angle <= 32]

cat("############ 1. Sample and the arithmetic of the two rates ############\n")
cat(sprintf("Balls in play with all swing dimensions: %s (2025-2026)\n",
            format(nrow(bip), big.mark=",")))
cat(sprintf("  hard-hit %%  = %.1f\n  barrel %%    = %.1f\n  conversion  = %.1f%% of hard-hit balls are barrels\n",
  100*mean(bip$hard), 100*mean(bip$barrel), 100*sum(bip$barrel)/sum(bip$hard)))
cat(sprintf("  barrels that are not hard-hit: %d (barrel is a strict subset: %s)\n",
  sum(bip$barrel & !bip$hard), all(bip$barrel <= bip$hard)))

# residualise depth the same way as the rest of this project, then demean per hitter
fitd <- lm(depth ~ ns(px_bat,5)*pgrp + ns(pz_rel,5) + ns(release_speed,4) + stand, data=bip)
bip[, r1 := residuals(fitd)]
bip[, nb := .N, by=batter]; bip <- bip[nb >= 50]
bip[, tdev := r1 - mean(r1), by=batter]

DIMS <- c("bat_speed","attack_angle","swing_path_tilt","attack_direction","swing_length",
          "tdev","latx")
LABEL <- c(bat_speed="Bat speed", attack_angle="Attack angle",
           swing_path_tilt="Swing path tilt", attack_direction="Attack direction",
           swing_length="Swing length", tdev="Timing deviation (depth)",
           latx="Lateral intercept")

# =============================================================================
# 2. Which dimension drives exit velocity, and which drives launch angle?
# =============================================================================
cat("\n############ 2. Exit velocity vs launch angle: separate drivers ############\n")
cat("Standardized betas from a joint linear model on all seven dimensions, and each\n")
cat("dimension's UNIQUE R2 contribution (drop-one). Balls in play only.\n\n")
uniq_r2 <- function(form_full, d, dims) {
  full <- summary(lm(form_full, data=d))$r.squared
  sapply(dims, function(v) {
    f <- as.formula(paste(deparse(form_full[[2]]), "~",
      paste(sprintf("scale(%s)", setdiff(dims, v)), collapse=" + ")))
    full - summary(lm(f, data=d))$r.squared
  })
}
rhs <- paste(sprintf("scale(%s)", DIMS), collapse=" + ")
report <- function(yv, lab) {
  f <- as.formula(paste0("scale(", yv, ") ~ ", rhs))
  m <- lm(f, data=bip); co <- summary(m)$coefficients
  u <- uniq_r2(f, bip, DIMS)
  cat(sprintf("  -- %s (total R2 = %.3f) --\n", lab, summary(m)$r.squared))
  o <- data.table(dim=LABEL[DIMS], beta=round(co[-1,1],3), p=signif(co[-1,4],2),
                  uniqR2=round(u,4))[order(-abs(beta))]
  print(o, row.names=FALSE)
  cat("\n")
}
report("launch_speed", "EXIT VELOCITY (the hard-hit channel)")
report("launch_angle", "LAUNCH ANGLE (the extra condition barrel imposes)")

cat("  -- Conditional: among HARD-HIT balls only, what turns one into a barrel? --\n")
hh <- bip[hard == TRUE]
fb <- as.formula(paste0("barrel ~ ", rhs))
mb <- glm(fb, data=hh, family=binomial())
cb <- summary(mb)$coefficients
ob <- data.table(dim=LABEL[DIMS], logodds=round(cb[-1,1],3), p=signif(cb[-1,4],2))[order(-abs(logodds))]
print(ob, row.names=FALSE)
cat(sprintf("  n = %s hard-hit balls, %.1f%% barrels\n", format(nrow(hh), big.mark=","),
            100*mean(hh$barrel)))

cat("\n  -- Same seven dimensions predicting each RATE directly (pseudo-R2 comparison) --\n")
for (yv in c("hard","barrel")) {
  m0 <- glm(as.formula(paste(yv, "~ 1")), data=bip, family=binomial())
  m1 <- glm(as.formula(paste0(yv, " ~ ", rhs)), data=bip, family=binomial())
  cat(sprintf("     %-7s McFadden pseudo-R2 = %.3f\n", yv,
              1 - as.numeric(logLik(m1))/as.numeric(logLik(m0))))
}

# =============================================================================
# 3. What a pitcher can actually move
# =============================================================================
cat("\n############ 3. How much of each dimension does the PITCHER control? ############\n")
cat("Share of variance in each dimension attributable to pitcher identity vs batter identity\n")
cat("(one-way R2 from a factor model on balls in play).\n\n")
cat("Reported as eta-squared, and also net of the group-size bias that inflates it when\n")
cat("units have few observations (expected eta2 under pure noise = (k-1)/(n-1)).\n\n")
eta2 <- function(v, key, d) {
  x <- d[[v]]; g <- d[[key]]
  ok <- is.finite(x); x <- x[ok]; g <- g[ok]
  gm <- mean(x)
  tb <- data.table(x = x, g = g)[, .(n = .N, m = mean(x)), by = g]
  ssb <- sum(tb$n * (tb$m - gm)^2); sst <- sum((x - gm)^2)
  k <- nrow(tb); n <- length(x)
  raw <- ssb/sst
  list(raw = raw, adj = max(0, (raw - (k-1)/(n-1)) / (1 - (k-1)/(n-1))), k = k)
}
cat("Contact depth uses the location/velo-adjusted residual BEFORE hitter demeaning, since\n")
cat("tdev subtracts each hitter's mean and would force his eta2 to zero.\n\n")
EDIMS <- c("bat_speed","attack_angle","swing_path_tilt","attack_direction","swing_length",
           "r1","latx","launch_speed","launch_angle")
ELAB <- c(LABEL, r1="Contact depth (adjusted)",
          launch_speed="EXIT VELOCITY", launch_angle="LAUNCH ANGLE")
vs <- rbindlist(lapply(EDIMS, function(v) {
  rp <- eta2(v, "pitcher", bip); rb <- eta2(v, "batter", bip)
  data.table(dim = ELAB[[v]], pitcher_eta2 = round(rp$adj, 4),
             batter_eta2 = round(rb$adj, 4),
             batter_over_pitcher = round(rb$adj/pmax(rp$adj, 1e-6), 1))
}))
print(vs[order(-pitcher_eta2)], row.names=FALSE)

# =============================================================================
# 4. Splitting pitchers: hard-hit suppression vs barrel conversion
# =============================================================================
cat("\n############ 4. Contact-manager types ############\n")
agg <- function(d, by, minb) {
  a <- d[, .(
    bip = .N,
    hardhit = 100*mean(hard),
    barrel  = 100*mean(barrel),
    conv    = 100*sum(barrel)/pmax(sum(hard),1),
    ev      = mean(launch_speed),
    la      = mean(launch_angle),
    sweet   = 100*mean(sweet),
    gb      = 100*mean(bb_type=="ground_ball"),
    fb      = 100*mean(bb_type=="fly_ball"),
    pop     = 100*mean(bb_type=="popup"),
    batspd  = mean(bat_speed),
    atkang  = mean(attack_angle),
    tdev    = mean(tdev),
    xw      = mean(estimated_woba_using_speedangle, na.rm=TRUE),
    rv      = mean(pit_rv, na.rm=TRUE)
  ), by=by]
  a[bip >= minb]
}
P <- agg(bip[game_year==2026], "pitcher", 120)
NM <- unique(bip[, .(pitcher, player_name)], by="pitcher")
P <- merge(P, NM, by="pitcher")
cat(sprintf("2026 pitchers with >=120 balls in play: %d\n", nrow(P)))

cat("\n  Variance decomposition of barrel%% between pitchers (log scale, barrel = hardhit x conv):\n")
Q <- P[barrel > 0 & conv > 0]
lb <- log(Q$barrel); lh <- log(Q$hardhit); lc <- log(Q$conv)
cat(sprintf("     var(log barrel%%)      = %.4f\n", var(lb)))
cat(sprintf("     var(log hardhit%%)     = %.4f  -> %.0f%% of barrel variance\n",
  var(lh), 100*var(lh)/var(lb)))
cat(sprintf("     var(log conversion)   = %.4f  -> %.0f%% of barrel variance\n",
  var(lc), 100*var(lc)/var(lb)))
cat(sprintf("     2*cov                 = %.4f  -> %.0f%%\n",
  2*cov(lh,lc), 100*2*cov(lh,lc)/var(lb)))
cat(sprintf("     cor(log hardhit, log conversion) = %+.3f\n", cor(lh, lc)))

hmed <- median(P$hardhit); cmed <- median(P$conv)
P[, type := fifelse(hardhit <= hmed & conv <= cmed, "Weak-contact (suppresses both)",
            fifelse(hardhit >  hmed & conv <= cmed, "Angle manager (allows EV, kills LA)",
            fifelse(hardhit <= hmed & conv >  cmed, "EV suppressor (soft but well-struck)",
                    "Vulnerable (allows both)")))]
cat(sprintf("\n  Split at the medians: hard-hit %.1f%%, conversion %.1f%%\n\n", hmed, cmed))
print(P[, .(pitchers=.N, hardhit=round(mean(hardhit),1), conv=round(mean(conv),1),
  barrel=round(mean(barrel),1), ev=round(mean(ev),1), la=round(mean(la),1),
  gb=round(mean(gb),1), fb=round(mean(fb),1), pop=round(mean(pop),1),
  batspd=round(mean(batspd),1), atkang=round(mean(atkang),1), tdev=round(mean(tdev),1),
  xwobacon=round(mean(xw),3), rv=round(mean(rv),4)), by=type][order(barrel)], row.names=FALSE)

show <- function(d, n=12) d[, .(player_name, bip, hardhit=round(hardhit,1),
  conv=round(conv,1), barrel=round(barrel,1), ev=round(ev,1), la=round(la,1),
  gb=round(gb,1), pop=round(pop,1), tdev=round(tdev,1), xwobacon=round(xw,3),
  rv=round(rv,4))][1:min(n,nrow(d))]
for (ty in c("Weak-contact (suppresses both)","Angle manager (allows EV, kills LA)",
             "EV suppressor (soft but well-struck)","Vulnerable (allows both)")) {
  cat(sprintf("\n  -- %s: lowest barrel%% examples --\n", ty))
  print(show(P[type==ty][order(barrel)]), row.names=FALSE)
}

# =============================================================================
# 5. Are these types repeatable?
# =============================================================================
cat("\n############ 5. Repeatability -- do contact-manager types persist? ############\n")
bip[, half := sample(rep_len(1:2, .N)), by=pitcher]
sb <- function(r) 2*r/(1+r)
MS <- c("hardhit","conv","barrel","ev","la","gb","fb","pop","batspd","atkang","tdev","xw")
LM <- c(hardhit="Hard-hit % allowed", conv="Barrel conversion (barrel|hard)",
        barrel="Barrel % allowed", ev="Mean exit velocity", la="Mean launch angle",
        gb="Ground-ball %", fb="Fly-ball %", pop="Popup %", batspd="Bat speed induced",
        atkang="Attack angle induced", tdev="Timing deviation induced",
        xw="xwOBAcon allowed")
pan <- function(A1, B1, lab) {
  m <- merge(A1, B1, by="pitcher", suffixes=c("_a","_b"))
  cat(sprintf("\n  -- %s (n=%d pitchers) --\n", lab, nrow(m)))
  o <- rbindlist(lapply(MS, function(v) {
    x <- m[[paste0(v,"_a")]]; y <- m[[paste0(v,"_b")]]
    ok <- is.finite(x) & is.finite(y)
    r <- cor(x[ok], y[ok])
    data.table(metric=LM[[v]], r=round(r,3),
               full_season=if (grepl("half", lab)) round(sb(r),3) else NA_real_)
  }))
  print(o[order(-r)], row.names=FALSE)
}
pan(agg(bip[game_year==2026 & half==1], "pitcher", 60),
    agg(bip[game_year==2026 & half==2], "pitcher", 60), "split-half within 2026")
pan(agg(bip[game_year==2025], "pitcher", 120),
    agg(bip[game_year==2026], "pitcher", 120), "2025 -> 2026")

cat("\n  Type persistence: of pitchers qualifying in both seasons, how many keep their type?\n")
P25 <- agg(bip[game_year==2025], "pitcher", 120); P26 <- agg(bip[game_year==2026], "pitcher", 120)
lab4 <- function(d) { h <- median(d$hardhit); c <- median(d$conv)
  fifelse(d$hardhit<=h & d$conv<=c, "Weak", fifelse(d$hardhit>h & d$conv<=c, "Angle",
  fifelse(d$hardhit<=h & d$conv>c, "EVsupp", "Vuln"))) }
P25[, ty := lab4(P25)]; P26[, ty := lab4(P26)]
mm <- merge(P25[, .(pitcher, ty25=ty)], P26[, .(pitcher, ty26=ty)], by="pitcher")
print(table(mm$ty25, mm$ty26))
cat(sprintf("\n  Same type both years: %.1f%% (chance = 25%%), n = %d\n",
  100*mean(mm$ty25==mm$ty26), nrow(mm)))
terc <- function(x) cut(x, quantile(x, c(0, 1/3, 2/3, 1)), labels=FALSE, include.lowest=TRUE)
for (v in c("hardhit","conv","la","gb")) {
  a <- terc(P25[[v]][match(mm$pitcher, P25$pitcher)])
  b <- terc(P26[[v]][match(mm$pitcher, P26$pitcher)])
  cat(sprintf("  %-8s quantile tercile held both years: %.1f%% (chance 33.3%%)\n",
              v, 100*mean(a == b)))
}

# =============================================================================
# 6. Payoff: the mechanisms repeat better than the rates -- do they forecast better?
# =============================================================================
cat("\n############ 6. Out-of-sample forecast of 2026 barrel% from 2025 ############\n")
cat("Everything on the right-hand side is measured in 2025 only. Same 190 pitchers.\n\n")
F <- merge(P25[, .(pitcher, hardhit, conv, barrel, la, gb, fb, pop, ev, atkang, tdev, batspd)],
           P26[, .(pitcher, barrel26 = barrel, xw26 = xw)], by="pitcher")
one <- rbindlist(lapply(c("barrel","conv","hardhit","la","gb","fb","ev","atkang","tdev","batspd","pop"),
  function(v) data.table(predictor_2025 = v, r_with_2026_barrel = round(cor(F[[v]], F$barrel26), 3))))
print(one[order(-abs(r_with_2026_barrel))], row.names=FALSE)
r2 <- function(f) summary(lm(as.formula(f), data=F))$r.squared
cat(sprintf("\n  R2 predicting 2026 barrel%%:\n"))
cat(sprintf("     2025 barrel%% alone                    = %.3f\n", r2("barrel26 ~ barrel")))
cat(sprintf("     2025 hard-hit%% alone                  = %.3f\n", r2("barrel26 ~ hardhit")))
cat(sprintf("     2025 launch angle + GB%% (mechanism)   = %.3f\n", r2("barrel26 ~ la + gb")))
cat(sprintf("     2025 hard-hit%% + launch angle + GB%%   = %.3f\n", r2("barrel26 ~ hardhit + la + gb")))
cat(sprintf("     everything above + 2025 barrel%%       = %.3f\n", r2("barrel26 ~ hardhit + la + gb + barrel")))

fwrite(P[order(barrel)], file.path("data","statcast_2026","contact_manager_types_2026.csv"))
cat("\nWrote data/statcast_2026/contact_manager_types_2026.csv\n")

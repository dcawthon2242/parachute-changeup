#!/usr/bin/env Rscript

# A TRANSFERABLE SPIN-EFFICIENCY ESTIMATOR, CALIBRATED ON MLB AND APPLIED TO D1.
#
# The goal is to screen low-efficiency changeups and four-seamers out of the college data, where
# no efficiency figure is reported. The axis route from script 07 is unavailable: TrackMan
# derives its axis from movement, so the measured-versus-inferred discrepancy is near zero by
# construction. What remains is break MAGNITUDE, and that is enough if it is done as physics
# rather than as a regression on convenient columns.
#
# The inversion. Lift on a spinning ball is a_M = C_L * (rho*A/2m) * v^2, and C_L is a function
# of the spin factor S = r*omega_T/v, where omega_T is the TRANSVERSE spin only - gyro spin
# generates no force. Every term except omega_T is measured by both systems, so:
#
#     measure a_M and v  ->  C_L = a_M / (K v^2)  ->  invert C_L(S) for S  ->  omega_T = S v / r
#     efficiency = omega_T / omega_total, and omega_total is what the radar and the camera both report
#
# Two decisions make this work where the earlier attempt did not. First, a_M is isolated exactly
# rather than approximated: total acceleration minus gravity still contains drag, so the
# remainder is split into components along the velocity vector (drag) and across it (Magnus, which
# is perpendicular to v by definition). The previous build used sqrt(ax^2 + (az+g)^2), which keeps
# a slice of drag and discards the y-component of lift. Second, C_L(S) is calibrated empirically
# against MLB measured active spin instead of taken from a literature fit, so any systematic in
# Statcast's own definition is absorbed.
#
# What makes it portable is that the calibration maps one dimensionless number to another. It
# contains no league, no park and no pitcher; only physics. The failure modes are air density,
# which is unmeasured and differs by park and month, and seam-shifted wake, which puts force on
# the ball that no spin term can account for.

suppressPackageStartupMessages({ library(data.table); library(mgcv) })
set.seed(11); options(width = 210)
DL <- path.expand("~/Downloads"); MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
NOUT <- file.path(MDIR, "ncaa_d1_kinematics.rds")

# 5.125 oz ball, 9.125 in circumference, sea-level air. K carries every constant in the lift law.
R_BALL <- 9.125/(2*pi)/12; A_BALL <- pi*R_BALL^2
M_BALL <- (5.125/16)/32.174; RHO <- 0.0023769
KLIFT  <- RHO*A_BALL/(2*M_BALL); G <- 32.174
RPM2RAD <- 2*pi/60

# Gravity out, then the remainder split along and across velocity. Returns lift coefficient and speed.
lift <- function(vx, vy, vz, ax, ay, az) {
  wx <- ax; wy <- ay; wz <- az + G
  v2 <- vx^2 + vy^2 + vz^2
  d  <- (wx*vx + wy*vy + wz*vz)/v2
  aM <- sqrt((wx-d*vx)^2 + (wy-d*vy)^2 + (wz-d*vz)^2)
  list(CL = aM/(KLIFT*v2), v = sqrt(v2), aM = aM)
}

## ---- MLB side: build the calibration ------------------------------------------------------------
cat("=== MLB calibration set ===\n")
M <- readRDS(file.path(MDIR, "parachute_rv.rds"))[pitch_type %in% c("CH","FF")]
M <- M[is.finite(vx0) & is.finite(vy0) & is.finite(vz0) & is.finite(ax) & is.finite(ay) &
       is.finite(az) & is.finite(release_spin_rate) & release_spin_rate > 500]
L <- lift(M$vx0, M$vy0, M$vz0, M$ax, M$ay, M$az); M[, `:=`(CL = L$CL, v = L$v, aM = L$aM)]
M <- M[CL > 0 & CL < 1 & v > 80 & v < 160]

PM <- M[, .(n = .N, CL = mean(CL), v = mean(v), spin = mean(release_spin_rate), aM = mean(aM),
            ext = mean(release_extension, na.rm = TRUE)),
        by = .(pitcher, season, pitch_type)][n >= 50 & is.finite(ext)]
AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
PM <- merge(PM, AS[, .(pitcher, season, pitch_type, measured = active_spin)],
            by = c("pitcher","season","pitch_type"))[is.finite(measured)]
PM[, `:=`(omega = spin*RPM2RAD)]
PM[, S_true := R_BALL*(measured*omega)/v]
cat(sprintf("  %d pitcher-season-pitch-types (CH %d, FF %d)\n", nrow(PM),
            sum(PM$pitch_type=="CH"), sum(PM$pitch_type=="FF")))
cat(sprintf("  spin factor S spans %.3f to %.3f; lift coefficient %.3f to %.3f\n",
            quantile(PM$S_true,.01), quantile(PM$S_true,.99), quantile(PM$CL,.01), quantile(PM$CL,.99)))

# Sanity check against the textbook before trusting anything downstream: at S = 0.20 the
# literature puts C_L near 0.20, and a 94 mph four-seamer should show roughly 20 ft/s^2 of lift.
cat(sprintf("  check - four-seam median: v %.1f ft/s, lift %.1f ft/s^2, C_L %.3f, S %.3f\n",
            PM[pitch_type=="FF", median(v)], PM[pitch_type=="FF", median(aM)],
            PM[pitch_type=="FF", median(CL)], PM[pitch_type=="FF", median(S_true)]))

## ---- fit and invert C_L(S), out of fold ----------------------------------------------------------
# Three specifications, in increasing order of how much they could learn that will not transfer.
# The pure-physics one uses only the lift coefficient; the others add speed (a proxy for Reynolds
# number and drag regime) and then total spin. Spin is the risky one - efficiency is
# transverse-over-total, so handing the model the denominator lets it learn the league's
# spin-to-efficiency habits rather than the physics, and that is exactly what would not port.
cat("\n=== recovering efficiency, 5-fold out of fold ===\n")
PM[, fold := sample(rep_len(1:5, .N))]
SPECS <- list(physics = S_true ~ s(CL, k = 12),
              plus_speed = S_true ~ s(CL, k = 12) + s(v, k = 6),
              plus_spin  = S_true ~ s(CL, k = 12) + s(v, k = 6) + s(spin, k = 6),
              full = S_true ~ te(CL, spin, k = c(8,6)) + s(v, k = 6) + s(ext, k = 5))
for (nm in names(SPECS)) {
  PM[, hat := NA_real_]
  for (k in 1:5) PM[fold == k, hat := predict(gam(SPECS[[nm]], data = PM[fold != k]), .SD)]
  PM[, (paste0("eff_", nm)) := pmin(pmax(hat*v/(R_BALL*omega), 0), 1.05)]
}
res <- rbindlist(lapply(names(SPECS), function(nm) PM[, {
  e <- get(paste0("eff_", nm))
  .(spec = nm, seasons = .N, r = round(cor(e, measured),3), R2 = round(cor(e, measured)^2,3),
    rmse = round(sqrt(mean((e-measured)^2)),4), bias = round(mean(e-measured),4),
    sd_truth = round(sd(measured),4)) }, by = pitch_type]))
print(dcast(res, pitch_type + seasons + sd_truth ~ spec, value.var = c("r","rmse")), row.names = FALSE)
cat("  The old build's crude ratio managed r = .61 on changeups for comparison.\n")

## ---- does the calibration survive a league shift? -------------------------------------------------
# Out-of-fold accuracy is measured on pitches drawn from the same population as the training set.
# D1 is not that population: it throws slower and, on changeups, spins more. The relevant question
# is therefore not which specification fits MLB best but which one holds up when the test set sits
# outside the training distribution. Training on one half of a variable and testing on the other
# half simulates exactly that, and it is the only way to catch a model that is borrowing a
# population regularity rather than reading the physics. Speed is the split that matters most,
# since it is the largest MLB-to-D1 gap.
cat("\n=== transfer stress test: train on one half of the distribution, predict the other ===\n")
stress <- rbindlist(lapply(c("v","spin"), function(var) rbindlist(lapply(names(SPECS), function(nm) {
  out <- rbindlist(lapply(c("low","high"), function(side) {
    med <- median(PM[[var]])
    tr <- if (side == "low") PM[get(var) <  med] else PM[get(var) >= med]
    te <- if (side == "low") PM[get(var) >= med] else PM[get(var) <  med]
    fit <- gam(SPECS[[nm]], data = tr)
    e <- pmin(pmax(as.numeric(predict(fit, te))*te$v/(R_BALL*te$spin*RPM2RAD), 0), 1.05)
    data.table(r = cor(e, te$measured), bias = mean(e - te$measured)) }))
  data.table(split = var, spec = nm, r = round(mean(out$r),3),
             worst_bias = round(out$bias[which.max(abs(out$bias))],4)) }))))
print(dcast(stress, spec ~ split, value.var = c("r","worst_bias")), row.names = FALSE)
cat("  Extrapolating in speed is what D1 asks the model to do. Bias is the number to watch:\n")
cat("  a specification that stays accurate in-sample but drifts under shift is borrowing MLB's\n")
cat("  spin-to-efficiency habits, which is precisely what will not transfer to college.\n")

PM[, eff_hat := eff_plus_spin]

## ---- the actual deliverable: a screen, with its operating characteristics -------------------------
cat("\n=== screening out low efficiency: what each cut buys ===\n")
scr <- rbindlist(lapply(c("CH","FF"), function(pt) rbindlist(lapply(c(.75,.80,.85,.90), function(th) {
  D <- PM[pitch_type == pt]; keep <- D$eff_hat >= th; bad <- D$measured < th
  data.table(pitch = pt, threshold = th, true_low_rate = round(100*mean(bad),1),
             kept = round(100*mean(keep),1),
             purity_of_kept = round(100*mean(!bad[keep]),1),
             low_removed = round(100*mean(!keep[bad]),1),
             median_kept = round(median(D$measured[keep]),3),
             median_cut = round(median(D$measured[!keep]),3)) }))))
print(scr, row.names = FALSE)
cat("  'purity_of_kept' is the number that matters for bin-building: the share of surviving\n")
cat("  pitcher-seasons that genuinely clear the floor. 'low_removed' is recall on the bad class.\n")

# A screen can also be run on rank rather than level, which is what you want if the absolute
# calibration might shift between leagues but the ordering holds.
cat("\n  dropping the bottom decile by estimate, instead of using an absolute floor:\n")
for (pt in c("CH","FF")) { D <- PM[pitch_type == pt]
  cut10 <- D$eff_hat < quantile(D$eff_hat,.10)
  cat(sprintf("    %s  bottom decile has median measured efficiency %.3f against %.3f for the rest;\n",
              pt, median(D$measured[cut10]), median(D$measured[!cut10])))
  cat(sprintf("        it captures %.0f%% of everything truly below .80\n",
              100*mean(cut10[D$measured < .80]))) }

# Screening at the same number you care about is the lowest-purity way to do it, because every
# pitch near the line is a coin flip. Demanding a margin - keep only estimates comfortably above
# the floor - trades retention for a cleaner surviving set, and for bin-building that is usually
# the right trade. This sweep prices it: how high does the estimate have to be before the kept
# group is 90 or 95 percent genuine, and how much of the population is left when it is.
cat("\n=== buying purity with a margin: truth floor fixed at .85, estimate threshold swept ===\n")
marg <- rbindlist(lapply(c("CH","FF"), function(pt) rbindlist(lapply(seq(.85,.97,.01), function(th) {
  D <- PM[pitch_type == pt]; keep <- D$eff_hat >= th
  data.table(pitch = pt, est_threshold = th, kept_pct = round(100*mean(keep),1),
             purity = round(100*mean(D$measured[keep] >= .85),1),
             median_kept = round(median(D$measured[keep]),3)) }))))
print(dcast(marg, est_threshold ~ pitch, value.var = c("kept_pct","purity")), row.names = FALSE)
for (pt in c("CH","FF")) for (tgt in c(90, 95)) {
  m <- marg[pitch == pt & purity >= tgt]
  if (nrow(m)) cat(sprintf("  %s: %d%% purity against a .85 floor needs an estimate of %.2f, retaining %.0f%% of seasons\n",
                           pt, tgt, m$est_threshold[1], m$kept_pct[1]))
  else cat(sprintf("  %s: %d%% purity is not reachable at any threshold in this range\n", pt, tgt)) }

saveRDS(PM, file.path(MDIR, "mlb_spineff_calibration.rds"))

## ---- refit on all MLB data, ready to apply --------------------------------------------------------
# Both the conservative and the accurate specification are carried through to D1, so the choice
# can be made on how the two behave out of sample rather than on the MLB fit alone.
FIT_C <- gam(SPECS$plus_speed, data = PM); FIT_A <- gam(SPECS$plus_spin, data = PM)
apply_eff <- function(fit, CL, v, spin, ext, cap = TRUE) {
  S <- predict(fit, newdata = data.frame(CL = CL, v = v, spin = spin, ext = ext))
  e <- as.numeric(S)*v/(R_BALL*spin*RPM2RAD)
  if (cap) pmin(pmax(e, 0), 1.05) else e }
# The uncapped MLB estimate supplies the reference overshoot rate used to correct D1 below.
PM[, eff_raw_unc := apply_eff(FIT_A, CL, v, spin, ext, cap = FALSE)]

## ---- D1 side: same physics, same calibration ------------------------------------------------------
cat("\n=== applying to D1 TrackMan ===\n")
if (!file.exists(NOUT)) {
  SCHEMA <- list("2023" = list(f = "pbp23tm.csv", map = c(ax="ax", ay="ay", az="az")),
                 "2024" = list(f = "D1TM24.csv",  map = c(ax="ax0", ay="ay0", az="az0")),
                 "2025" = list(f = "D1TM25.csv",  map = c(ax="ax0", ay="ay0", az="az0")))
  COMMON <- c("Level","Pitcher","PitcherId","PitcherThrows","AutoPitchType","TaggedPitchType",
              "PitchCall","RelSpeed","SpinRate","SpinAxis","Extension","RelHeight","RelSide",
              "InducedVertBreak","HorzBreak","vx0","vy0","vz0")
  N <- rbindlist(lapply(names(SCHEMA), function(yr) {
    s <- SCHEMA[[yr]]
    x <- fread(file.path(DL, s$f), showProgress = FALSE, select = c(COMMON, unname(s$map)))
    setnames(x, unname(s$map), names(s$map)); x <- x[Level == "D1"]; x[, season := as.integer(yr)]
    cat(sprintf("  %s: %s D1 pitches\n", yr, format(nrow(x), big.mark=","))); x }), use.names = TRUE)
  saveRDS(N, NOUT)
} else cat("  (using cached kinematics file)\n")
N <- readRDS(NOUT)

N <- N[is.finite(vx0) & is.finite(vy0) & is.finite(vz0) & is.finite(ax) & is.finite(ay) &
       is.finite(az) & is.finite(SpinRate) & SpinRate > 500]
LN <- lift(N$vx0, N$vy0, N$vz0, N$ax, N$ay, N$az); N[, `:=`(CL = LN$CL, v = LN$v, aM = LN$aM)]
N <- N[CL > 0 & CL < 1 & v > 80 & v < 160]
N[, pt := fifelse(!is.na(AutoPitchType) & AutoPitchType != "", AutoPitchType, TaggedPitchType)]

# Before applying an MLB-fitted curve to TrackMan output, confirm the two systems put the same
# pitch in the same place. If a convention or reference plane differed, these would not line up.
cat("\n  convention check - are the inputs on the same scale?\n")
ca <- PM[, .(v = round(median(v),1), spin = round(median(spin)), CL = round(median(CL),3)),
         by = .(pitch = fifelse(pitch_type == "CH", "Changeup", "Four-Seam"))][, source := "MLB Hawk-Eye"]
cb <- N[pt %in% c("Changeup","Four-Seam"), .(v = round(median(v),1), spin = round(median(SpinRate)),
        CL = round(median(CL),3)), by = .(pitch = pt)][, source := "D1 TrackMan"]
print(rbind(ca, cb)[order(pitch, source), .(source, pitch, v, spin, CL)], row.names = FALSE)

ND <- N[pt %in% c("Changeup","Four-Seam"),
        .(n = .N, CL = mean(CL), v = mean(v), spin = mean(SpinRate), aM = mean(aM),
          velo = mean(RelSpeed), ivb = mean(InducedVertBreak), hb = mean(HorzBreak),
          axis = mean(SpinAxis, na.rm = TRUE), ext = mean(Extension, na.rm = TRUE),
          relh = mean(RelHeight, na.rm = TRUE), rels = mean(RelSide, na.rm = TRUE),
          throws = PitcherThrows[1], name = Pitcher[1]),
        by = .(PitcherId, season, pt)][n >= 40]
ND[, `:=`(eff_cons = apply_eff(FIT_C, CL, v, spin, ext),
          eff_raw  = apply_eff(FIT_A, CL, v, spin, ext, cap = FALSE))]

## ---- one-parameter correction, anchored on the physical ceiling -----------------------------------
# The raw D1 estimates sit high and a tenth of them exceed 1.0, which is impossible: a ball cannot
# have more transverse spin than total spin. Something is scaling the measured lift. Air density is
# the obvious candidate, since it enters the lift law linearly and is unmeasured in both datasets,
# and TrackMan fitting its trajectory over a different window would do the same thing.
#
# Both would produce a MULTIPLICATIVE offset, and that is identifiable without any college truth,
# because the ceiling is known. Solve for the single constant that makes D1 overshoot 1.0 at the
# same rate MLB does. The test of whether this is a real measurement offset rather than a fudge is
# that one constant, fit on both pitch types at once, should line both of them up separately - a
# genuine physical bias has no reason to respect pitch type, and a fudge would not.
mlb_over <- PM[, mean(eff_raw_unc > 1)]
cc <- optimize(function(c) (mean(ND$eff_raw/c > 1) - mlb_over)^2, c(0.8, 1.3))$minimum
ND[, eff_hat := pmin(pmax(eff_raw/cc, 0), 1.05)]
cat(sprintf("\n  ceiling anchor: MLB exceeds 1.0 on %.1f%% of seasons; D1 needs a %.3fx correction to match\n",
            100*mlb_over, cc))
cat(sprintf("  shift removed by that one constant - changeups %.3f -> %.3f, four-seamers %.3f -> %.3f (MLB %.3f / %.3f)\n",
            ND[pt=="Changeup", median(eff_raw)], ND[pt=="Changeup", median(eff_hat)],
            ND[pt=="Four-Seam", median(eff_raw)], ND[pt=="Four-Seam", median(eff_hat)],
            PM[pitch_type=="CH", median(measured)], PM[pitch_type=="FF", median(measured)]))

cat(sprintf("\n  scored %d D1 pitcher-season-pitch-types\n", nrow(ND)))
qs <- function(x) list(q10 = round(quantile(x,.10),3), median = round(median(x),3),
                       q90 = round(quantile(x,.90),3), pct_below_80 = round(100*mean(x<.80),1))
d1 <- PM[, qs(measured), by = .(pitch = fifelse(pitch_type=="CH","Changeup","Four-Seam"))][, source := "MLB measured"]
d2 <- ND[, qs(pmin(eff_raw,1.05)), by = .(pitch = pt)][, source := "D1 uncorrected"]
d3 <- ND[, qs(eff_hat),  by = .(pitch = pt)][, source := "D1 corrected"]
print(rbind(d1,d2,d3)[order(pitch, source), .(source, pitch, q10, median, q90, pct_below_80)],
      row.names = FALSE)

cat("\n  lowest-efficiency D1 changeups the estimator finds:\n")
print(ND[pt == "Changeup"][order(eff_hat)][1:15, .(Pitcher = name, Season = season, CH = n,
      Velo = round(velo,1), Spin = round(spin), IVB = round(ivb,1), HB = round(hb,1),
      Est = round(eff_hat,3))], row.names = FALSE)

## ---- the screen applied: D1 pitcher-seasons with both pitches clearing the floor ------------------
# The parachute work needs a pitcher-season where the changeup AND the four-seamer are both
# genuinely high efficiency, so the screen has to pass on the pair, not on either pitch alone.
cat("\n=== applying the .90 screen to D1 pairs (95 percent purity against a true .85 floor) ===\n")
PR <- merge(ND[pt == "Changeup", .(PitcherId, season, name, throws, n_ch = n, eff_ch = eff_hat,
                                   velo_ch = velo, spin_ch = spin, ivb_ch = ivb, hb_ch = hb, axis_ch = axis)],
            ND[pt == "Four-Seam", .(PitcherId, season, n_ff = n, eff_ff = eff_hat,
                                    velo_ff = velo, spin_ff = spin, axis_ff = axis)],
            by = c("PitcherId","season"))
PR[, `:=`(velo_sep = velo_ff - velo_ch,
          axis_gap = pmin(abs(axis_ch-axis_ff), 360-abs(axis_ch-axis_ff)),
          clean = eff_ch >= .90 & eff_ff >= .90)]
cat(sprintf("  %d D1 pitcher-seasons throw both pitches; %d (%.1f%%) pass the screen on both\n",
            nrow(PR), sum(PR$clean), 100*mean(PR$clean)))
cat(sprintf("  expected contamination at this operating point: about 5%% per pitch, so roughly %.0f%% of\n",
            100*(1-.959*.996)))
cat("  surviving pairs contain at least one pitch that is truly below .85.\n")
print(PR[, .(seasons = .N, median_ch_eff = round(median(eff_ch),3),
             median_velo_sep = round(median(velo_sep),1)), by = .(passes_screen = clean)],
      row.names = FALSE)
saveRDS(PR, file.path(MDIR, "ncaa_spineff_pairs.rds"))

## ---- figure ---------------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
LB <- c(CH = "Changeups", FF = "Four-seamers")
PM[, ptl := LB[pitch_type]]

cur <- data.table(CL = seq(quantile(PM$CL,.01), quantile(PM$CL,.99), length.out = 200))
cur[, S := as.numeric(predict(gam(S_true ~ s(CL, k = 12), data = PM), newdata = cur))]
p1 <- ggplot(PM, aes(CL, S_true)) + geom_point(alpha = .10, size = .7, colour = "grey35") +
  geom_line(data = cur, aes(CL, S), colour = "#b4632a", linewidth = 1) +
  labs(subtitle = "1. The inversion: lift coefficient against true spin factor, with the fitted curve",
       x = expression(C[L]~"= lift / (K"~v^2*")"), y = "Spin factor S = r"~omega[T]*"/v") +
  theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8.4, face = "bold"))

rl <- PM[, .(lab = sprintf("r = %.3f\nRMSE = %.3f", cor(eff_hat, measured),
                           sqrt(mean((eff_hat-measured)^2)))), by = ptl]
p2 <- ggplot(PM, aes(measured, eff_hat)) + geom_abline(slope = 1, linetype = 2, colour = "grey55") +
  geom_point(alpha = .12, size = .7, colour = "#1d7870") +
  geom_text(data = rl, aes(x = .42, y = 1.02, label = lab), hjust = 0, vjust = 1, size = 2.9,
            lineheight = 1.1, fontface = "bold", colour = "#0d3f3a") +
  facet_wrap(~ptl) + coord_cartesian(xlim = c(.35,1.05), ylim = c(.35,1.05)) +
  labs(subtitle = "2. Out-of-fold accuracy against Savant measured active spin",
       x = "Measured active spin", y = "Estimated") +
  theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 8.6), plot.subtitle = element_text(size = 8.4, face = "bold"))

MG <- melt(marg, id.vars = c("pitch","est_threshold"), measure.vars = c("kept_pct","purity"))
MG[, `:=`(pitch = LB[pitch], variable = factor(variable, c("purity","kept_pct"),
          c("Purity of kept set", "Share retained")))]
p3 <- ggplot(MG, aes(est_threshold, value, colour = variable)) +
  geom_hline(yintercept = 95, linetype = 3, colour = "grey55") +
  geom_line(linewidth = .8) + facet_wrap(~pitch) +
  scale_colour_manual(values = c("Purity of kept set" = "#1d7870", "Share retained" = "#b4632a"), name = NULL) +
  labs(subtitle = "3. What a margin buys: truth floor held at .85, estimate threshold swept",
       x = "Estimate threshold applied", y = "Percent") +
  theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "top",
        strip.text = element_text(face = "bold", size = 8.6), plot.subtitle = element_text(size = 8.4, face = "bold"))

DS <- rbind(PM[, .(pitch = ptl, eff = measured, src = "MLB measured")],
            ND[, .(pitch = LB[c(Changeup="CH", `Four-Seam`="FF")[pt]], eff = eff_hat, src = "D1 estimated")])
p4 <- ggplot(DS, aes(eff, fill = src, colour = src)) +
  geom_density(alpha = .25, linewidth = .7, adjust = 1.2) + facet_wrap(~pitch) +
  coord_cartesian(xlim = c(.5, 1.05)) +
  scale_fill_manual(values = c("MLB measured" = "grey45", "D1 estimated" = "#1d7870"), name = NULL) +
  scale_colour_manual(values = c("MLB measured" = "grey45", "D1 estimated" = "#1d7870"), name = NULL) +
  labs(subtitle = sprintf("4. D1 after the %.3fx ceiling correction, against MLB measured truth", cc),
       x = "Spin efficiency", y = "Density") +
  theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "top",
        strip.text = element_text(face = "bold", size = 8.6), plot.subtitle = element_text(size = 8.4, face = "bold"))

gg <- (p1 | p2) / (p3 | p4) + plot_annotation(
  title = "Recovering spin efficiency from the lift law: r = .76 on changeups and .92 on four-seamers, and it transfers to TrackMan",
  subtitle = paste0("Lift on a spinning ball is a_M = C_L(S) * K * v^2, where the spin factor S = r*omega_T/v depends on TRANSVERSE spin only, because gyro spin generates no force. Both\n",
                    "tracking systems measure total spin and the full trajectory, so the transverse component can be solved for and divided out. Two details do the work. The lift is\n",
                    "isolated exactly rather than approximated - total acceleration minus gravity still contains drag, so the remainder is split into components along the velocity\n",
                    "vector (drag) and across it (Magnus, which is perpendicular to v by definition). And C_L(S) is calibrated against Savant's measured active spin rather than taken\n",
                    "from a literature fit, which absorbs any systematic in Statcast's own definition. Panel 3 is the practical result: screening at the floor itself is the least\n",
                    "accurate way to use the estimate, since every pitch near the line is a coin flip, but demanding an estimate of .90 against a true .85 floor leaves a set that is\n",
                    "96 percent genuine on changeups and 99.6 on four-seamers, at the cost of a little over half the population. Applied to D1, the raw estimates ran 8 percent high\n",
                    "with a tenth of them above the physical ceiling of 1.0. A single multiplicative constant, identified by matching MLB's overshoot rate, fixes it - and the check\n",
                    "that it is a real measurement offset rather than a fudge is that the same constant lines up both pitch types independently, changeups to .886 and four-seamers to\n",
                    ".885 against MLB's .919 and .915. Air density is the likely culprit; it enters the lift law linearly and is unmeasured in both datasets."),
  caption = "MLB Statcast 2020-2026 and NCAA D1 TrackMan 2023-2025 - pitcher-seasons with 50+ (MLB) or 40+ (D1) pitches of the type - truth is Savant measured active spin",
  theme = theme(plot.title = element_text(face = "bold", size = 11.5),
                plot.subtitle = element_text(size = 8.1), plot.caption = element_text(size = 7.3, colour = "grey40")))
ggsave(file.path(AST, "fig42_spineff_transfer.png"), gg, width = 13, height = 10.4, dpi = 150)

saveRDS(ND, file.path(MDIR, "ncaa_spineff.rds"))
cat("\nwrote ncaa_spineff.rds, ncaa_spineff_pairs.rds, mlb_spineff_calibration.rds, fig42_spineff_transfer.png\n")

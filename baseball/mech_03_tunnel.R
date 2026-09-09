#!/usr/bin/env Rscript

# M3: DOES THE ADVANTAGE CONCENTRATE ON WELL-TUNNELED CHANGEUPS?
#
# The third way to separate mimicry from pitcher quality. path_to_location_ratio is the integrated
# separation between this pitch's flight path and the previous pitch's, over the window before the
# hitter has to commit, divided by how far apart the two finish at the plate. Low means the two
# pitches looked alike while they mattered and diverged afterwards.
#
# If the bin's changeup works by being mistaken for the four-seamer, the advantage should be
# concentrated on the well-tunneled ones. Pitcher-season fixed effects absorb the bin main effect,
# so the identifying variation is within a single pitcher's own season: his best-tunneled changeups
# against his worst.
#
# The D1 trajectory columns live in a separate file from the pitch table. They are row-aligned by
# construction, both being reads of the same three CSVs in the same order, but that is asserted and
# checked rather than assumed - a silent off-by-one would produce a plausible-looking ratio computed
# from another pitcher's pitch.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 200); set.seed(17)
MDIR <- "data/statcast_model"
L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L); L[, id := as.character(id)]

WX <- 1.2; WY <- 0.4; WZ <- 1.4; REACT <- 0.150; NSTEP <- 40
pos <- function(r, v, a, t) r + v*t + 0.5*a*t^2

# The path integral, weighted to match how a hitter actually resolves a pitch: vertical separation
# reads most strongly, depth least. Weights and the 150 ms commit window are carried over unchanged
# from the tunneling project so this is the same quantity used everywhere else in the analysis.
tunnel_int <- function(d, sfx = "p_") {
  g <- function(n) d[[n]]; p <- function(n) d[[paste0(sfx, n)]]
  Tmax <- pmax(g("t_react"), p("t_react")); acc <- numeric(nrow(d))
  for (k in 1:NSTEP) {
    tk <- (k - 0.5)/NSTEP * Tmax
    dx <- pos(g("rx"), g("vx0"), g("ax"), tk) - pos(p("rx"), p("vx0"), p("ax"), tk)
    dy <- pos(g("ry"), g("vy0"), g("ay"), tk) - pos(p("ry"), p("vy0"), p("ay"), tk)
    dz <- pos(g("rz"), g("vz0"), g("az"), tk) - pos(p("rz"), p("vz0"), p("az"), tk)
    acc <- acc + sqrt(WX*dx^2 + WY*dy^2 + WZ*dz^2) * (Tmax/NSTEP)
  }
  acc
}

## ---- MLB tunnel ratios ---------------------------------------------------------------------------
MT <- file.path(MDIR, "mech_tunnel_mlb.rds")
if (!file.exists(MT)) {
  X <- readRDS(file.path(MDIR, "parachute_rv.rds")); setDT(X)
  setnames(X, c("release_pos_x","release_pos_y","release_pos_z"), c("rx","ry","rz"))
  X[, t_plate := (-vy0 - sqrt(vy0^2 - 2*ay*(ry - 17/12)))/ay]
  X[, t_react := pmax(t_plate - REACT, 0.05)]
  setorder(X, game_pk, at_bat_number, pitch_number)
  lag <- c("rx","ry","rz","vx0","vy0","vz0","ax","ay","az","t_react","pitch_number",
           "plate_x","plate_z","pitch_type")
  for (c in lag) X[, (paste0("p_", c)) := shift(get(c)), by = .(game_pk, at_bat_number)]
  X <- X[pitch_type == "CH" & !is.na(p_pitch_number) & pitch_number - p_pitch_number == 1L]
  X[, tun := tunnel_int(X)]
  X[, sep := pmax(sqrt((plate_x - p_plate_x)^2 + (plate_z - p_plate_z)^2), 0.1)]
  X[, pr := tun/sep]
  saveRDS(X[, .(game_pk, at_bat_number, pitch_number, prev_pt = p_pitch_type, pr)], MT)
}
MP <- readRDS(MT); setDT(MP)

## ---- D1 tunnel ratios ------------------------------------------------------------------------------
DT <- file.path(MDIR, "mech_tunnel_d1.rds")
if (!file.exists(DT)) {
  D <- readRDS(file.path(MDIR, "ncaa_d1_seq_pitches.rds")); setDT(D)
  K <- readRDS(file.path(MDIR, "ncaa_d1_kinematics.rds")); setDT(K)
  stopifnot(nrow(D) == nrow(K))
  align <- mean(abs(D$RelSpeed - K$RelSpeed) < 1e-6, na.rm = TRUE)
  ralign <- mean(abs(D$RelHeight - K$RelHeight) < 1e-6, na.rm = TRUE)
  cat(sprintf("D1 kinematics row alignment: RelSpeed %.5f, RelHeight %.5f of rows identical\n",
              align, ralign))
  stopifnot(align > 0.999, ralign > 0.999)
  D[, `:=`(vx0 = K$vx0, vy0 = K$vy0, vz0 = K$vz0, ay = K$ay)]
  rm(K); gc()
  # TrackMan reports no release y, but extension off the 60.5 ft rubber gives it directly.
  D[, `:=`(rx = RelSide, rz = RelHeight, ry = 60.5 - Extension)]
  D[, t_plate := (-vy0 - sqrt(vy0^2 - 2*ay*(ry - 17/12)))/ay]
  D[, t_react := pmax(t_plate - REACT, 0.05)]
  setorder(D, game_id, Inning, PAofInning, PitchofPA, PitchNo)
  lag <- c("rx","ry","rz","vx0","vy0","vz0","ax","ay","az","t_react","PitchofPA","px","pz","pt")
  for (c in lag) D[, (paste0("p_", c)) := shift(get(c)),
                   by = .(game_id, Inning, PAofInning, BatterId)]
  D <- D[pt == "Changeup" & !is.na(p_PitchofPA) & is.finite(t_react) & is.finite(p_t_react)]
  D[, tun := tunnel_int(D)]
  D[, sep := pmax(sqrt((px - p_px)^2 + (pz - p_pz)^2), 0.1)]
  D[, pr := tun/sep]
  saveRDS(D[, .(PitchUID, prev_pt = p_pt, pr)], DT); rm(D); gc()
}
DP <- readRDS(DT); setDT(DP)

## ---- join to residuals -----------------------------------------------------------------------------
MS <- readRDS(file.path(MDIR, "mlb_whiff_resid_seq.rds")); setDT(MS)
MS <- merge(MS, MP[, .(game_pk, at_bat_number, pitch_number, pr)],
            by = c("game_pk","at_bat_number","pitch_number"))
MS[, `:=`(league = "MLB", id = as.character(pitcher))]
DS <- readRDS(file.path(MDIR, "ncaa_whiff_resid_seq.rds")); setDT(DS)
DS <- merge(DS, DP[, .(PitchUID, pr)], by = "PitchUID")
DS[, `:=`(league = "D1", id = as.character(PitcherId), prev_ff = prev_pt == "Four-Seam")]
K <- rbind(MS[, .(league, id, season, prev_ff, whiff, r_all, pr)],
           DS[, .(league, id, season, prev_ff, whiff, r_all, pr)])
K <- merge(K, L[, .(league, id, season, bin)], by = c("league","id","season"))
K <- K[!is.na(id) & is.finite(pr) & pr > 0 & pr < quantile(pr, .999, na.rm = TRUE)]

cat(sprintf("\npooled changeup swings with a tunnel ratio: %s\n", format(nrow(K), big.mark = ",")))
print(K[, .(swings = .N, arms = uniqueN(id), median_pr = round(median(pr),2),
            p10 = round(quantile(pr,.10),2), p90 = round(quantile(pr,.90),2)),
        by = .(league, bin)][order(league, -bin)], row.names = FALSE)

## ---- the interaction ---------------------------------------------------------------------------------
# path_ratio is standardised within league because the two tracking systems put it on different
# scales, and sign-flipped so that positive means better tunnelled. That way a positive interaction
# is the mechanism's prediction in both leagues without further translation.
K[, tz := -scale(log(pr))[,1], by = league]
K[, `:=`(y = 100*r_all, g = paste(league, id, season))]
K[, `:=`(yd = y - mean(y), td = tz - mean(tz)), by = g]
K[, bt := as.numeric(bin) * td]
K[, lgt := as.numeric(league == "MLB") * td]

clus <- function(fit, cl) {
  u <- residuals(fit); Xm <- model.matrix(fit); bread <- solve(crossprod(Xm))
  meat <- crossprod(rowsum(Xm * u, cl)); nc <- length(unique(cl)); k <- ncol(Xm)
  sqrt(diag(bread %*% meat %*% bread) * (nc/(nc-1)) * ((nrow(Xm)-1)/(nrow(Xm)-k))) }

fit <- lm(yd ~ 0 + td + lgt + bt, data = K); se <- clus(fit, K$id); co <- coef(fit)
cat("\n=== M3: bin x tunnel quality, pitcher-season fixed effects, clustered on pitcher ===\n")
for (r in names(co)) { t <- co[r]/se[r]
  cat(sprintf("  %-4s %+7.3f  se %5.3f  t %+5.2f  p %.4f\n", r, co[r], se[r], t, 2*pnorm(-abs(t)))) }
cat(sprintf("\n  bt is the estimate: %+.2f whiff points per standard deviation of better tunnelling,\n",
            co["bt"]))
cat(sprintf("  over and above the same slope in every non-bin arm. Detectable at 80%% power: %.1f\n",
            2.8*se["bt"]))

cat("\n=== per league ===\n")
for (lgv in c("MLB","D1")) { z <- K[league == lgv]
  f2 <- lm(yd ~ 0 + td + bt, data = z); s2 <- clus(f2, z$id)
  cat(sprintf("  %-4s n=%s (bin %s)  base slope %+6.3f  interaction %+6.3f se %.3f  p %.3f\n", lgv,
              format(nrow(z), big.mark = ","), format(sum(z$bin), big.mark = ","),
              coef(f2)["td"], coef(f2)["bt"], s2["bt"], 2*pnorm(-abs(coef(f2)["bt"]/s2["bt"])))) }

## ---- and restricted to the case the mechanism actually describes ----------------------------------
# The tunnel story is about a changeup masquerading as the fastball, so the sharpest version of the
# test uses only changeups that did follow a four-seamer.
cat("\n=== the same interaction, restricted to changeups thrown after a four-seamer ===\n")
Z <- K[prev_ff == TRUE]
Z[, `:=`(yd = y - mean(y), td2 = tz - mean(tz)), by = g]
Z[, bt2 := as.numeric(bin) * td2][, lgt2 := as.numeric(league == "MLB") * td2]
f3 <- lm(yd ~ 0 + td2 + lgt2 + bt2, data = Z); s3 <- clus(f3, Z$id)
cat(sprintf("  n=%s (bin %s, %d arms)  interaction %+.3f  se %.3f  p %.3f\n",
            format(nrow(Z), big.mark = ","), format(sum(Z$bin), big.mark = ","),
            uniqueN(Z[bin == TRUE]$id), coef(f3)["bt2"], s3["bt2"],
            2*pnorm(-abs(coef(f3)["bt2"]/s3["bt2"]))))
saveRDS(K, file.path(MDIR, "mech_tunnel.rds"))
cat("\nwrote mech_tunnel.rds\n")

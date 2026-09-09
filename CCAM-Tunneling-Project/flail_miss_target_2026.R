#!/usr/bin/env Rscript

# Adaptation of the CCAM tunneling project: instead of the swing-decision run
# value, use FLAIL rate and MISS DISTANCE as the target.
#
# The original tunneling metric integrates the (perception-weighted) 3D distance
# between consecutive pitch trajectories from release to the hitter's reaction
# point. Here we reconstruct each pitch's trajectory analytically from the 9
# Statcast kinematic params (release_pos + v0 + a0), rebuild the same metric with
# the project's optimal weights (x=1.2, y=0.4, z=1.4), then test whether better
# tunneling (smaller trajectory distance) predicts bigger misses / more flail.

suppressPackageStartupMessages({ library(data.table) })

pbp_csv   <- file.path("..", "data", "statcast_2026", "statcast_2026_all.csv")
flail_csv <- file.path("..", "data", "statcast_2026", "savant_swing_timing_2026.csv")

WX <- 1.2; WY <- 0.4; WZ <- 1.4      # perception weights from the project README
REACT <- 0.150                        # reaction point: 150 ms before plate
NSTEP <- 60                           # trajectory integration steps
FASTBALLS  <- c("FF","SI","FC")
BREAKING   <- c("SL","ST","CU","KC","SV","CS")
OFFSPEED   <- c("CH","FS","FO")

dt <- fread(pbp_csv, showProgress = FALSE, select = c(
  "game_pk","at_bat_number","pitch_number","pitcher","pitch_type","description",
  "balls","strikes","stand","p_throws","miss_distance","plate_x","plate_z",
  "release_pos_x","release_pos_y","release_pos_z","vx0","vy0","vz0","ax","ay","az","game_type"))
dt <- dt[game_type=="R" & !is.na(vx0) & pitch_type!="" & !is.na(release_pos_y)]
setorder(dt, game_pk, at_bat_number, pitch_number)

# time to plate (from release y) and reaction-point time
yf <- 17/12
dt[, t_plate := (-vy0 - sqrt(vy0^2 - 2*ay*(release_pos_y - yf)))/ay]
dt[, t_react := pmax(t_plate - REACT, 0.05)]

# previous pitch in the same at-bat (only consecutive pitch numbers)
lagcols <- c("release_pos_x","release_pos_z","release_pos_y","vx0","vy0","vz0",
             "ax","ay","az","t_react","pitch_type","pitch_number","plate_x","plate_z")
for (c in lagcols) dt[, (paste0("p_",c)) := shift(get(c)), by=.(game_pk, at_bat_number)]
dt[, has_prev := !is.na(p_pitch_number) & (pitch_number - p_pitch_number == 1)]

pos <- function(r,v,a,t) r + v*t + 0.5*a*t^2

# integrate weighted trajectory distance between this pitch and the previous one
dt[, tunnel := NA_real_]
idx <- which(dt$has_prev)
Tmax <- pmax(dt$t_react[idx], dt$p_t_react[idx])
acc <- numeric(length(idx))
for (k in 1:NSTEP) {
  tt <- (k-0.5)/NSTEP * Tmax
  dx <- pos(dt$release_pos_x[idx], dt$vx0[idx], dt$ax[idx], tt) -
        pos(dt$p_release_pos_x[idx], dt$p_vx0[idx], dt$p_ax[idx], tt)
  dy <- pos(dt$release_pos_y[idx], dt$vy0[idx], dt$ay[idx], tt) -
        pos(dt$p_release_pos_y[idx], dt$p_vy0[idx], dt$p_ay[idx], tt)
  dz <- pos(dt$release_pos_z[idx], dt$vz0[idx], dt$az[idx], tt) -
        pos(dt$p_release_pos_z[idx], dt$p_vz0[idx], dt$p_az[idx], tt)
  acc <- acc + sqrt(WX*dx^2 + WY*dy^2 + WZ*dz^2) * (Tmax/NSTEP)
}
dt$tunnel[idx] <- acc
# lower tunnel = tighter overlap = better disguise

# ---- path_to_location_ratio (reconstructed) --------------------------------
# path  = integrated trajectory distance up to reaction point (how alike along the way)
# location = Pythagorean gap between the two pitches AS THEY CROSS THE PLATE (plate_x, plate_z)
# ratio = path / location  -> LOW ratio = looked alike but ended far apart = elite tunnel/deception
dt[, location := sqrt((plate_x - p_plate_x)^2 + (plate_z - p_plate_z)^2)]
dt[, location := pmax(location, 0.1)]                    # guard against divide-by-~0
dt[, path_to_location_ratio := tunnel / location]
# BP-style cross-check: reconstructed commit-point separation / same plate location
wdist_commit <- function(i, tt) {
  dx <- pos(dt$release_pos_x[i], dt$vx0[i], dt$ax[i], tt) -
        pos(dt$p_release_pos_x[i], dt$p_vx0[i], dt$p_ax[i], tt)
  dz <- pos(dt$release_pos_z[i], dt$vz0[i], dt$az[i], tt) -
        pos(dt$p_release_pos_z[i], dt$p_vz0[i], dt$p_az[i], tt)
  sqrt(dx^2 + dz^2)
}
dt[, sep_commit := NA_real_]
dt$sep_commit[idx] <- wdist_commit(idx, dt$t_react[idx])
dt[, commit_to_plate_ratio := sep_commit / location]
# lower ratio = better tunnel

whiff_desc <- c("swinging_strike","swinging_strike_blocked","foul_tip","missed_bunt")
dt[, is_whiff := description %in% whiff_desc]
dt[, seq_class := fifelse(p_pitch_type %in% FASTBALLS & pitch_type %in% BREAKING, "breaking-after-FB",
                  fifelse(p_pitch_type %in% FASTBALLS & pitch_type %in% OFFSPEED, "offspeed-after-FB",
                  fifelse(p_pitch_type %in% FASTBALLS & pitch_type %in% FASTBALLS, "FB-after-FB", "other")))]
dt[, putaway := strikes == 2]

cat(sprintf("Pitches: %d | with tunnel metric: %d | whiffs w/ miss_distance & tunnel: %d\n",
            nrow(dt), sum(!is.na(dt$tunnel)),
            dt[!is.na(tunnel) & is_whiff & !is.na(miss_distance), .N]))

# ============================================================================
# TARGET 1: MISS DISTANCE (per whiff)
# ============================================================================
md <- dt[!is.na(path_to_location_ratio) & is.finite(path_to_location_ratio) &
         is_whiff & !is.na(miss_distance) & putaway==TRUE]     # 2-STRIKE COUNTS ONLY
cat("\n=== 2-STRIKE COUNTS ONLY ===\n")
cat("=== TARGET = MISS DISTANCE (whiffs). PREDICTOR = path_to_location_ratio ===\n")
cat("    (lower ratio = better tunnel; hypothesis: better tunnel -> bigger miss = NEGATIVE r)\n")
report <- function(d, lab, pred="path_to_location_ratio"){
  if (nrow(d) < 40) { cat(sprintf("%-30s n=%d (too few)\n", lab, nrow(d))); return(invisible()) }
  x <- d[[pred]]
  r <- cor(x, d$miss_distance, use="complete.obs")
  fit <- lm(d$miss_distance ~ x)
  p <- summary(fit)$coefficients[2,4]; b <- coef(fit)[2]
  cat(sprintf("%-30s n=%5d | r=%+.3f | slope=%+.3f | p=%.2g\n", lab, nrow(d), r, b, p))
}
report(md, "ALL whiffs")
for (sc in c("breaking-after-FB","offspeed-after-FB","FB-after-FB")) report(md[seq_class==sc], sc)
cat("\nCross-check with commit_to_plate_ratio (BP-style):\n")
report(md, "ALL whiffs", pred="commit_to_plate_ratio")
report(md[seq_class=="breaking-after-FB"], "breaking-after-FB", pred="commit_to_plate_ratio")
cat("\nControlling for pitch type (miss ~ ratio + pitch_type), breaking/offspeed after FB:\n")
sub <- md[seq_class %in% c("breaking-after-FB","offspeed-after-FB")]
fit <- lm(miss_distance ~ path_to_location_ratio + factor(pitch_type), data=sub)
print(round(summary(fit)$coefficients[1:2,], 4))

# ============================================================================
# TARGET 2: FLAIL rate (aggregate pitcher x pitch type)
# ============================================================================
fl <- fread(flail_csv, showProgress=FALSE)[, .(pitcher=id, pitch_type=api_pitch_type,
             flail=flailed_percent, st_n=n_swings)]
# mean tunnel per pitcher x pitch type, restricted to secondary-after-fastball seqs
agg <- dt[is.finite(path_to_location_ratio) & putaway==TRUE &
          seq_class %in% c("breaking-after-FB","offspeed-after-FB"),
          .(ratio=mean(path_to_location_ratio), n=.N), by=.(pitcher, pitch_type)][n>=15]
fa <- merge(agg, fl, by=c("pitcher","pitch_type"))
cat(sprintf("\n=== TARGET = FLAIL%% vs path_to_location_ratio (2-strike seqs, n=%d). Better tunnel (lower ratio) -> more flail = NEGATIVE r ===\n", nrow(fa)))
cat("    (note: flail leaderboard is all-count; ratio measured on 2-strike sequences)\n")
cat(sprintf("  overall: r = %+.3f (p=%.2g)\n",
            cor(fa$ratio, fa$flail), cor.test(fa$ratio, fa$flail)$p.value))
for (pt in c("CH","FS","SL","ST","CU")) {
  s <- fa[pitch_type==pt]
  if (nrow(s) >= 25) cat(sprintf("  %-3s: r = %+.3f (n=%d, p=%.2g)\n", pt,
      cor(s$ratio, s$flail), nrow(s), cor.test(s$ratio, s$flail)$p.value))
}
fwrite(dt[!is.na(tunnel), .(game_pk, at_bat_number, pitch_number, pitcher, pitch_type,
        p_pitch_type, seq_class, putaway, tunnel, is_whiff, miss_distance)],
       file.path("..","data","statcast_2026","tunnel_metric_2026.csv"))
cat("\nWrote ../data/statcast_2026/tunnel_metric_2026.csv\n")

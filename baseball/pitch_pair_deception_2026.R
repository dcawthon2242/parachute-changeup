#!/usr/bin/env Rscript

# Generalize the "invisiball" idea to all MLB arsenals (2026).
#
# Hypothesis: a secondary pitch that shares its FASTBALL's spin signature
# (similar spin axis + spin rate) but diverges in movement/velocity fools
# hitters -> whiffs ABOVE what a stuff model would predict. The mirror case
# (spin axis ~180 deg opposite, e.g. Glasnow's backspin FF vs topspin curve)
# is the north-south analogue.
#
# Method:
#   1. Anchor each pitcher on their PRIMARY FASTBALL (most-thrown of FF/SI/FC).
#   2. For every other pitch (secondary) with enough usage, measure its
#      relationship to that fastball: spin-axis diff, spin-rate diff, velo
#      separation, movement separation, release/tunnel similarity.
#   3. Fit a "stuff" model of Whiff% from PHYSICAL traits + separations from the
#      FB (velocity, movement, spin, release) -- but NOT the spin-axis
#      relationship. The residual (actual - predicted whiff) is whiff value the
#      stuff model can't see.
#   4. Test whether the spin-axis relationship (matched ~0 deg, or mirrored
#      ~180 deg) predicts that residual, and rank the biggest overperformers.
#
# Also reports: vertical separation from the FB (the geometric "swing over the
# top" magnitude) and bat-tracking swing geometry on whiffs (attack angle,
# swing path tilt) for the example pitches.

suppressPackageStartupMessages({ library(data.table) })

in_csv   <- file.path("data", "statcast_2026", "statcast_2026_all.csv")
out_dir  <- file.path("data", "statcast_2026")

REGULAR_ONLY <- TRUE
MIN_PRIMARY  <- 100   # min fastballs to anchor
MIN_SECOND   <- 40    # min secondary pitches to evaluate
FASTBALLS    <- c("FF", "SI", "FC")
AXIS_MATCH   <- 25    # deg; <= this = "same spin as FB" (matched)
AXIS_MIRROR  <- 155   # deg; >= this = "mirror" (opposite spin)

# ---- Load & flags ----------------------------------------------------------

dt <- data.table::fread(in_csv, showProgress = FALSE, select = c(
  "pitch_type","game_type","pitcher","player_name","p_throws",
  "release_speed","release_spin_rate","spin_axis","pfx_x","pfx_z",
  "release_pos_x","release_pos_z","release_extension",
  "description","zone","bb_type","launch_speed",
  "estimated_woba_using_speedangle","delta_run_exp",
  "bat_speed","attack_angle","swing_path_tilt",
  "intercept_ball_minus_batter_pos_y_inches"))

if (REGULAR_ONLY) dt <- dt[game_type == "R"]
dt <- dt[!is.na(pfx_x) & !is.na(pfx_z) & pitch_type != ""]

swing_desc <- c("swinging_strike","swinging_strike_blocked","foul","foul_tip",
                "hit_into_play","foul_bunt","missed_bunt","bunt_foul_tip")
whiff_desc <- c("swinging_strike","swinging_strike_blocked","foul_tip","missed_bunt")
dt[, is_swing  := description %in% swing_desc]
dt[, is_whiff  := description %in% whiff_desc]
dt[, is_called := description == "called_strike"]
dt[, out_zone  := !is.na(zone) & zone >= 11]

cmean <- function(a){ a<-a[!is.na(a)]; if(!length(a)) return(NA_real_)
  r<-a*pi/180; ((atan2(mean(sin(r)),mean(cos(r)))*180/pi)+360)%%360 }
circd <- function(a,b){ d<-abs(a-b)%%360; pmin(d,360-d) }

# ---- Per pitcher x pitch type ---------------------------------------------

agg <- dt[, .(
  n     = .N,
  velo  = mean(release_speed, na.rm=TRUE),
  spin  = mean(release_spin_rate, na.rm=TRUE),
  axis  = cmean(spin_axis),
  ivb   = mean(pfx_z)*12,
  hb    = mean(pfx_x)*12,
  relx  = mean(release_pos_x, na.rm=TRUE),
  relz  = mean(release_pos_z, na.rm=TRUE),
  ext   = mean(release_extension, na.rm=TRUE),
  swings = sum(is_swing), whiffs = sum(is_whiff), called = sum(is_called),
  oz = sum(out_zone), ozsw = sum(out_zone & is_swing),
  xwobacon = mean(estimated_woba_using_speedangle[!is.na(estimated_woba_using_speedangle)]),
  # whiff-only swing geometry
  wh_attack = mean(attack_angle[is_whiff], na.rm=TRUE),
  wh_tilt   = mean(swing_path_tilt[is_whiff], na.rm=TRUE),
  wh_batspeed = mean(bat_speed[is_whiff], na.rm=TRUE)
), by = .(pitcher, player_name, p_throws, pitch_type)]

# ---- Primary fastball per pitcher -----------------------------------------

fb <- agg[pitch_type %in% FASTBALLS & n >= MIN_PRIMARY]
setorder(fb, pitcher, -n)
primary <- fb[, .SD[1], by = pitcher]
setnames(primary, c("pitch_type","n","velo","spin","axis","ivb","hb","relx","relz","ext"),
         c("fb_type","fb_n","fb_velo","fb_spin","fb_axis","fb_ivb","fb_hb",
           "fb_relx","fb_relz","fb_ext"))
primary <- primary[, .(pitcher, fb_type, fb_n, fb_velo, fb_spin, fb_axis,
                       fb_ivb, fb_hb, fb_relx, fb_relz, fb_ext)]

# ---- Secondary pitches vs their primary fastball ---------------------------

sec <- agg[n >= MIN_SECOND]
pairs <- merge(sec, primary, by = "pitcher")
pairs <- pairs[pitch_type != fb_type]   # exclude the anchor itself

pairs[, axis_diff := circd(axis, fb_axis)]
pairs[, spin_diff := fb_spin - spin]          # + = secondary spins slower
pairs[, velo_diff := fb_velo - velo]          # + = secondary slower
pairs[, vert_sep  := fb_ivb - ivb]            # + = secondary drops below FB tunnel
pairs[, horz_sep  := hb - fb_hb]
pairs[, move_sep  := sqrt((fb_ivb-ivb)^2 + (fb_hb-hb)^2)]
pairs[, rel_diff  := sqrt((relx-fb_relx)^2 + (relz-fb_relz)^2)*12]  # inches
pairs[, ext_diff  := fb_ext - ext]
pairs[, whiff_pct := whiffs / swings * 100]
pairs[, chase_pct := ozsw / oz * 100]
pairs[, csw_pct   := (called + whiffs) / n * 100]

# ---- Stuff model: Whiff% ~ physical traits (NOT spin-axis relationship) ----

model_types <- c("FF","SI","FC","SL","ST","CU","KC","CH","FS","SV","CS")
mp <- pairs[pitch_type %in% model_types & is.finite(whiff_pct) & swings >= 15 &
            is.finite(rel_diff) & is.finite(ext_diff)]
mp[, pt := factor(pitch_type)]
stuff <- lm(whiff_pct ~ pt + velo + ivb + hb + spin + velo_diff +
              I(abs(vert_sep)) + I(abs(horz_sep)) + rel_diff + ext_diff,
            data = mp, weights = n)
mp[, whiff_resid := residuals(stuff)]
pairs <- merge(pairs, mp[, .(pitcher, pitch_type, whiff_resid)],
               by = c("pitcher","pitch_type"), all.x = TRUE)

# ---- Does the spin-axis relationship explain the residual? -----------------

valid <- pairs[is.finite(whiff_resid)]
matched <- valid[axis_diff <= AXIS_MATCH]
mirror  <- valid[axis_diff >= AXIS_MIRROR]
middle  <- valid[axis_diff > AXIS_MATCH & axis_diff < AXIS_MIRROR]

cat("=== Whiff-above-stuff residual by spin-axis relationship to fastball ===\n")
cat(sprintf("  matched spin (<=%d deg): n=%3d  mean whiff resid = %+.2f\n",
            AXIS_MATCH, nrow(matched), weighted.mean(matched$whiff_resid, matched$n)))
cat(sprintf("  middle       (%d-%d):    n=%3d  mean whiff resid = %+.2f\n",
            AXIS_MATCH, AXIS_MIRROR, nrow(middle), weighted.mean(middle$whiff_resid, middle$n)))
cat(sprintf("  mirror spin  (>=%d deg): n=%3d  mean whiff resid = %+.2f\n",
            AXIS_MIRROR, nrow(mirror), weighted.mean(mirror$whiff_resid, mirror$n)))
cat(sprintf("  corr(whiff_resid, axis_diff) = %+.3f\n\n",
            cov.wt(valid[,.(whiff_resid, axis_diff)], wt=valid$n, cor=TRUE)$cor[1,2]))

show_cols <- function(d) d[, .(player_name, p_throws, pitch_type, sec_n=n, fb_type,
  velo_diff=round(velo_diff,1), axis_diff=round(axis_diff,1), spin_diff=round(spin_diff),
  vert_sep=round(vert_sep,1), whiff_pct=round(whiff_pct,1), csw_pct=round(csw_pct,1),
  whiff_resid=round(whiff_resid,1))]

cat("=== Top MATCHED-SPIN look-alikes (axis<=25 deg), by whiff-above-stuff ===\n")
print(show_cols(matched[velo_diff > 6][order(-whiff_resid)][1:20]))

cat("\n=== Top MIRROR pairs (axis>=155 deg), by whiff-above-stuff ===\n")
print(show_cols(mirror[order(-whiff_resid)][1:20]))

# ---- Changeup example set (Luzardo removed, Vesia added) -------------------

examples <- c("Cease, Dylan","Ribalta, Orlando","Dion, Will","Cantillo, Joey",
              "Ragans, Cole","Vesia, Alex","Boyd, Matthew","Bibee, Tanner",
              "Weathers, Ryan","Detmers, Reid","Cecconi, Slade","Sullivan, Sean")
ex <- pairs[player_name %in% examples & pitch_type == "CH"]
setorder(ex, -whiff_resid)
cat("\n=== Changeup examples (CH vs primary FB) ===\n")
print(ex[, .(player_name, p_throws, ch_n=n, velo_diff=round(velo_diff,1),
  axis_diff=round(axis_diff,1), spin_diff=round(spin_diff), vert_sep=round(vert_sep,1),
  whiff_pct=round(whiff_pct,1), whiff_resid=round(whiff_resid,1),
  wh_attack=round(wh_attack,1), wh_tilt=round(wh_tilt,1))])

# ---- Glasnow (mirror exemplar) --------------------------------------------

cat("\n=== Glasnow arsenal vs his primary FB (mirror exemplar) ===\n")
print(pairs[grepl("Glasnow", player_name),
  .(pitch_type, sec_n=n, fb_type, velo_diff=round(velo_diff,1),
    axis_diff=round(axis_diff,1), vert_sep=round(vert_sep,1),
    whiff_pct=round(whiff_pct,1), whiff_resid=round(whiff_resid,1))][order(-axis_diff)])

fwrite(pairs, file.path(out_dir, "pitch_pair_deception_2026.csv"))
cat(sprintf("\nWrote %s\n", file.path(out_dir, "pitch_pair_deception_2026.csv")))

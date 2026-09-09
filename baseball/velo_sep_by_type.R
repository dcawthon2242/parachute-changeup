#!/usr/bin/env Rscript
# Velocity separation by pitch type, from raw per-pitch Statcast (no bat tracking needed).
# Cell = pitcher x season x pitch type x batter side. Within cell-group (type, side, season)
# demeaning; outcomes regressed on the pitch's own velo holding the pitcher's FB velo,
# location and movement fixed. Reported per 1 mph SLOWER (= +1 mph separation).
# Usage: Rscript baseball/velo_sep_by_type.R 2020 2021 ... 2026
suppressPackageStartupMessages(library(data.table))
options(width = 220)
args <- commandArgs(trailingOnly = TRUE)
yrs <- if (length(args)) as.integer(args) else 2020:2026
MIN_SW <- 50; MIN_BIP <- 25

cols <- c("game_year","game_type","pitcher","p_throws","stand","pitch_type","description",
          "release_speed","pfx_x","pfx_z","plate_x","plate_z","sz_top","sz_bot",
          "launch_speed","estimated_woba_using_speedangle","delta_run_exp")
pp <- rbindlist(lapply(yrs, function(yr) {
  f <- sprintf("data/statcast_%d/statcast_%d_all.csv", yr, yr)
  if (!file.exists(f)) { cat("missing", f, "\n"); return(NULL) }
  x <- fread(f, select = cols, showProgress = FALSE)
  x[game_type == "R" & pitch_type != "" & is.finite(release_speed)]
}), fill = TRUE)
cat(sprintf("pitches %s, seasons %s\n", format(nrow(pp), big.mark=","), paste(sort(unique(pp$game_year)), collapse=",")))

pp[, ptype := pitch_type]
pp[ptype %in% c("KC","CS"), ptype := "CU"]
pp[, fam := fifelse(ptype %in% c("SL","ST","SV"), "SL+ST", ptype)]   # sweeper label only exists from 2023
WH <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SW <- c(WH, "foul", "hit_into_play")
pp <- pp[description %in% SW]                                        # swings only
pp[, whiff := description %in% WH]
pp[, bip := description == "hit_into_play"]
pp[, pfx_x_bat := fifelse(stand == "R", -pfx_x, pfx_x)]
pp[, plate_x_bat := fifelse(stand == "R", -plate_x, plate_x)]
pp[, pz_rel := (plate_z - sz_bot) / pmax(sz_top - sz_bot, 0.1)]
pp[, pit_rv := -delta_run_exp]

fb <- pp[pitch_type %in% c("FF","SI"), .(fb_velo = mean(release_speed)), by = .(pitcher, game_year)]

build <- function(key) {
  a <- pp[, .(n_swings = .N, whiff = mean(whiff), n_bip = sum(bip),
              xwobacon = mean(estimated_woba_using_speedangle[bip], na.rm = TRUE),
              ev = mean(launch_speed[bip], na.rm = TRUE),
              rv_swing = mean(pit_rv, na.rm = TRUE),
              velo = mean(release_speed), hb = mean(pfx_x_bat, na.rm=TRUE), ivb = mean(pfx_z, na.rm=TRUE),
              loc_x = mean(plate_x_bat, na.rm=TRUE), loc_z = mean(pz_rel, na.rm=TRUE)),
          by = c("pitcher", "game_year", key, "stand")]
  setnames(a, key, "pt")
  a <- merge(a, fb, by = c("pitcher","game_year"))
  a <- a[n_swings >= MIN_SW & is.finite(fb_velo)]
  a[, velo_kill := fb_velo - velo]
  dem <- function(v) v - mean(v, na.rm = TRUE)
  for (v in c("velo","fb_velo","velo_kill","loc_x","loc_z","hb","ivb","xwobacon","ev","whiff","rv_swing"))
    a[, (paste0("d_", v)) := dem(get(v)), by = .(pt, stand, game_year)]
  a
}

slope <- function(s, y, w) {
  s <- s[is.finite(get(y))]
  f <- lm(as.formula(paste0(y, " ~ 0 + d_velo + d_fb_velo + d_loc_x + d_loc_z + d_hb + d_ivb")), data = s, weights = s[[w]])
  co <- summary(f)$coefficients; c(-co["d_velo",1], -co["d_velo",3])
}
report <- function(a, types) {
  res <- rbindlist(lapply(types, function(t) {
    s <- a[pt == t]; sb <- s[n_bip >= MIN_BIP]
    xw <- slope(sb,"d_xwobacon","n_bip"); ev <- slope(sb,"d_ev","n_bip")
    wh <- slope(s,"d_whiff","n_swings"); rv <- slope(s,"d_rv_swing","n_swings")
    data.table(pitch = t, cells = nrow(s), pitchers = uniqueN(s$pitcher), swings = sum(s$n_swings),
               sep = round(weighted.mean(s$velo_kill, s$n_swings),1),
               xwOBACON = round(xw[1],4), t_xw = round(xw[2],1), EV = round(ev[1],2), t_ev = round(ev[2],1),
               whiff_pp = round(100*wh[1],2), t_wh = round(wh[2],1), RV = round(rv[1],5), t_rv = round(rv[2],1))
  }))
  print(res)
}
curve <- function(a, t, nb = 8) {
  s <- a[pt == t]
  for (y in c("d_whiff","d_xwobacon","d_rv_swing","d_ev")) {
    w <- if (y %in% c("d_xwobacon","d_ev")) s$n_bip else s$n_swings
    ok <- is.finite(s[[y]]) & w > 0; r <- rep(NA_real_, nrow(s))
    f <- lm(as.formula(paste0(y,"~0+d_fb_velo+d_loc_x+d_loc_z+d_hb+d_ivb")), data=s[ok], weights=w[ok])
    r[ok] <- resid(f); s[, (paste0("r",y)) := r]
  }
  s[, bin := cut(velo_kill, unique(quantile(velo_kill, seq(0,1,1/nb))), include.lowest=TRUE)]
  tab <- s[, .(sep_lo=round(min(velo_kill),1), sep_hi=round(max(velo_kill),1), cells=.N,
               whiff_pp=round(100*weighted.mean(rd_whiff,n_swings),2),
               xwOBACON=round(weighted.mean(rd_xwobacon,n_bip,na.rm=TRUE),4),
               EV=round(weighted.mean(rd_ev,n_bip,na.rm=TRUE),2),
               RV=round(weighted.mean(rd_rv_swing,n_swings),5)), by=bin][order(sep_lo)]
  s[, d_vk2 := d_velo_kill^2]
  q <- function(y, w) { ss <- s[is.finite(get(y)) & get(w) > 0]
    f <- lm(as.formula(paste0(y,"~0+d_velo_kill+d_vk2+d_fb_velo+d_loc_x+d_loc_z+d_hb+d_ivb")), data=ss, weights=ss[[w]])
    co <- summary(f)$coefficients; sprintf("lin %+.5f (t %+.1f) quad %+.6f (t %+.1f)", co[1,1],co[1,3],co[2,1],co[2,3]) }
  cat("\n---", t, "by separation bin (residualized) ---\n"); print(tab[, -"bin"])
  cat("  quad  whiff:", q("d_whiff","n_swings"), "\n        xwOBACON:", q("d_xwobacon","n_bip"), "\n        RV:", q("d_rv_swing","n_swings"), "\n")
}

cat("\n=== Pooled, slider family merged (SL+ST+SV), per 1 mph slower given same FB/location/movement ===\n")
A <- build("fam"); report(A, c("FF","SI","FC","CH","FS","SL+ST","CU"))
for (t in c("SL+ST","CH","CU","FS")) curve(A, t)
if (any(yrs >= 2023)) {
  cat("\n=== 2023+ only, SL and ST separate ===\n")
  B <- build("ptype")[game_year >= 2023]; report(B, c("SL","ST"))
  for (t in c("SL","ST")) curve(B, t)
}
fwrite(A, sprintf("data/swing_timing/velo_sep_cells_%d_%d.csv", min(yrs), max(yrs)))

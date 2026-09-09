#!/usr/bin/env Rscript
# Join OpenCommand (tomdoyo/open-command) per-pitch inferred targets to Statcast and build pitcher-season command metrics.
suppressPackageStartupMessages(library(data.table)); options(width=200)
out <- list(); val <- list()
for (yr in 2024:2026) {
  tg <- fread(cmd=sprintf("gzcat data/open_command/%d/targets.csv.gz", yr)); pb <- fread(cmd=sprintf("gzcat data/open_command/%d/pbp_info.csv.gz", yr), select=c("game_pk","play_id","game_type","pitcher_id","batter_id","pitch_type","vx0","vy0","vz0","plate_x","plate_z"))
  tg <- tg[status=="ok" & plausible==TRUE & is.finite(inferred_x_in) & is.finite(inferred_z_in)]
  oc <- merge(tg, pb, by=c("game_pk","play_id"))
  oc[, miss_in := sqrt((plate_x_in-inferred_x_in)^2 + (plate_z_in-inferred_z_in)^2)]
  oc[, miss_naive := sqrt((plate_x_in-naive_x_in)^2 + (plate_z_in-naive_z_in)^2)]
  oc[, `:=`(kx=round(vx0,2), ky=round(vy0,2), kz=round(vz0,2))]
  sc <- fread(sprintf("data/statcast_%d/statcast_%d_all.csv", yr, yr), select=c("game_pk","game_year","pitcher","batter","pitch_type","vx0","vy0","vz0","plate_x","plate_z","sz_top","sz_bot","description","balls","strikes","p_throws","stand","at_bat_number","pitch_number","game_type"), showProgress=FALSE)
  sc <- sc[game_type=="R" & is.finite(vx0)]; sc[, `:=`(kx=round(vx0,2), ky=round(vy0,2), kz=round(vz0,2))]
  sc[, dup := .N, by=.(game_pk,pitcher,batter,kx,ky,kz)]; sc <- sc[dup==1]
  J <- merge(sc, oc[, .(game_pk, pitcher=pitcher_id, batter=batter_id, kx, ky, kz, play_id, miss_in, miss_naive, target_x_in=inferred_x_in, target_z_in=inferred_z_in, plate_x_in, plate_z_in)], by=c("game_pk","pitcher","batter","kx","ky","kz"))
  cat(sprintf("%d: OpenCommand ok pitches %s | statcast regular-season pitches %s | joined %s (%.1f%% of OC)\n", yr, format(nrow(oc),big.mark=","), format(nrow(sc),big.mark=","), format(nrow(J),big.mark=","), 100*nrow(J)/nrow(oc)))
  cat(sprintf("   plate_x agreement check: mean |statcast plate_x*12 - OC plate_x_in| = %.2f in\n", mean(abs(J$plate_x*12 - J$plate_x_in))))
  out[[as.character(yr)]] <- J
  # validate vs their pitcher-level command_scores
  cs <- fread(sprintf("data/open_command/%d/command_scores.csv", yr))[pitch_type=="ALL"]
  nm <- fread(sprintf("data/statcast_%d/statcast_%d_all.csv", yr, yr), select=c("pitcher","player_name"), showProgress=FALSE)[, .SD[1], by=pitcher]
  mine <- J[, .(n=.N, med=median(miss_in)), by=pitcher]; mine <- merge(mine, nm, by="pitcher")
  # their names are "First Last"; statcast is "Last, First"
  mine[, nm2 := sapply(strsplit(player_name, ", "), function(z) if (length(z)==2) paste(z[2], z[1]) else z[1])]
  v <- merge(mine, cs[, .(nm2=pitcher, n_oc=n, inferred_in)], by="nm2"); v <- v[n>=200]
  cat(sprintf("   validation vs command_scores.csv (n>=200 pitches, %d pitchers): cor(median miss) = %.3f, mean abs diff = %.2f in, my n / their n = %.2f\n", nrow(v), cor(v$med, v$inferred_in), mean(abs(v$med-v$inferred_in)), mean(v$n/v$n_oc)))
}
J <- rbindlist(out); FB <- c("FF","SI","FC")
J[, pg2 := fifelse(pitch_type %in% FB, "FB", "OFF")]; J[, zone := abs(plate_x)<=0.83 & plate_z>=sz_bot & plate_z<=sz_top]
J[, tz := abs(target_x_in/12)<=0.83 & (target_z_in/12)>=sz_bot & (target_z_in/12)<=sz_top]   # was the TARGET in the zone
J[, pz_rel := (plate_z-sz_bot)/pmax(sz_top-sz_bot,0.1)]
C <- J[, .(n_cmd=.N, miss_med=median(miss_in), miss_mean=mean(miss_in), miss_fb=median(miss_in[pg2=="FB"]), miss_off=median(miss_in[pg2=="OFF"]),
           miss_naive_med=median(miss_naive), big_miss=mean(miss_in>18), tight=mean(miss_in<6),
           target_zone=mean(tz), target_edge=mean(tz & (abs(target_x_in/12)>0.55 | (target_z_in/12) < sz_bot+0.2*(sz_top-sz_bot) | (target_z_in/12) > sz_top-0.2*(sz_top-sz_bot))),
           hit_zone_when_targeted=mean(zone[tz]), miss_ahead=median(miss_in[strikes>balls]), miss_behind=median(miss_in[balls>strikes]), miss_two_k=median(miss_in[strikes==2])), by=.(pitcher, game_year)]
saveRDS(J[, .(game_pk, game_year, pitcher, batter, at_bat_number, pitch_number, pitch_type, pg2, miss_in, miss_naive, target_x_in, target_z_in, plate_x, plate_z, sz_top, sz_bot, balls, strikes, description)], "data/swing_timing/opencommand_pitches_2024_2026.rds")
fwrite(C, "data/swing_timing/command_pitchers_2024_2026.csv")
cat(sprintf("\npitcher-seasons %d | with >=500 tracked pitches %d\n", nrow(C), nrow(C[n_cmd>=500])))
a <- C[n_cmd>=500]; b <- copy(a)[, game_year := game_year-1L]; q <- merge(a,b,by=c("pitcher","game_year"),suffixes=c("","_n"))
cat("year-over-year reliability (>=500 tracked pitches both seasons, n=", nrow(q), "):\n")
for (v in c("miss_med","miss_mean","miss_fb","miss_off","miss_naive_med","big_miss","tight","target_zone","target_edge","hit_zone_when_targeted")) cat(sprintf("  %-24s r = %.3f\n", v, cor(q[[v]], q[[paste0(v,"_n")]], use="complete.obs")))

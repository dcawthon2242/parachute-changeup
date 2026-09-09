#!/usr/bin/env Rscript

# A SPEC BOTH LEAGUES CAN BUILD.
#
# The reason the imaged-axis spec does not port is now precise: D1's SpinAxis reproduces break
# direction at r = .999, and break direction is already inside both whiff models via ax_diff and
# az_diff. So D1 cannot supply a variable the model has not already seen through the movement
# channel. Searching harder among movement quantities is searching for something that cannot exist.
#
# What IS absent from BOTH models - the MLB one (whiff_parallel) and the D1 one (ncaa_13) - and
# present in both datasets:
#     arm angle          neither model carries it; D1 has a backcast version
#     spin efficiency    neither model carries it; D1 has an inferred version (r = .92 FF, .75 CH)
# Everything else on the D1 menu is either already a model feature or a deterministic function of
# one. So the portable spec has to be built from those two, alone or crossed with separation.
#
# The screen below runs identically in both leagues. That is the point: this is not an MLB search
# followed by a hoped-for replication, it is the same regression in both populations, reported
# together, so a feature that only works in one is visible immediately as such. Movement quantities
# are included as negative controls - they should be null in both, since both models price them.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 205); MDIR <- "data/statcast_model"
L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L)

## ---- MLB features ---------------------------------------------------------------------------------
P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(ax) & is.finite(az) & is.finite(ax_diff) & is.finite(az_diff)]
P[, `:=`(mx = fifelse(p_throws == "L", -ax, ax), mz = az + 32.174,
         fx = fifelse(p_throws == "L", -(ax - ax_diff), ax - ax_diff), fz = (az - az_diff) + 32.174)]
AM <- P[, .(n = .N, ivb = mean(mz)/32.174*12, hb = mean(mx)/32.174*12,
            fivb = mean(fz)/32.174*12, fhb = mean(fx)/32.174*12,
            spin = mean(release_spin_rate, na.rm = TRUE), velo = mean(release_speed),
            ext = mean(release_extension)), by = .(pitcher, season)][n >= 40]
AM[, id := as.character(pitcher)]
MLB <- merge(L[league == "MLB"], AM, by = c("id","season"))

## ---- D1 features ------------------------------------------------------------------------------------
SH <- readRDS(file.path(MDIR, "ncaa_pitcher_season_shapes.rds")); setDT(SH)
SH <- SH[pt %in% c("Changeup","Four-Seam") & throws %in% c("Right","Left") & n >= 30]
SH[, hbm := fifelse(throws == "Left", -hb, hb)]
WD <- dcast(SH, PitcherId + season ~ pt, value.var = c("ivb","hbm","spin","velo","ext"))
setnames(WD, gsub("-", "_", names(WD)))
WD <- WD[is.finite(ivb_Changeup) & is.finite(ivb_Four_Seam)]
WD[, `:=`(id = as.character(PitcherId), ivb = ivb_Changeup, hb = hbm_Changeup,
          fivb = ivb_Four_Seam, fhb = hbm_Four_Seam, spin = spin_Changeup,
          velo = velo_Changeup, ext = ext_Changeup)]
D1 <- merge(L[league == "D1" & !is.na(id) & id != "0" & id != ""],
            WD[, .(id, season, ivb, hb, fivb, fhb, spin, velo, ext)], by = c("id","season"))

## ---- shared derived quantities ------------------------------------------------------------------
prep <- function(D) {
  D <- copy(D)
  D[, `:=`(mv_mag = sqrt(ivb^2 + hb^2), d_ivb = ivb - fivb, d_hb = hb - fhb,
           mv_per_1k = sqrt(ivb^2 + hb^2)/(spin/1000),
           mv_axis = abs(((atan2(hb, ivb)*180/pi - atan2(fhb, fivb)*180/pi + 180) %% 360) - 180))]
  D[, hi := vs >= quantile(vs, 2/3)]
  D
}
MLB <- prep(MLB); D1 <- prep(D1)
cat(sprintf("MLB %d seasons / %d arms | D1 %d seasons / %d arms\n",
            nrow(MLB), uniqueN(MLB$id), nrow(D1), uniqueN(D1$id)))

crob <- function(D, f, k) {
  D <- D[complete.cases(D[, c(all.vars(f), "id"), with = FALSE])]
  m <- lm(f, D); u <- residuals(m); X <- model.matrix(m); nc <- uniqueN(D$id)
  if (!k %in% colnames(X)) return(c(NA,NA,NA))
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k]); c(e, s, 2*pt(-abs(e/s), nc-1))
}

## ---- 1. continuous screen, both leagues, inside high separation -------------------------------------
NOVEL <- c(arm_angle = "arm", ch_efficiency = "ec", ff_efficiency = "ef")
CTRL  <- c(mv_axis = "mv_axis", drop_vs_fb = "d_ivb", run_vs_fb = "d_hb",
           mv_magnitude = "mv_mag", mv_per_1k_spin = "mv_per_1k", spin_rate = "spin")
scan1 <- function(vars, tag) {
  cat(sprintf("\n  -- %s --\n", tag))
  R <- rbindlist(lapply(names(vars), function(nm) {
    row <- list(feature = nm)
    for (lg in c("MLB","D1")) {
      D <- copy(if (lg == "MLB") MLB else D1)[hi == TRUE]
      D[, z := as.numeric(scale(get(vars[[nm]])))]
      r <- crob(D, y ~ z + vs, "z")
      row[[paste0(lg, "_est")]] <- round(r[1],3); row[[paste0(lg, "_p")]] <- round(r[3],4)
    }
    as.data.table(row) }))
  R[, agree := sign(MLB_est) == sign(D1_est)]
  print(R, row.names = FALSE)
}
cat("\n=== 1. continuous screen inside the top separation tercile (per SD of the feature) ===")
scan1(NOVEL, "absent from both whiff models - the only genuine candidates")
scan1(CTRL,  "negative controls - already priced by both models, should be null")

## ---- 2. the same features crossed with separation ---------------------------------------------------
cat("\n=== 2. interaction with velocity separation, both leagues ===\n")
R2 <- rbindlist(lapply(names(NOVEL), function(nm) {
  row <- list(feature = nm)
  for (lg in c("MLB","D1")) {
    D <- copy(if (lg == "MLB") MLB else D1)
    D[, `:=`(z = as.numeric(scale(get(NOVEL[[nm]]))), zs = as.numeric(scale(vs)))]
    r <- crob(D, y ~ z * zs, "z:zs")
    row[[paste0(lg, "_est")]] <- round(r[1],3); row[[paste0(lg, "_p")]] <- round(r[3],4)
  }
  as.data.table(row) }))
print(R2, row.names = FALSE)

## ---- 3. gated specs, identical rules in both leagues ------------------------------------------------
cat("\n=== 3. candidate gates, same rule applied in both leagues ===\n")
cat("    thresholds are PERCENTILES of each league's own distribution, so the rule transfers\n\n")
gates <- list(
  "arm top third"                       = quote(arm >= q(arm, 2/3)),
  "arm top quartile"                    = quote(arm >= q(arm, .75)),
  "arm top decile"                      = quote(arm >= q(arm, .90)),
  "arm top third + CH eff top third"    = quote(arm >= q(arm, 2/3) & ec >= q(ec, 2/3)),
  "arm top third + CH eff bottom third" = quote(arm >= q(arm, 2/3) & ec <= q(ec, 1/3)),
  "CH eff top decile"                   = quote(ec >= q(ec, .90)),
  "arm top third + mv_axis bottom third"= quote(arm >= q(arm, 2/3) & mv_axis <= q(mv_axis, 1/3)))
G <- rbindlist(lapply(names(gates), function(k) {
  row <- list(gate = k)
  for (lg in c("MLB","D1")) {
    D <- copy(if (lg == "MLB") MLB else D1)
    q <- function(v, p) quantile(v, p, na.rm = TRUE)
    D[, g := eval(gates[[k]], D)]
    Dh <- D[hi == TRUE]
    r <- crob(Dh, y ~ g + vs, "gTRUE")
    row[[paste0(lg,"_n")]] <- sum(Dh$g, na.rm = TRUE)
    row[[paste0(lg,"_est")]] <- round(r[1],3); row[[paste0(lg,"_p")]] <- round(r[3],4)
  }
  as.data.table(row) }))
G[, agree := sign(MLB_est) == sign(D1_est)]
print(G, row.names = FALSE)
cat(sprintf("\n  %d of %d gates agree in sign across the two leagues\n", sum(G$agree, na.rm=TRUE), nrow(G)))

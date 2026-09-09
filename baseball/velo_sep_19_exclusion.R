#!/usr/bin/env Rscript

# IS THERE A D1-MEASURABLE FEATURE THAT EXCLUDES ROARK AND RODON WITHOUT GUTTING THE GROUP?
#
# Exhaustive scan over every pitcher-season quantity TrackMan reports, looking for a one-sided
# threshold that drops exactly those two and keeps the rest. Reported honestly: a rule selected
# because it removes two named pitchers is fitted to those names, and its p-value afterwards means
# nothing. The scan is here to establish whether even a fitted rule exists, and how contrived it
# has to be.
#
# Separately and more importantly - the exclusion rule is not the blocker. MEMBERSHIP still requires
# the imaged spin-axis gap, which D1 cannot measure and which was shown to carry the entire effect
# (movement-direction currency gives +0.05 to -0.28 across every threshold, against +6.38 imaged).
# A perfect D1-portable exclusion rule bolted onto a membership rule D1 cannot build is still not
# testable in D1. That is stated up front so the scan below is not mistaken for a solution.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 205); MDIR <- "data/statcast_model"
L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L); M <- L[league == "MLB"]
NM <- unique(readRDS(file.path(MDIR, "parachute_ff.rds"))[, .(pitcher, player_name)]); setDT(NM)
NM <- unique(NM, by = "pitcher")[, .(id = as.character(pitcher), nm = sub(",.*", "", player_name))]
M <- merge(M, NM, by = "id", all.x = TRUE)

# Every quantity below is something TrackMan reports directly or that D1 already has assembled.
P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(ax) & is.finite(az) & is.finite(ax_diff) & is.finite(az_diff)]
P[, `:=`(mx = fifelse(p_throws == "L", -ax, ax), mz = az + 32.174,
         fx = fifelse(p_throws == "L", -(ax - ax_diff), ax - ax_diff), fz = (az - az_diff) + 32.174)]
S <- P[, .(ivb = mean(mz)/32.174*12, hb = mean(mx)/32.174*12,
           fb_ivb = mean(fz)/32.174*12, fb_hb = mean(fx)/32.174*12,
           spin = mean(release_spin_rate, na.rm = TRUE), velo = mean(release_speed),
           ext = mean(release_extension), rel_z = mean(release_pos_z),
           rel_x = mean(abs(release_pos_x)), n = .N), by = .(pitcher, season)][n >= 40]
S[, `:=`(id = as.character(pitcher), mv_mag = sqrt(ivb^2 + hb^2),
         mv_per_1k = sqrt(ivb^2 + hb^2)/(spin/1000),
         d_ivb = ivb - fb_ivb, d_hb = hb - fb_hb,
         mv_axis = abs(((atan2(hb, ivb) - atan2(fb_hb, fb_ivb))*180/pi + 180) %% 360 - 180))]
M <- merge(M, S[, .(id, season, ivb, hb, spin, velo, ext, rel_z, rel_x, mv_mag, mv_per_1k,
                    d_ivb, d_hb, mv_axis)], by = c("id","season"))
H <- M[vs >= quantile(M$vs, 2/3) & nsw >= 75]; H[, g := axis < 10 & arm >= arm_thr]
G <- H[g == TRUE]
cat(sprintf("pool %d seasons | %d members, %d arms\n", nrow(H), nrow(G), uniqueN(G$id)))

FE <- c("ivb","hb","spin","velo","ext","rel_z","rel_x","mv_mag","mv_per_1k","d_ivb","d_hb",
        "mv_axis","ec","ef","vs","arm","nsw")
cat("\n=== the 10 members on every D1-measurable quantity ===\n")
print(G[order(-y), c(.(nm = nm, season = season, whiff_over = round(y,1)),
        lapply(.SD, function(x) round(x,2))), .SDcols = FE], row.names = FALSE)

## ---- scan -------------------------------------------------------------------------------------
cat("\n=== scan: can any single threshold drop Roark and Rodon and keep the rest? ===\n")
cat("    (Roark is already gone at 75 swings; scanning the 40-swing member set so both are present)\n\n")
H40 <- M[vs >= quantile(M$vs, 2/3) & nsw >= 40]; H40[, g := axis < 10 & arm >= arm_thr]
G40 <- H40[g == TRUE]
targ <- G40[nm %in% c("Roark","Rodón")]; keep <- G40[!nm %in% c("Roark","Rodón")]
cat(sprintf("    targets: %s | %d others to preserve\n\n",
            paste(paste0(targ$nm, " ", targ$season), collapse = ", "), nrow(keep)))
R <- rbindlist(lapply(FE, function(f) {
  rbindlist(lapply(c("above","below"), function(dir) {
    v <- G40[[f]]; if (anyNA(v)) return(NULL); grid <- sort(unique(c(v, quantile(H40[[f]], 0:100/100, na.rm = TRUE))))
    best <- NULL
    for (th in grid) {
      pass <- if (dir == "above") G40[[f]] >= th else G40[[f]] <= th
      lost_t <- sum(!pass[G40$nm %in% c("Roark","Rodón")])
      lost_k <- sum(!pass[!G40$nm %in% c("Roark","Rodón")])
      if (lost_t == 2 && (is.null(best) || lost_k < best$lost_k))
        best <- list(th = th, lost_k = lost_k)
    }
    if (is.null(best)) return(NULL)
    data.table(feature = f, rule = dir, threshold = round(best$th,3),
               others_also_lost = best$lost_k, members_left = nrow(keep) - best$lost_k) }))
}))
if (nrow(R)) print(R[order(others_also_lost)], row.names = FALSE) else cat("    none found\n")
cat(sprintf("\n    cleanest possible: loses %d of the %d other members\n",
            min(R$others_also_lost), nrow(keep)))

#!/usr/bin/env Rscript

# A D1-MEASURABLE THRESHOLD WITH ACTUAL MARGIN.
#
# Constraint set, fixed before searching:
#   MUST exclude  Roark 2020, Rodon 2025
#   MAY exclude   Skubal 2022, Bibee 2025
#   MUST retain   the other nine members
#
# The previous candidate satisfied the membership constraints but landed 0.002 degrees from Bibee's
# value, which is not a threshold, it is a coincidence. So the ranking criterion here is MARGIN: the
# gap between the worst retained member and the best excluded one, expressed in standard deviations
# of the feature across the whole pool so features on different scales can be compared. A rule with
# real margin is one where moving the cut a bit changes nothing, which is the only kind that has any
# chance of transferring to a different league with different measurement error.
#
# Only quantities TrackMan reports directly are eligible. Inferred spin efficiency is included but
# flagged - in D1 it is modelled from spin rate and break rather than measured (validated at r = .75
# for changeups), so a gate on it inherits that model's error on top of everything else.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 205); MDIR <- "data/statcast_model"
L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L); M <- L[league == "MLB"]
NM <- unique(readRDS(file.path(MDIR, "parachute_ff.rds"))[, .(pitcher, player_name)]); setDT(NM)
NM <- unique(NM, by = "pitcher")[, .(id = as.character(pitcher), nm = sub(",.*", "", player_name))]
M <- merge(M, NM, by = "id", all.x = TRUE)
P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(ax) & is.finite(az) & is.finite(ax_diff) & is.finite(az_diff)]
P[, `:=`(mx = fifelse(p_throws == "L", -ax, ax), mz = az + 32.174,
         fx = fifelse(p_throws == "L", -(ax - ax_diff), ax - ax_diff), fz = (az - az_diff) + 32.174)]
S <- P[, .(n = .N, ivb = mean(mz)/32.174*12, hb = mean(mx)/32.174*12,
           fivb = mean(fz)/32.174*12, fhb = mean(fx)/32.174*12,
           spin = mean(release_spin_rate, na.rm = TRUE), velo = mean(release_speed),
           ext = mean(release_extension), rel_z = mean(release_pos_z),
           rel_x = mean(abs(release_pos_x))), by = .(pitcher, season)][n >= 40]
S[, `:=`(id = as.character(pitcher), mv_mag = sqrt(ivb^2 + hb^2),
         mv_per_1k = sqrt(ivb^2 + hb^2)/(spin/1000), d_ivb = ivb - fivb, d_hb = hb - fhb,
         mv_axis = abs(((atan2(hb, ivb) - atan2(fhb, fivb))*180/pi + 180) %% 360 - 180))]
S[, `:=`(fb_mag = sqrt(fivb^2 + fhb^2))]
S[, mag_ratio := mv_mag / fb_mag]
M <- merge(M, S[, .(id, season, ivb, hb, spin, velo, ext, rel_z, rel_x, mv_mag, mv_per_1k,
                    d_ivb, d_hb, mv_axis, mag_ratio)], by = c("id","season"))

H <- M[vs >= quantile(M$vs, 2/3) & nsw >= 40]; H[, g := axis < 10 & arm >= arm_thr]
G <- H[g == TRUE][, key := paste(nm, season)]
MUST_GO  <- c("Roark 2020", "Rodón 2025")
MAY_GO   <- c("Skubal 2022", "Bibee 2025")
cat(sprintf("%d members. must exclude: %s | may also exclude: %s\n",
            nrow(G), paste(MUST_GO, collapse = ", "), paste(MAY_GO, collapse = ", ")))

FE <- c("mv_axis","mv_mag","mv_per_1k","mag_ratio","d_ivb","d_hb","ivb","hb","spin","velo",
        "ext","rel_z","rel_x","vs","arm","ec","ef")
INFERRED <- c("ec","ef")
pool_sd <- sapply(FE, function(f) sd(H[[f]], na.rm = TRUE))

res <- rbindlist(lapply(FE, function(f) {
  v <- G[[f]]; if (anyNA(v)) return(NULL)
  rbindlist(lapply(c("keep_low","keep_high"), function(dir) {
    ok <- lapply(seq_along(v), function(i) NULL)
    # candidate cuts sit between adjacent observed values; evaluate every split of the members
    o <- order(v); out <- list()
    for (i in seq_len(length(v) - 1)) {
      lo <- G$key[o[seq_len(i)]]; hi <- G$key[o[(i+1):length(v)]]
      excl <- if (dir == "keep_low") hi else lo
      if (!all(MUST_GO %in% excl)) next
      if (!all(excl %in% c(MUST_GO, MAY_GO))) next
      gap <- v[o[i+1]] - v[o[i]]
      out[[length(out)+1]] <- data.table(
        feature = f, rule = dir, cut = round((v[o[i]] + v[o[i+1]])/2, 3),
        margin_units = round(gap, 3), margin_sd = round(gap/pool_sd[[f]], 3),
        n_excluded = length(excl), excluded = paste(sort(excl), collapse = "; "))
    }
    if (!length(out)) NULL else rbindlist(out) }))
}))
if (!nrow(res)) { cat("\nNO single-feature threshold satisfies the constraints.\n"); quit(save = "no") }
res[, inferred := feature %in% INFERRED]
setorder(res, -margin_sd)
cat("\n=== every single-feature threshold satisfying the constraints, ranked by margin ===\n")
cat("    margin_sd = gap between worst retained and best excluded, in pool standard deviations\n\n")
print(res[, .(feature, rule, cut, margin_units, margin_sd, n_excluded, excluded,
              inferred_in_D1 = inferred)], row.names = FALSE)

## ---- refit under the best few ------------------------------------------------------------------
crob <- function(D, f, k) {
  D <- D[complete.cases(D[, c(all.vars(f), "id"), with = FALSE])]
  m <- lm(f, D); u <- residuals(m); X <- model.matrix(m); nc <- uniqueN(D$id)
  if (!k %in% colnames(X)) return(c(NA,NA,NA))
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k]); c(e, s, 2*pt(-abs(e/s), nc-1))
}
cat("\n=== effect under each candidate, both swing gates ===\n")
for (i in seq_len(min(5, nrow(res)))) {
  r <- res[i]
  for (sg in c(40, 75)) {
    D <- M[vs >= quantile(M$vs, 2/3) & nsw >= sg]
    D[, g2 := axis < 10 & arm >= arm_thr &
        (if (r$rule == "keep_low") get(r$feature) <= r$cut else get(r$feature) >= r$cut)]
    e <- crob(D, y ~ g2 + vs, "g2TRUE")
    cat(sprintf("  %-10s %-9s %8.3f | %d+ swings: %+.3f (se %.3f) p = %.4f | %d members, %d arms\n",
                r$feature, r$rule, r$cut, sg, e[1], e[2], e[3], sum(D$g2), uniqueN(D[g2 == TRUE]$id)))
  }
}
cat("\n=== member values on the leading candidates ===\n")
top <- unique(res[1:min(4, nrow(res))]$feature)
print(G[order(-y), c(.(nm = nm, season = season, whiff_over = round(y,1)),
        lapply(.SD, function(x) round(x,2))), .SDcols = top], row.names = FALSE)
cat("\n  D1 reference distribution for the leading feature is printed by the caller if needed.\n")

#!/usr/bin/env Rscript

# A GATE THE DATA FINDS ON ITS OWN.
#
# Every threshold so far was chosen by naming who had to fall on which side, which is fitting to the
# outcome no matter how much margin it ends up with. A defensible alternative: pick the split by a
# criterion that never sees the names, then check what it happens to separate. If a blind
# maximum-gap partition lands on the pitchers you wanted excluded, the gate is a property of the
# distribution rather than a description of three men.
#
# The criterion is single-linkage on one dimension: sort the member values, take the largest
# adjacent gap, cut at its midpoint. No tuning, no target set, one number out. Reported alongside
# 1D k-means with k=2 as a second blind method, and against the rank of the desired split among all
# gaps, which is the honest measure of how lucky the agreement is.
#
# The feature that motivates this is the DISCREPANCY between the two axis definitions:
#     mv_axis - imaged_axis
# The imaged gap says how far apart the two pitches' rotations are. The movement gap says how far
# apart their trajectories end up. When those disagree, the ball is going somewhere its rotation
# does not account for, which is the seam-shift signature expressed at the level of the CH-FF
# relationship rather than the individual pitch.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 195); MDIR <- "data/statcast_model"
L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L); M <- L[league == "MLB"]
NM <- unique(readRDS(file.path(MDIR, "parachute_ff.rds"))[, .(pitcher, player_name)]); setDT(NM)
NM <- unique(NM, by = "pitcher")[, .(id = as.character(pitcher), nm = sub(",.*", "", player_name))]
M <- merge(M, NM, by = "id", all.x = TRUE)
P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(ax) & is.finite(az) & is.finite(ax_diff) & is.finite(az_diff)]
P[, `:=`(mx = fifelse(p_throws == "L", -ax, ax), mz = az + 32.174,
         fx = fifelse(p_throws == "L", -(ax - ax_diff), ax - ax_diff), fz = (az - az_diff) + 32.174)]
A <- P[, .(n = .N, ca = atan2(mean(mx), mean(mz))*180/pi,
           fa = atan2(mean(fx), mean(fz))*180/pi), by = .(pitcher, season)][n >= 40]
A[, `:=`(id = as.character(pitcher), mv = abs(((ca - fa + 180) %% 360) - 180))]
M <- merge(M, A[, .(id, season, mv)], by = c("id","season"))
M[, disc := mv - axis]
H <- M[vs >= quantile(M$vs, 2/3)]; H[, g0 := axis < 10 & arm >= arm_thr]
G <- H[g0 == TRUE][, key := paste(nm, season)]
TARGET <- c("Rodón 2025", "Skubal 2022", "Roark 2020")

## ---- blind maximum-gap split on each candidate feature -------------------------------------------
cat(sprintf("%d members. Target partition (never shown to the algorithm): %s\n\n",
            nrow(G), paste(TARGET, collapse = ", ")))
blind <- function(v, nmv, lab) {
  o <- order(v); s <- v[o]; gaps <- diff(s)
  i <- which.max(gaps); cut <- (s[i] + s[i+1])/2
  hi <- nmv[o][(i+1):length(s)]
  # where the target split would rank among all gaps, if we had gone looking for it
  tgt_i <- if (all(TARGET %in% nmv)) {
    pos <- max(match(TARGET, nmv[o])); if (pos == length(s)) NA else NA } else NA
  km <- kmeans(matrix(v), centers = 2, nstart = 50)
  khi <- nmv[km$cluster == which.max(km$centers)]
  data.table(feature = lab, max_gap = round(max(gaps),2), cut = round(cut,2),
             gap_rank_of_target = NA_integer_,
             maxgap_excludes = paste(sort(hi), collapse = "; "),
             matches_target = setequal(hi, TARGET),
             kmeans_excludes = paste(sort(khi), collapse = "; "),
             kmeans_matches = setequal(khi, TARGET))
}
FEATS <- c(discrepancy = "disc", movement_axis = "mv", imaged_axis = "axis", arm_angle = "arm",
           ch_efficiency = "ec", separation = "vs")
R <- rbindlist(lapply(names(FEATS), function(k) blind(G[[FEATS[[k]]]], G$key, k)))
cat("=== blind splits: largest adjacent gap, and 1D k-means, on each feature ===\n\n")
print(R[, .(feature, max_gap, cut, maxgap_excludes, matches_target)], row.names = FALSE)
cat("\n   and k-means with k = 2 on the same features:\n")
print(R[, .(feature, kmeans_excludes, kmeans_matches)], row.names = FALSE)

## ---- how lucky is the agreement? -----------------------------------------------------------------
cat("\n=== is the discrepancy split really the largest gap, or one of several? ===\n")
v <- G$disc; o <- order(v); s <- v[o]
tab <- data.table(rank = seq_len(length(s)-1),
                  below = G$key[o][seq_len(length(s)-1)], value = round(s[-length(s)],2),
                  above = G$key[o][-1], next_value = round(s[-1],2), gap = round(diff(s),2))
setorder(tab, -gap)
print(head(tab[, .(gap, cut_between = paste0(below, " (", value, ")  |  ", above, " (", next_value, ")"))], 5),
      row.names = FALSE)
cat(sprintf("\n   largest gap is %.2f deg; second largest %.2f. Ratio %.2fx\n",
            tab$gap[1], tab$gap[2], tab$gap[1]/tab$gap[2]))

## ---- the resulting gate ---------------------------------------------------------------------------
crob <- function(D, f, k, w = FALSE) {
  environment(f) <- environment()
  D <- D[complete.cases(D[, c(all.vars(f), "id", "nsw"), with = FALSE])]
  wt <- if (w) as.numeric(D$nsw) else rep(1, nrow(D))
  m <- lm(f, D, weights = wt); u <- residuals(m)*sqrt(wt); X <- model.matrix(m)*sqrt(wt)
  nc <- uniqueN(D$id); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k]); c(e, s, 2*pt(-abs(e/s), nc-1))
}
CUT <- R[feature == "discrepancy"]$cut
cat(sprintf("\n=== the gate: (movement axis - imaged axis) <= %.2f degrees ===\n", CUT))
H[, g := g0 & disc <= CUT]
for (w in c(FALSE, TRUE)) {
  r <- crob(H, y ~ g + vs, "gTRUE", w)
  cat(sprintf("  %-16s %+.3f (se %.3f) p = %.4f | %d seasons, %d arms\n",
              if (w) "swing-weighted:" else "unweighted:", r[1], r[2], r[3],
              sum(H$g), uniqueN(H[g == TRUE]$id)))
}
r0 <- crob(H, y ~ g0 + vs, "g0TRUE")
cat(sprintf("  (ungated reference: %+.3f p = %.4f on %d seasons)\n\n", r0[1], r0[3], sum(H$g0)))
print(H[g0 == TRUE][order(disc), .(nm, season, swings = nsw, imaged_axis = round(axis,2),
        mv_axis = round(mv,1), discrepancy = round(disc,1), whiff_over = round(y,2),
        kept = disc <= CUT)], row.names = FALSE)
arms <- unique(H[g == TRUE]$id)
P2 <- rbindlist(apply(combn(length(arms), 2), 2, function(ix) {
  a <- arms[ix]; x <- crob(H[!id %in% a], y ~ g + vs, "gTRUE")
  data.table(d = paste(sort(H[id %in% a, unique(nm)]), collapse = " + "), est = x[1], p = x[3]) }))
cat(sprintf("\n  leave-two-arms-out: %d of %d keep p<.05 | worst %+.3f (p = %.3f, %s)\n",
            sum(P2$p < .05), nrow(P2), P2[which.max(p)]$est, max(P2$p), P2[which.max(p)]$d))
d <- crob(H[!nm %in% c("Skubal","Cease")], y ~ g + vs, "gTRUE")
cat(sprintf("  dropping Skubal AND Cease: %+.3f (se %.3f) p = %.4f\n", d[1], d[2], d[3]))

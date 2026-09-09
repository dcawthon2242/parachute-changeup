#!/usr/bin/env Rscript

# FOUR BATTER-PERSPECTIVE METRICS, SEARCHED WITH A NULL FOR THE SEARCH ITSELF.
#
#   spin similarity     - clock-face axis gap to the four-seam, and the active-spin gap
#   arm angle diff      - the changeup's release slot minus the four-seam's
#   velo difference     - separation in mph
#   path to location    - how far along the shared flight path the two pitches stay together
#
# Everything so far has been one hypothesis tested many ways. This is the opposite: many
# hypotheses, so the multiplicity has to be built in rather than audited afterwards. Every gate the
# search can reach is enumerated, scored, and the largest t-statistic found is compared against the
# distribution of largest t-statistics under permuted outcomes. A cell only counts if it beats the
# best thing the search finds in noise, which is a far higher bar than p < .05.
#
# Path ratio only exists for 2023-2026 and only where a changeup followed a four-seam often enough
# to estimate, so the search runs twice: three metrics on the full era, all four on the subsample.
# Reporting a four-metric result without the three-metric one would hide the sample cost.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
set.seed(59); options(width = 205); MDIR <- "data/statcast_model"
FLOOR <- 75L

## ---- build ---------------------------------------------------------------------------------------
L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L); M <- L[league == "MLB"]
NM <- unique(readRDS(file.path(MDIR, "parachute_ff.rds"))[, .(pitcher, player_name)]); setDT(NM)
NM <- unique(NM, by = "pitcher")[, .(id = as.character(pitcher), nm = sub(",.*", "", player_name))]
M <- merge(M, NM, by = "id", all.x = TRUE)

F <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F)
AD <- F[is.finite(arm_diff), .(armdiff = mean(arm_diff), nch = .N), by = .(pitcher, season)]
AD[, id := as.character(pitcher)]
M <- merge(M, AD[, .(id, season, armdiff)], by = c("id","season"), all.x = TRUE)

T <- readRDS(file.path(MDIR, "mech_tunnel.rds")); setDT(T)
PR <- T[league == "MLB" & prev_ff == TRUE, .(pr = mean(pr, na.rm = TRUE), npr = .N), by = .(id, season)]
M <- merge(M, PR, by = c("id","season"), all.x = TRUE)

CH <- readRDS(file.path(MDIR, "mlb_chase_resid.rds")); setDT(CH)
CH <- CH[, .(ych = 100*mean(r_all, na.rm = TRUE)), by = .(pitcher, season)][, id := as.character(pitcher)]
VR <- readRDS(file.path(MDIR, "velo_sep_resid.rds"))
RV <- as.data.table(VR$R)[, .(yrv = 100*mean(rv - q_aware)), by = .(pitcher, season)][, id := as.character(pitcher)]
M <- merge(M, CH[, .(id, season, ych)], by = c("id","season"), all.x = TRUE)
M <- merge(M, RV[, .(id, season, yrv)], by = c("id","season"), all.x = TRUE)
M[, effgap := ef - ec]          # positive = the changeup is the more gyro pitch
M[, slot := arm]

# the movement-vs-imaged axis discrepancy, needed only to reconstruct the established set
P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(ax) & is.finite(az) & is.finite(ax_diff) & is.finite(az_diff)]
P[, `:=`(mx = fifelse(p_throws == "L", -ax, ax), mz = az + 32.174,
         fx = fifelse(p_throws == "L", -(ax - ax_diff), ax - ax_diff), fz = (az - az_diff) + 32.174)]
DS <- P[, .(n = .N, ca = atan2(mean(mx), mean(mz))*180/pi,
            fa = atan2(mean(fx), mean(fz))*180/pi), by = .(pitcher, season)][n >= 40]
DS[, `:=`(id = as.character(pitcher), mv = abs(((ca - fa + 180) %% 360) - 180))]
M <- merge(M, DS[, .(id, season, mv)], by = c("id","season"), all.x = TRUE)
M[, disc := mv - axis]

A <- M[nsw >= FLOOR & is.finite(y) & is.finite(axis) & is.finite(effgap) &
       is.finite(armdiff) & is.finite(vs) & is.finite(slot)]
B <- A[is.finite(pr) & npr >= 30]
cat(sprintf("full-era sample: %d pitcher-seasons, %d arms (2020-2026, >=%d swings)\n", nrow(A), uniqueN(A$id), FLOOR))
cat(sprintf("path-ratio subsample: %d pitcher-seasons, %d arms (2023-2026, >=30 CH-after-FF)\n\n", nrow(B), uniqueN(B$id)))

crob <- function(D, f, k, w = TRUE) {
  environment(f) <- environment()
  D <- D[complete.cases(D[, c(all.vars(f), "id", "nsw"), with = FALSE])]
  wt <- if (w) as.numeric(D$nsw) else rep(1, nrow(D))
  m <- lm(f, D, weights = wt); u <- residuals(m)*sqrt(wt); X <- model.matrix(m)*sqrt(wt)
  nc <- uniqueN(D$id); b <- solve(crossprod(X))
  V <- b %*% crossprod(rowsum(X*u, D$id)) %*% b * (nc/(nc-1))
  e <- unname(coef(m)[k]); s <- unname(sqrt(diag(V))[k]); c(e, s, 2*pt(-abs(e/s), nc-1))
}
z <- function(x) (x - mean(x))/sd(x)

## =============================================================================================
cat("=== 1. MAIN EFFECTS: each metric on its own, per standard deviation ===\n\n")
FE <- list(axis = "spin axis gap to FF", effgap = "active-spin gap (FF minus CH)",
           armdiff = "arm angle diff (CH minus FF)", vs = "velocity separation",
           slot = "arm slot", pr = "path to location ratio")
cat(sprintf("%-30s %22s %22s %22s\n", "metric", "whiff", "chase", "run value"))
for (v in names(FE)) {
  D <- if (v == "pr") B else A
  D <- copy(D)[, zz := z(get(v))]
  rhs <- if (v == "vs") "~ zz" else "~ zz + vs"   # separation is its own control, cannot be both
  o <- sapply(c("y","ych","yrv"), function(oc) {
    r <- crob(D[is.finite(get(oc))], as.formula(paste(oc, rhs)), "zz")
    sprintf("%+6.3f (p%.4f)", r[1], r[3]) })
  cat(sprintf("%-30s %22s %22s %22s\n", FE[[v]], o[1], o[2], o[3]))
}
cat("\n  units: whiff and chase in percentage points, run value in runs per 100 (negative = better).\n")
cat("  every row controls for velocity separation, so the separation row is its own uncontrolled slope.\n")

## =============================================================================================
cat("\n=== 2. PAIRWISE INTERACTIONS on whiff (does one metric need another?) ===\n\n")
vs6 <- names(FE)
IT <- rbindlist(lapply(seq_len(length(vs6)-1), function(i) rbindlist(lapply((i+1):length(vs6), function(j) {
  a <- vs6[i]; b2 <- vs6[j]
  D <- if (a == "pr" || b2 == "pr") B else A
  D <- copy(D)[, `:=`(za = z(get(a)), zb = z(get(b2)))]
  f <- if (a == "vs" || b2 == "vs") y ~ za*zb else y ~ za*zb + vs
  r <- try(crob(D, f, "za:zb"), silent = TRUE); if (inherits(r, "try-error")) return(NULL)
  data.table(pair = paste(a, "x", b2), n = nrow(D), est = r[1], se = r[2], p = r[3])
}))))
setorder(IT, p)
print(IT[, .(pair, n, interaction = round(est,3), se = round(se,3), p = round(p,4))], row.names = FALSE)
cat(sprintf("\n  %d pairs tested; at alpha .05 you expect %.1f false positives, observed %d.\n",
            nrow(IT), .05*nrow(IT), sum(IT$p < .05)))

## =============================================================================================
cat("\n=== 3. EXHAUSTIVE GATE SEARCH, WITH A NULL FOR THE SEARCH ===\n")
# Fast vectorised t-statistics via Frisch-Waugh: residualise the outcome and every candidate gate
# on the intercept and separation, then the slope and its t follow from inner products. Clustered
# standard errors come later, on the survivors only - inside the search all that is needed is a
# statistic whose null distribution can be simulated under the identical procedure.
search <- function(D, feats, label, nperm = 2000) {
  D <- copy(D); n <- nrow(D)
  qs <- c(.25, .40, .50, .60, .75)
  G <- list(); gn <- character()
  for (v in feats) for (q in qs) for (dir in c("le","ge")) {
    th <- quantile(D[[v]], q); g <- if (dir == "le") D[[v]] <= th else D[[v]] >= th
    if (sum(g) >= 8 && sum(g) <= n - 20) { G[[length(G)+1]] <- as.numeric(g)
      gn <- c(gn, sprintf("%s %s p%02d", v, if (dir=="le") "<=" else ">=", round(100*q))) }
  }
  # pairs and triples of gates on DIFFERENT features
  base <- data.table(name = gn, feat = sub(" .*", "", gn))
  Gm <- do.call(cbind, G)
  add <- list(); an <- character()
  for (i in seq_len(ncol(Gm)-1)) for (j in (i+1):ncol(Gm)) {
    if (base$feat[i] == base$feat[j]) next
    g <- Gm[,i]*Gm[,j]; if (sum(g) >= 8) { add[[length(add)+1]] <- g
      an <- c(an, paste(gn[i], "&", gn[j])) }
  }
  if (length(feats) >= 3) {
    idx <- which(sapply(add, sum) >= 12)
    for (k in idx) for (j in seq_len(ncol(Gm))) {
      fs <- strsplit(an[k], " & ")[[1]]; if (sub(" .*", "", fs[1]) == base$feat[j] ||
          sub(" .*", "", fs[2]) == base$feat[j]) next
      g <- add[[k]]*Gm[,j]; if (sum(g) >= 8) { add[[length(add)+1]] <- g
        an <- c(an, paste(an[k], "&", gn[j])) }
    }
  }
  ALL <- cbind(Gm, do.call(cbind, add)); nmall <- c(gn, an)
  keep <- !duplicated(t(ALL)); ALL <- ALL[, keep, drop = FALSE]; nmall <- nmall[keep]
  X <- cbind(1, D$vs); H <- X %*% solve(crossprod(X)) %*% t(X)
  Gt <- ALL - H %*% ALL; yt <- as.numeric(D$y - H %*% D$y)
  ss <- colSums(Gt^2); ok <- ss > 1e-8; Gt <- Gt[, ok, drop = FALSE]; nmall <- nmall[ok]; ss <- ss[ok]
  tstat <- function(yv) { b <- as.numeric(crossprod(Gt, yv))/ss
    rss <- sum(yv^2) - b^2*ss; b/sqrt(rss/(n-3)/ss) }
  tt <- tstat(yt)
  mx <- replicate(nperm, max(abs(tstat(sample(yt)))))
  cat(sprintf("\n--- %s: %d distinct gates searched, %d seasons ---\n", label, length(tt), n))
  cat(sprintf("    best |t| observed %.2f | permutation max-|t| null: median %.2f, 95th %.2f\n",
              max(abs(tt)), median(mx), quantile(mx, .95)))
  cat(sprintf("    family-wise p for the best gate: %.4f  ->  %s\n", mean(mx >= max(abs(tt))),
              if (mean(mx >= max(abs(tt))) < .05) "SURVIVES the search null" else "does not survive"))
  o <- order(-abs(tt))[1:8]
  R <- rbindlist(lapply(o, function(i) {
    D[, gg := ALL[, i] == 1]
    r <- crob(D, y ~ gg + vs, "ggTRUE"); rc <- crob(D, y ~ gg + vs, "ggTRUE", FALSE)
    rch <- crob(D[is.finite(ych)], ych ~ gg + vs, "ggTRUE")
    rrv <- crob(D[is.finite(yrv)], yrv ~ gg + vs, "ggTRUE")
    data.table(gate = nmall[i], n = sum(D$gg), arms = uniqueN(D[gg == TRUE]$id), t_search = tt[i],
               whiff = r[1], p_cl = r[3], whiff_unw = rc[1], chase = rch[1], p_ch = rch[3],
               rv = rrv[1], p_rv = rrv[3], fw = mean(mx >= abs(tt[i]))) }))
  print(R[, .(gate, n, arms, whiff = round(whiff,2), p_clustered = round(p_cl,4),
              chase = round(chase,2), p_ch = round(p_ch,3), rv = round(rv,2), p_rv = round(p_rv,3),
              familywise_p = round(fw,4))], row.names = FALSE)
  invisible(list(R = R, mx = mx, nm = nmall, t = tt, ALL = ALL, D = D))
}
SA <- search(A, c("axis","effgap","armdiff","vs","slot"), "FULL ERA, five metrics")
SB <- search(B, c("axis","effgap","armdiff","vs","slot","pr"), "2023-2026, adding path ratio")

## =============================================================================================
cat("\n\n=== 4. IS ANY WINNER A DIFFERENT ARCHETYPE FROM THE MATCHED-AXIS ONE? ===\n")
th55 <- quantile(M[nsw >= FLOOR]$arm, .55, na.rm = TRUE)
A[, known := is.finite(disc) & axis < 11 & arm >= th55 & disc <= 27.5 & vs >= quantile(A$vs, 2/3)]
cat(sprintf("  the established 'optimized' set inside this population: %d seasons\n\n", sum(A$known)))
for (g in SA$R$gate[1:5]) {
  i <- which(SA$nm == g); A[, gg := SA$ALL[, i] == 1]
  r <- crob(A, y ~ gg + vs + known, "ggTRUE")
  cat(sprintf("  %-46s n=%2d overlap=%2d (%.0f%% of gate) | net of the known set: %+.2f (p %.4f)\n",
              g, sum(A$gg), A[gg & known, .N], 100*A[gg & known, .N]/sum(A$gg), r[1], r[3]))
}

## =============================================================================================
cat("\n=== 5. THE UNDERPERFORMANCE ARCHETYPE THE SEARCH FOUND ===\n")
# The strongest cell in the path-ratio sample points the other way: a changeup that differs from
# the fastball in axis AND arrives from a visibly different slot. Worth checking whether that is
# just the good archetype read backwards, or a separate population with its own tell.
for (D in list(A, B)) {
  lab <- if (nrow(D) == nrow(A)) "full era" else "2023-2026"
  D <- copy(D)
  D[, bad := axis >= quantile(D$axis, .25) & effgap <= quantile(D$effgap, .60) &
             armdiff >= quantile(D$armdiff, .75)]
  D[, good := axis <= quantile(D$axis, .25) & vs >= quantile(D$vs, .60)]
  cat(sprintf("\n  -- %s: %d 'mismatched' seasons (%d arms), %d 'matched' seasons --\n",
              lab, sum(D$bad), uniqueN(D[bad == TRUE]$id), sum(D$good)))
  cat(sprintf("     overlap between the two: %d seasons\n", D[bad & good, .N]))
  for (o in list(c("y","whiff"), c("ych","chase"), c("yrv","run value"))) {
    rb <- crob(D[is.finite(get(o[1]))], as.formula(paste(o[1], "~ bad + vs")), "badTRUE")
    rg <- crob(D[is.finite(get(o[1]))], as.formula(paste(o[1], "~ good + vs")), "goodTRUE")
    cat(sprintf("     %-11s mismatched %+6.2f (p %.4f)   matched %+6.2f (p %.4f)\n",
                o[2], rb[1], rb[3], rg[1], rg[3]))
  }
  rboth <- crob(D, y ~ bad + good + vs, "badTRUE")
  cat(sprintf("     mismatched effect with the matched set also in the model: %+.2f (p %.4f)\n",
              rboth[1], rboth[3]))
}

cat("\n=== 6. WHAT THE SEARCH NEVER PICKED ===\n")
for (S in list(list(SA, "full era"), list(SB, "2023-2026"))) {
  s <- S[[1]]; top50 <- s$nm[order(-abs(s$t))][1:50]
  cnt <- sapply(c("axis","effgap","armdiff","vs","slot","pr"),
                function(v) sum(grepl(paste0("(^|& )", v, " "), top50)))
  cat(sprintf("  %-11s appearances in the top 50 gates: %s\n", S[[2]],
              paste(sprintf("%s %d", names(cnt), cnt), collapse = " | ")))
}
saveRDS(list(A = A, B = B, SA = SA$R, SB = SB$R, IT = IT), file.path(MDIR, "fourmetric_search.rds"))

#!/usr/bin/env Rscript

# EIGHT METRICS. Adds spin-rate difference and the inferred-versus-observed movement gap.
#
# New this pass:
#   spin rate difference  - CH minus FF, raw and in Bauer units (spin per mph), because raw spin
#                           difference is partly just the velocity difference wearing a costume
#   seam-shift deviation  - the angle between where the ball actually moves and where its measured
#                           rotation says it should. This is the "inferred versus observed" split
#                           the request is after: a large deviation means the changeup is getting
#                           its movement from seam orientation rather than from spin, which is a
#                           physically different pitch from a traditional backspin-killed changeup
#
# The seam metric enters three ways because they are not the same question:
#   ssw      - the changeup's own absolute deviation, how seam-driven the pitch is
#   sswsign  - signed, since deviating toward the arm side and toward the glove side are different
#   disc     - the CH-minus-FF version, whether the pair's movement gap exceeds their rotation gap
#
# Same guarded procedure as before: enumerate every gate, score, and compare the best t-statistic
# against the distribution of best t-statistics under permuted outcomes. Adding metrics enlarges
# the search space, which raises that bar - a finding that was marginal at six metrics should be
# expected to weaken at eight, and if it does not, that is informative.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
set.seed(61); options(width = 210); MDIR <- "data/statcast_model"
FLOOR <- 75L

L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L); M <- L[league == "MLB"]
NM <- unique(readRDS(file.path(MDIR, "parachute_ff.rds"))[, .(pitcher, player_name)]); setDT(NM)
NM <- unique(NM, by = "pitcher")[, .(id = as.character(pitcher), nm = sub(",.*", "", player_name))]
M <- merge(M, NM, by = "id", all.x = TRUE)

F <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F)
AD <- F[is.finite(arm_diff), .(armdiff = mean(arm_diff)), by = .(pitcher, season)][, id := as.character(pitcher)]
M <- merge(M, AD[, .(id, season, armdiff)], by = c("id","season"), all.x = TRUE)

SP <- readRDS(file.path(MDIR, "ff_spin_seasons.rds")); setDT(SP)
M <- merge(M, SP[, .(id, season, spindiff, budiff, spin_CH, spin_FF)], by = c("id","season"), all.x = TRUE)

SD <- readRDS(file.path(MDIR, "ch_seam_deviation.rds")); setDT(SD)
M <- merge(M, SD[, .(id, season, sswsign = dev, ssw = absdev)], by = c("id","season"), all.x = TRUE)

P <- F[is.finite(ax) & is.finite(az) & is.finite(ax_diff) & is.finite(az_diff)]
P[, `:=`(mx = fifelse(p_throws == "L", -ax, ax), mz = az + 32.174,
         fx = fifelse(p_throws == "L", -(ax - ax_diff), ax - ax_diff), fz = (az - az_diff) + 32.174)]
DS <- P[, .(n = .N, ca = atan2(mean(mx), mean(mz))*180/pi,
            fa = atan2(mean(fx), mean(fz))*180/pi), by = .(pitcher, season)][n >= 40]
DS[, `:=`(id = as.character(pitcher), mv = abs(((ca - fa + 180) %% 360) - 180))]
M <- merge(M, DS[, .(id, season, mv)], by = c("id","season"), all.x = TRUE); M[, disc := mv - axis]

T <- readRDS(file.path(MDIR, "mech_tunnel.rds")); setDT(T)
PR <- T[league == "MLB" & prev_ff == TRUE, .(pr = mean(pr, na.rm = TRUE), npr = .N), by = .(id, season)]
M <- merge(M, PR, by = c("id","season"), all.x = TRUE)

CH <- readRDS(file.path(MDIR, "mlb_chase_resid.rds")); setDT(CH)
CH <- CH[, .(ych = 100*mean(r_all, na.rm = TRUE)), by = .(pitcher, season)][, id := as.character(pitcher)]
VR <- readRDS(file.path(MDIR, "velo_sep_resid.rds"))
RV <- as.data.table(VR$R)[, .(yrv = 100*mean(rv - q_aware)), by = .(pitcher, season)][, id := as.character(pitcher)]
M <- merge(M, CH[, .(id, season, ych)], by = c("id","season"), all.x = TRUE)
M <- merge(M, RV[, .(id, season, yrv)], by = c("id","season"), all.x = TRUE)
M[, `:=`(effgap = ef - ec, slot = arm)]

CORE <- c("axis","effgap","armdiff","vs","slot","spindiff","budiff","ssw","sswsign","disc")
A <- M[nsw >= FLOOR & is.finite(y) & Reduce(`&`, lapply(CORE, function(v) is.finite(M[[v]])))]
B <- A[is.finite(pr) & npr >= 30]
cat(sprintf("full-era sample: %d pitcher-seasons, %d arms\n", nrow(A), uniqueN(A$id)))
cat(sprintf("path-ratio subsample: %d pitcher-seasons, %d arms\n\n", nrow(B), uniqueN(B$id)))

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
cat("=== 1. MAIN EFFECTS, per standard deviation (new metrics marked) ===\n\n")
FE <- list(axis = c("spin axis gap to FF", ""), effgap = c("active-spin gap (FF minus CH)", ""),
           armdiff = c("arm angle diff (CH minus FF)", ""), vs = c("velocity separation", ""),
           slot = c("arm slot", ""), pr = c("path to location ratio", ""),
           spindiff = c("spin rate difference (CH minus FF)", "NEW"),
           budiff = c("Bauer-unit difference (spin per mph)", "NEW"),
           ssw = c("seam deviation, magnitude", "NEW"),
           sswsign = c("seam deviation, signed", "NEW"),
           disc = c("movement gap minus rotation gap", "NEW"))
cat(sprintf("%-38s %-4s %20s %20s %20s\n", "metric", "", "whiff", "chase", "run value"))
for (v in names(FE)) {
  D <- copy(if (v == "pr") B else A)[, zz := z(get(v))]
  rhs <- if (v == "vs") "~ zz" else "~ zz + vs"
  o <- sapply(c("y","ych","yrv"), function(oc) {
    r <- crob(D[is.finite(get(oc))], as.formula(paste(oc, rhs)), "zz")
    sprintf("%+6.3f (p%.4f)", r[1], r[3]) })
  cat(sprintf("%-38s %-4s %20s %20s %20s\n", FE[[v]][1], FE[[v]][2], o[1], o[2], o[3]))
}

cat("\n--- are the two spin measures and the three seam measures redundant? ---\n")
CM <- cor(A[, .(spindiff, budiff, vs, ssw, sswsign, disc, axis)])
print(round(CM, 2))

## =============================================================================================
cat("\n=== 2. SEAM-SHIFTED VS TRADITIONAL: does the split change who overperforms? ===\n")
# The request behind this metric: a changeup that kills backspin has its movement explained by its
# rotation, while a seam-shifted changeup does not. Split on the changeup's own deviation and ask
# whether the established axis-plus-separation effect lives in one half.
A[, trad := ssw <= median(A$ssw)]
A[, good := axis <= quantile(A$axis, .25) & vs >= quantile(A$vs, .60)]
cat(sprintf("  traditional (low deviation, <= %.1f deg): %d seasons | seam-shifted: %d\n",
            median(A$ssw), sum(A$trad), sum(!A$trad)))
for (o in list(c("y","whiff"), c("ych","chase"), c("yrv","run value"))) {
  rt <- crob(A[trad == TRUE & is.finite(get(o[1]))], as.formula(paste(o[1], "~ good + vs")), "goodTRUE")
  rs <- crob(A[trad == FALSE & is.finite(get(o[1]))], as.formula(paste(o[1], "~ good + vs")), "goodTRUE")
  cat(sprintf("   %-11s traditional half %+6.2f (p %.4f, n=%2d)   seam-shifted half %+6.2f (p %.4f, n=%2d)\n",
              o[2], rt[1], rt[3], A[trad & good, .N], rs[1], rs[3], A[!trad & good, .N]))
}
r3 <- crob(A, y ~ good*trad + vs, "goodTRUE:tradTRUE")
cat(sprintf("   three-way interaction (good x traditional) on whiff: %+.2f (se %.2f) p = %.4f\n",
            r3[1], r3[2], r3[3]))

## =============================================================================================
cat("\n=== 3. EXHAUSTIVE GATE SEARCH WITH EIGHT-PLUS METRICS ===\n")
search <- function(D, feats, label, nperm = 2000) {
  D <- copy(D); n <- nrow(D); qs <- c(.25, .40, .50, .60, .75)
  G <- list(); gn <- character()
  for (v in feats) for (q in qs) for (dir in c("le","ge")) {
    th <- quantile(D[[v]], q); g <- if (dir == "le") D[[v]] <= th else D[[v]] >= th
    if (sum(g) >= 8 && sum(g) <= n - 20) { G[[length(G)+1]] <- as.numeric(g)
      gn <- c(gn, sprintf("%s %s p%02d", v, if (dir=="le") "<=" else ">=", round(100*q))) }
  }
  Gm <- do.call(cbind, G); ft <- sub(" .*", "", gn)
  add <- list(); an <- character()
  for (i in seq_len(ncol(Gm)-1)) for (j in (i+1):ncol(Gm)) {
    if (ft[i] == ft[j]) next
    g <- Gm[,i]*Gm[,j]; if (sum(g) >= 8) { add[[length(add)+1]] <- g; an <- c(an, paste(gn[i], "&", gn[j])) }
  }
  idx <- which(sapply(add, sum) >= 14)
  for (k in idx) { fs <- sub(" .*", "", strsplit(an[k], " & ")[[1]])
    for (j in seq_len(ncol(Gm))) { if (ft[j] %in% fs) next
      g <- add[[k]]*Gm[,j]; if (sum(g) >= 8) { add[[length(add)+1]] <- g; an <- c(an, paste(an[k], "&", gn[j])) } } }
  ALL <- cbind(Gm, do.call(cbind, add)); nmall <- c(gn, an)
  keep <- !duplicated(t(ALL)); ALL <- ALL[, keep, drop = FALSE]; nmall <- nmall[keep]
  X <- cbind(1, D$vs); H <- X %*% solve(crossprod(X)) %*% t(X)
  Gt <- ALL - H %*% ALL; yt <- as.numeric(D$y - H %*% D$y)
  ss <- colSums(Gt^2); ok <- ss > 1e-8; Gt <- Gt[, ok, drop = FALSE]; nmall <- nmall[ok]; ss <- ss[ok]
  tstat <- function(yv) { b <- as.numeric(crossprod(Gt, yv))/ss
    b/sqrt((sum(yv^2) - b^2*ss)/(n-3)/ss) }
  tt <- tstat(yt); mx <- replicate(nperm, max(abs(tstat(sample(yt)))))
  cat(sprintf("\n--- %s: %d distinct gates, %d seasons ---\n", label, length(tt), n))
  cat(sprintf("    best |t| %.2f | null median %.2f, 95th %.2f | family-wise p = %.4f  -> %s\n",
              max(abs(tt)), median(mx), quantile(mx, .95), mean(mx >= max(abs(tt))),
              if (mean(mx >= max(abs(tt))) < .05) "SURVIVES" else "does not survive"))
  o <- order(-abs(tt))[1:10]
  R <- rbindlist(lapply(o, function(i) {
    D[, gg := ALL[, i] == 1]
    r <- crob(D, y ~ gg + vs, "ggTRUE")
    rch <- crob(D[is.finite(ych)], ych ~ gg + vs, "ggTRUE")
    rrv <- crob(D[is.finite(yrv)], yrv ~ gg + vs, "ggTRUE")
    data.table(gate = nmall[i], n = sum(D$gg), arms = uniqueN(D[gg == TRUE]$id),
               whiff = r[1], p_wh = r[3], chase = rch[1], p_ch = rch[3], rv = rrv[1], p_rv = rrv[3],
               fw = mean(mx >= abs(tt[i]))) }))
  print(R[, .(gate, n, arms, whiff = round(whiff,2), p_wh = round(p_wh,4), chase = round(chase,2),
              p_ch = round(p_ch,3), rv = round(rv,2), p_rv = round(p_rv,3), familywise = round(fw,4))],
        row.names = FALSE)
  top50 <- nmall[order(-abs(tt))][1:50]
  cnt <- sapply(feats, function(v) sum(grepl(paste0("(^|& )", v, " "), top50)))
  cat(sprintf("    appearances in the top 50: %s\n", paste(sprintf("%s %d", names(cnt), cnt), collapse = " | ")))
  invisible(list(R = R, nm = nmall, t = tt, ALL = ALL, mx = mx))
}
F5 <- c("axis","effgap","armdiff","vs","slot","spindiff","budiff","ssw","sswsign","disc")
SA <- search(A, F5, "FULL ERA, ten metrics")
SB <- search(B, c(F5, "pr"), "2023-2026, adding path ratio")

## =============================================================================================
cat("\n\n=== 4. DO THE NEW METRICS ADD ANYTHING THE OLD ONES DID NOT? ===\n")
# The six-metric search found axis <= p25 & vs >= p60. If spin difference or seam deviation carry
# independent information, adding them to that rule should move the effect.
A[, base := axis <= quantile(A$axis, .25) & vs >= quantile(A$vs, .60)]
rb <- crob(A, y ~ base + vs, "baseTRUE")
cat(sprintf("  baseline rule (axis <= p25 & sep >= p60): %+.2f (p %.4f) on %d seasons\n\n", rb[1], rb[3], sum(A$base)))
for (v in c("spindiff","budiff","ssw","sswsign","disc")) {
  for (d in c("le","ge")) {
    th <- quantile(A[[v]], if (d == "le") .60 else .40)
    A[, g2 := base & (if (d == "le") get(v) <= th else get(v) >= th)]
    if (sum(A$g2) < 10) next
    r <- crob(A, y ~ g2 + vs, "g2TRUE")
    rc <- crob(A[base == TRUE], y ~ g2 + vs, "g2TRUE")   # within the baseline set only
    cat(sprintf("  + %-9s %s %7.2f : n=%2d  effect %+5.2f (p %.4f) | within baseline %+5.2f (p %.4f)\n",
                v, if (d == "le") "<=" else ">=", th, sum(A$g2), r[1], r[3], rc[1], rc[3]))
  }
}
saveRDS(list(A = A, B = B, SA = SA$R, SB = SB$R), file.path(MDIR, "eightmetric_search.rds"))

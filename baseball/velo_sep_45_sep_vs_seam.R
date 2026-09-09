#!/usr/bin/env Rscript

# DOES VELOCITY SEPARATION BUY WHIFFS ON ITS OWN, OR ONLY THROUGH SEAM SHIFT?
#
# The article leaves this ambiguous. The traditional cell has the separation gate (8.7+ mph) and no
# seam shift and posts a whiff residual of +0.88 at p = .44, which reads like separation buying
# nothing. But separation also carries a significant +1.45 per SD in the gradient table. Both cannot
# be the headline, so this runs the horse race directly.
#
# Four questions:
#   1. Univariate, does each trait predict whiff above the model?
#   2. Head to head, does either survive the other?
#   3. Is the relationship multiplicative - does separation only pay when the seams are moving it?
#   4. Is "separation" even about the gap? Split it into the changeup's own speed and the fastball's
#      speed. If only the changeup term matters, this is "throw it slower" and the fastball anchor is
#      decorative. If the fastball term carries weight with changeup speed held fixed, the gap itself
#      is doing work and the effect is genuinely relational.
#
# The whiff residual is against a LightGBM model carrying the pitch's own velocity, spin, movement,
# extension and release point. It does NOT carry the fastball anchor, which is why separation can
# show up here at all.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 220); MDIR <- "data/statcast_model"

M <- readRDS(file.path(MDIR, "archetype_roster.rds")); setDT(M)
M <- M[league == "MLB"]

P <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(P)
P <- P[is.finite(release_speed) & is.finite(speed_diff)]
P[, id := as.character(pitcher)]
V <- P[, .(chvelo = mean(release_speed), ffvelo = mean(release_speed - speed_diff),
           nsw2 = sum(is_swing == 1, na.rm = TRUE),
           whiff_raw = 100*mean(whiff[is_swing == 1] == 1, na.rm = TRUE)), by = .(id, season)]
D <- merge(M, V, by = c("id","season"))
D <- D[is.finite(vs) & is.finite(ssw) & is.finite(axis) & is.finite(slot) & is.finite(y) &
       is.finite(chvelo) & is.finite(ffvelo)]
D <- D[nsw >= 75]
cat(sprintf("%d MLB pitcher-seasons, %d arms, 75-swing floor.\n\n", nrow(D), uniqueN(D$id)))

clus <- function(dat, f, k, wt = "nsw") {
  environment(f) <- environment()
  dat <- dat[complete.cases(dat[, c(all.vars(f), "id", wt), with = FALSE])]
  w <- as.numeric(dat[[wt]]); m <- lm(f, dat, weights = w)
  u <- residuals(m)*sqrt(w); X <- model.matrix(m)*sqrt(w); nc <- uniqueN(dat$id)
  keep <- !is.na(coef(m)); X <- X[, keep, drop = FALSE]
  b <- solve(crossprod(X)); V <- b %*% crossprod(rowsum(X*u, dat$id)) %*% b * (nc/(nc-1))
  j <- which(names(coef(m))[keep] == k)
  e <- unname(coef(m)[keep][j]); s <- unname(sqrt(diag(V))[j]); t <- qt(.975, nc-1)
  c(est = e, se = s, lo = e-t*s, hi = e+t*s, p = 2*pt(-abs(e/s), nc-1))
}
st <- function(r, d = 2) paste0(formatC(r["est"], width = d+5, digits = d, format = "f", flag = "+"),
  ifelse(r["p"] < .01, " **", ifelse(r["p"] < .05, " * ", "   ")))
ci <- function(r, d = 2) sprintf("[%+.*f,%+.*f]", d, r["lo"], d, r["hi"])

# sign so that "more" is always "more toward the profile"
D[, `:=`(z_sep = scale(vs)[,1], z_seam = scale(ssw)[,1], z_match = scale(-axis)[,1],
         z_slot = scale(slot)[,1], z_ch = scale(-chvelo)[,1], z_ff = scale(ffvelo)[,1])]
NICE <- c(z_sep = "Velocity separation", z_seam = "Seam deviation", z_match = "Axis match (less gap)",
          z_slot = "Arm slot")

## =============================================================================================
cat("=== 1. ONE AT A TIME, ON WHIFF ABOVE MODEL ===\n\n")
cat(sprintf("  %-24s %10s %18s %9s\n", "trait, per 1 SD", "whiff pts", "95% interval", "p"))
for (v in names(NICE)) {
  r <- clus(D, as.formula(paste("y ~", v)), v)
  cat(sprintf("  %-24s %10s %18s %9.4f\n", NICE[[v]], st(r), ci(r), r["p"]))
}

## =============================================================================================
cat("\n=== 2. HEAD TO HEAD ===\n\n")
SPEC <- list(
  "separation alone"              = y ~ z_sep,
  "seam deviation alone"          = y ~ z_seam,
  "both together"                 = y ~ z_sep + z_seam,
  "both + axis match"             = y ~ z_sep + z_seam + z_match,
  "both + axis match + arm slot"  = y ~ z_sep + z_seam + z_match + z_slot)
cat(sprintf("  %-30s %14s %14s %14s %14s\n", "specification", "separation", "seam dev",
            "axis match", "arm slot"))
for (nm in names(SPEC)) {
  f <- SPEC[[nm]]; vars <- all.vars(f)[-1]
  cells <- sapply(c("z_sep","z_seam","z_match","z_slot"), function(v)
    if (v %in% vars) st(clus(D, f, v)) else "       -      ")
  cat(sprintf("  %-30s %14s %14s %14s %14s\n", nm, cells[1], cells[2], cells[3], cells[4]))
}

## =============================================================================================
cat("\n=== 3. IS IT MULTIPLICATIVE? ===\n\n")
r <- clus(D, y ~ z_sep * z_seam, "z_sep:z_seam")
cat(sprintf("  separation x seam deviation interaction: %s  %s  p = %.4f\n\n",
            st(r), ci(r), r["p"]))
D[, `:=`(hs = vs >= median(vs), hd = ssw >= median(ssw))]
CELL <- D[, .(seasons = .N, arms = uniqueN(id), msep = mean(vs), mseam = mean(ssw),
              whiff_resid = weighted.mean(y, nsw), whiff_raw = weighted.mean(whiff_raw, nsw),
              chase_resid = weighted.mean(ych, nsw), rv_resid = weighted.mean(yrv, nsw)),
          by = .(separation = fifelse(hs, "high", "low"), seamdev = fifelse(hd, "high", "low"))]
setorder(CELL, -separation, -seamdev)
print(CELL[, .(separation, seamdev, seasons, arms, `mph` = round(msep,1), `deg` = round(mseam,1),
               `whiff above model` = round(whiff_resid,2), `raw whiff %` = round(whiff_raw,1),
               `chase above model` = round(chase_resid,2),
               `rv above model` = round(rv_resid,2))], row.names = FALSE)
cat("\n  separation effect measured separately inside each half of seam deviation:\n\n")
for (h in c(TRUE, FALSE)) {
  S <- copy(D[hd == h]); S[, z_sep := scale(vs)[,1]]
  r <- clus(S, y ~ z_sep, "z_sep")
  cat(sprintf("  %-28s %s  %s  p = %.4f  (n = %d)\n",
              paste(if (h) "high" else "low", "seam deviation"), st(r), ci(r), r["p"], nrow(S)))
}
cat("\n  and the seam effect inside each half of separation:\n\n")
for (h in c(TRUE, FALSE)) {
  S <- copy(D[hs == h]); S[, z_seam := scale(ssw)[,1]]
  r <- clus(S, y ~ z_seam, "z_seam")
  cat(sprintf("  %-28s %s  %s  p = %.4f  (n = %d)\n",
              paste(if (h) "high" else "low", "separation"), st(r), ci(r), r["p"], nrow(S)))
}

## =============================================================================================
cat("\n=== 4. IS 'SEPARATION' ACTUALLY ABOUT THE GAP? ===\n\n")
cat("  separation is fastball velocity minus changeup velocity. entering the two sides separately\n")
cat("  asks whether the gap matters or whether this is just 'throw the changeup slower'.\n")
cat("  a true gap effect needs the fastball term to carry weight with changeup speed held fixed.\n\n")
cat(sprintf("  %-34s %14s %14s\n", "specification", "slower CH", "harder FB"))
for (nm in c("velocity sides alone", "+ seam deviation", "+ seam deviation and axis match")) {
  f <- switch(nm, "velocity sides alone" = y ~ z_ch + z_ff,
                  "+ seam deviation" = y ~ z_ch + z_ff + z_seam,
                  y ~ z_ch + z_ff + z_seam + z_match)
  cat(sprintf("  %-34s %14s %14s\n", nm, st(clus(D, f, "z_ch")), st(clus(D, f, "z_ff"))))
}
r1 <- clus(D, y ~ z_ch + z_ff, "z_ch"); r2 <- clus(D, y ~ z_ch + z_ff, "z_ff")
cat(sprintf("\n  slower changeup  %s per SD  %s\n", st(r1), ci(r1)))
cat(sprintf("  harder fastball  %s per SD  %s\n", st(r2), ci(r2)))
cat(sprintf("  correlation between the two sides: r = %+.3f\n", D[, cor(-chvelo, ffvelo)]))

## =============================================================================================
cat("\n=== 5. SWINGS OR MISSES? ===\n\n")
# whiff rate is misses per swing. a pitch can raise it by drawing worse swings (chase) or by being
# harder to touch once swung at. these are different products and separation may only buy one.
cat(sprintf("  %-24s %14s %14s %14s\n", "trait, per 1 SD", "whiff above", "chase above", "rv above"))
for (v in c("z_sep","z_seam")) {
  cells <- sapply(c("y","ych","yrv"), function(o)
    st(clus(D, as.formula(paste(o, "~ z_sep + z_seam")), v)))
  cat(sprintf("  %-24s %14s %14s %14s\n", NICE[[v]], cells[1], cells[2], cells[3]))
}

## =============================================================================================
cat("\n=== 6. THE ARTICLE'S TWO CELLS, RECONCILED ===\n\n")
# the traditional archetype has the separation gate and no seam shift. if separation worked on its
# own that cell should show a whiff edge, and it does not. this checks whether that is because
# separation does nothing or because the cell is small.
for (tag in c("A7_trad","A5_seam","A4_broad")) {
  S <- copy(D); S[, g := get(tag)]
  r <- clus(S, y ~ g, "gTRUE")
  cat(sprintf("  %-32s %s  %s  p = %.4f  (n = %d seasons)\n",
              tag, st(r), ci(r), r["p"], S[g == TRUE, .N]))
}
cat("\n  inside the high-separation half only, what the seam split is worth:\n\n")
H <- D[vs >= 8.7]
cat(sprintf("  seasons at 8.7+ mph separation: %d\n", nrow(H)))
H[, sg := fifelse(ssw >= 7.8, "seam shifted", "not")]
print(H[, .(seasons = .N, arms = uniqueN(id), sep = round(mean(vs),1), seam = round(mean(ssw),1),
            `whiff above model` = round(weighted.mean(y, nsw),2),
            `raw whiff %` = round(weighted.mean(whiff_raw, nsw),1),
            `rv above model` = round(weighted.mean(yrv, nsw),2)), by = sg], row.names = FALSE)
r <- clus(copy(H)[, g := sg == "seam shifted"], y ~ g, "gTRUE")
cat(sprintf("\n  seam shifted vs not, within high separation: %s  %s  p = %.4f\n",
            st(r), ci(r), r["p"]))
## =============================================================================================
cat("\n=== 7. RAW WHIFF, BEFORE THE MODEL TAKES ITS CUT ===\n\n")
# section 1 found nothing above the model. that is not the same as finding nothing. the model prices
# the changeup's own velocity, and separation is mostly a function of that velocity, so the model
# absorbs the effect before the residual is formed. asking the raw question separates "no effect"
# from "already counted".
cat(sprintf("  %-24s %14s %14s\n", "trait, per 1 SD", "raw whiff %", "above model"))
for (v in c("z_sep","z_seam","z_match","z_slot")) {
  cat(sprintf("  %-24s %14s %14s\n", NICE[[v]],
              st(clus(D, as.formula(paste("whiff_raw ~", v)), v)),
              st(clus(D, as.formula(paste("y ~", v)), v))))
}
cat(sprintf("\n  correlation, separation with the changeup's own velocity: r = %+.3f\n",
            D[, cor(vs, chvelo)]))
cat("  that correlation is why the residual is quiet: the model already knows the pitch is slow.\n")
cat("\n  raw whiff by fifth of separation:\n\n")
D[, qs := cut(vs, quantile(vs, 0:5/5), include.lowest = TRUE, labels = paste0("Q", 1:5))]
print(D[, .(seasons = .N, sep = round(mean(vs),1), `CH velo` = round(mean(chvelo),1),
            `FB velo` = round(mean(ffvelo),1), `raw whiff %` = round(weighted.mean(whiff_raw, nsw),1),
            `above model` = round(weighted.mean(y, nsw),2)), by = qs][order(qs)], row.names = FALSE)

## =============================================================================================
cat("\n=== 8. WHY THE CELLS BEAT THE SLOPES ===\n\n")
# the archetype indicators run +3 to +4 while every continuous trait runs under +0.4 per SD. either
# the effect is a genuine conjunction that linear terms cannot express, or the gates were selected
# on the outcome and the indicator is carrying that selection. controlling the flag for the linear
# terms tells the two apart in the only direction the data can: if the flag survives, it is not
# merely a repackaging of the slopes.
cat(sprintf("  %-42s %14s\n", "archetype flag, whiff above model", "estimate"))
for (tag in c("A4_broad","A5_seam")) {
  S <- copy(D); S[, g := get(tag)]
  a <- clus(S, y ~ g, "gTRUE")
  b <- clus(S, y ~ g + z_sep + z_seam + z_match, "gTRUE")
  c2 <- clus(S, y ~ g + z_sep + z_seam + z_match + z_slot + z_ch + z_ff, "gTRUE")
  cat(sprintf("  %-42s %14s\n", paste(tag, "alone"), st(a)))
  cat(sprintf("  %-42s %14s\n", paste(tag, "+ the three linear traits"), st(b)))
  cat(sprintf("  %-42s %14s\n", paste(tag, "+ traits, slot and both velo sides"), st(c2)))
}
cat("\n  the flags survive the slopes, so the effect is conjunctive rather than linear.\n")
cat("  it is worth stating plainly that the gates were chosen from a threshold search on this same\n")
cat("  data, so part of that margin is selection and the honest estimate is below the printed one.\n")

cat("\n  whiff above model by fifth of each trait, to see where the effect actually sits:\n\n")
for (v in c("vs","ssw")) {
  D[, qq := cut(get(v), quantile(get(v), 0:5/5), include.lowest = TRUE, labels = paste0("Q", 1:5))]
  R <- D[, .(seasons = .N, mean = round(mean(get(v)),1),
             `above model` = round(weighted.mean(y, nsw),2)), by = qq][order(qq)]
  cat(sprintf("  %s: %s\n", ifelse(v == "vs", "separation    ", "seam deviation"),
              paste(sprintf("%s %+.2f", R$qq, R$`above model`), collapse = "   ")))
}
saveRDS(D, file.path(MDIR, "sep_vs_seam.rds"))

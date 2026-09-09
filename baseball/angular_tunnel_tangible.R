#!/usr/bin/env Rscript

# Can the angular tunnel metric be restated in units a pitching coach already uses --
# velocity separation, IVB separation, horizontal-break separation?
#
# Three questions, in order:
#
#   1 DESCRIPTIVE   how strongly does the break fraction track the three differentials,
#                   and how much of its variance do they jointly explain? If the answer is
#                   most of it, the metric is movement in disguise and can simply be
#                   restated. If the answer is some of it, the residual is the interesting
#                   part.
#   2 HORSE RACE    for predicting chase, does the break fraction survive controlling for
#                   the three differentials, and do the differentials survive controlling
#                   for it? Either can be redundant, both, or neither.
#   3 EXCHANGE RATE if it survives, express one SD of tunnel in mph and inches, so the
#                   effect has a size a coach can act on.
#
# Differentials are taken against the ACTUAL setup fastball on the previous pitch, not the
# pitcher's season average, because the tunnel metric is defined on that specific pair.
# Horizontal quantities are mirrored so positive is always toward the pitcher's arm side.

suppressPackageStartupMessages(library(data.table))
options(width = 215)
PRIMARY <- "brk_any_005"
GRP <- list(Slider = "SL", Sweeper = "ST", Curveball = c("CU","KC"),
            Changeup = "CH", Splitter = "FS")

d <- readRDS("data/statcast_model/angular_tunnel_2026.rds")
cf <- list.files("data/statcast_2026/chunks", pattern = "csv$", full.names = TRUE)
sh <- unique(rbindlist(lapply(cf, function(f) fread(f, select = c(
  "game_pk","at_bat_number","pitch_number","release_speed","pfx_x","pfx_z",
  "sz_top","sz_bot","game_type"), showProgress = FALSE)))[game_type == "R"],
  by = c("game_pk","at_bat_number","pitch_number"))
sh[, game_type := NULL]
lag <- sh[, .(game_pk, at_bat_number, pitch_number = pitch_number + 1L,
              p_speed = release_speed, p_pfx_x = pfx_x, p_pfx_z = pfx_z)]
K <- c("game_pk","at_bat_number","pitch_number")
d <- merge(d, sh, by = K, all.x = TRUE)
d <- merge(d, lag, by = K, all.x = TRUE)

mir <- function(x, L) fifelse(L, -x, x)
d[, L := p_throws == "L"]
# All three are "how much the secondary differs from the fastball that set it up".
d[, `:=`(d_velo = p_speed - release_speed,                       # mph slower
         d_ivb  = (p_pfx_z - pfx_z) * 12,                        # inches less rise
         d_hb   = (mir(p_pfx_x, L) - mir(pfx_x, L)) * 12)]       # inches, arm side positive
d[, d_mov := sqrt(d_ivb^2 + d_hb^2)]
d[, zdist := sqrt(pmax(abs(plate_x) - 0.95, 0)^2 + pmax(plate_z - sz_top, sz_bot - plate_z, 0)^2)]
d[, `:=`(out_zone = zdist > 0 & is.finite(zdist), chase = as.numeric(swing),
         count = paste(balls, strikes, sep = "-"))]
d[, pit_count := paste(pitcher, count)]
d[, grp := NA_character_]
for (g in names(GRP)) d[pitch_type %in% GRP[[g]], grp := g]
d <- d[!is.na(grp) & is.finite(d_velo) & is.finite(d_ivb) & is.finite(d_hb)]
d[, brk := get(PRIMARY)]

## ---- 1. descriptive -----------------------------------------------------
cat("=== 1. WHAT THE BREAK FRACTION IS MADE OF ===\n")
cat("Pearson r of the break fraction with each differential, then the joint R2 from\n")
cat("regressing the break fraction on all three (linear, and with squares/interactions\n")
cat("to allow curvature). d_velo in mph, d_ivb and d_hb in inches.\n\n")
desc <- rbindlist(lapply(names(GRP), function(g) {
  s <- d[grp == g]
  lin <- summary(lm(brk ~ d_velo + d_ivb + d_hb, s))$r.squared
  flx <- summary(lm(brk ~ poly(d_velo,2) + poly(d_ivb,2) + poly(d_hb,2) +
                      d_velo:d_ivb + d_velo:d_hb + d_ivb:d_hb, s))$r.squared
  data.table(group = g, n = nrow(s),
             mean_velo = mean(s$d_velo), mean_ivb = mean(s$d_ivb), mean_hb = mean(s$d_hb),
             r_velo = cor(s$brk, s$d_velo), r_ivb = cor(s$brk, s$d_ivb),
             r_hb = cor(s$brk, s$d_hb), r_absmov = cor(s$brk, s$d_mov),
             R2_linear = lin, R2_flexible = flx) }), fill = TRUE)
print(desc[, lapply(.SD, function(z) if (is.numeric(z)) round(z,3) else z)], row.names = FALSE)

cat("\nSo the three differentials explain between",
    sprintf("%.0f%% and %.0f%%", 100*min(desc$R2_flexible), 100*max(desc$R2_flexible)),
    "of the break fraction.\n")

## ---- 2. horse race ------------------------------------------------------
fe_fit <- function(dt, yvar, xvars, controls, fe) {
  vars <- c(yvar, xvars, controls)
  dt <- dt[stats::complete.cases(dt[, vars, with = FALSE])]
  if (nrow(dt) < 400) return(NULL)
  dm <- dt[, lapply(.SD, as.numeric), .SDcols = vars]
  dm[, `:=`(g2 = dt[[fe]], cl = dt$pitcher)]
  for (v in vars) dm[, (v) := get(v) - mean(get(v)), by = g2]
  X <- cbind(1, as.matrix(dm[, c(xvars, controls), with = FALSE])); y <- dm[[yvar]]
  XtXi <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XtXi)) return(NULL)
  b <- XtXi %*% crossprod(X, y); e <- as.vector(y - X %*% b)
  meat <- matrix(0, ncol(X), ncol(X))
  for (ix in split(seq_len(nrow(X)), dm$cl)) {
    u <- crossprod(X[ix, , drop = FALSE], e[ix]); meat <- meat + u %*% t(u) }
  V <- XtXi %*% meat %*% XtXi; nc <- uniqueN(dm$cl); V <- V*(nc/(nc-1))
  out <- data.table(term = c("(int)", xvars, controls), beta = as.vector(b),
                    se = sqrt(diag(V)))
  out[, `:=`(t = beta/se, p = 2*pt(-abs(beta/se), nc-1))]
  out[, sd_x := c(NA, sapply(c(xvars, controls), function(v) sd(dt[[v]])))]
  out[, n := nrow(dt)][, arms := nc][]
}

LOC <- c("plate_x","plate_z","zdist","release_speed","balls","strikes")
OZ <- d[out_zone == TRUE]
cat("\n\n=== 2. HORSE RACE FOR CHASE: TUNNEL vs THE THREE DIFFERENTIALS ===\n")
cat("Same fixed-effects spec as before. per_sd is the chase change in percentage points\n")
cat("for a 1 SD move in that variable, holding the others fixed.\n")
for (g in c("POOLED", names(GRP))) {
  s <- if (g == "POOLED") OZ else OZ[grp == g]
  a <- fe_fit(s, "chase", "brk", LOC, "pitcher")
  b <- fe_fit(s, "chase", c("d_velo","d_ivb","d_hb"), LOC, "pitcher")
  cc <- fe_fit(s, "chase", c("brk","d_velo","d_ivb","d_hb"), LOC, "pitcher")
  if (is.null(cc)) next
  cat(sprintf("\n--- %s (n=%s, base chase %.1f%%) ---\n", g,
              format(nrow(s), big.mark = ","), 100*mean(s$chase)))
  tab <- rbind(
    cbind(model = "tunnel alone",  a[term == "brk"]),
    cbind(model = "diffs alone",   b[term %in% c("d_velo","d_ivb","d_hb")]),
    cbind(model = "both together", cc[term %in% c("brk","d_velo","d_ivb","d_hb")]))
  print(tab[, .(model, term, beta = round(beta,5), per_sd_pp = round(100*beta*sd_x,2),
                t = round(t,2), p = signif(p,3))], row.names = FALSE)
}

## ---- 3. exchange rate ---------------------------------------------------
cat("\n\n=== 3. EXCHANGE RATE: ONE SD OF TUNNEL, IN MPH AND INCHES ===\n")
cat("From the joint model, the tunnel effect divided by each differential's per-unit\n")
cat("effect. Read as: 1 SD more tunnel buys the same chase as this much more separation.\n")
cat("Blank means that differential had the wrong sign or was too weak to convert against.\n\n")
xr <- rbindlist(lapply(c("POOLED", names(GRP)), function(g) {
  s <- if (g == "POOLED") OZ else OZ[grp == g]
  cc <- fe_fit(s, "chase", c("brk","d_velo","d_ivb","d_hb"), LOC, "pitcher")
  if (is.null(cc)) return(NULL)
  gb <- cc[term == "brk"]; eff <- gb$beta * gb$sd_x
  conv <- function(tm) { r <- cc[term == tm]
    if (!nrow(r) || !is.finite(r$beta) || abs(r$t) < 1.5) return(NA_real_)
    eff/r$beta }
  data.table(group = g, n = gb$n, tunnel_per_sd_pp = 100*eff,
             eq_mph = conv("d_velo"), eq_ivb_in = conv("d_ivb"), eq_hb_in = conv("d_hb"),
             sd_velo = cc[term=="d_velo"]$sd_x, sd_ivb = cc[term=="d_ivb"]$sd_x,
             sd_hb = cc[term=="d_hb"]$sd_x) }), fill = TRUE)
print(xr[, lapply(.SD, function(z) if (is.numeric(z)) round(z,2) else z)], row.names = FALSE)

cat("\n\n=== 4. WHICH DIFFERENTIAL IS THE TUNNEL ACTUALLY TRADING AGAINST? ===\n")
cat("Partial correlation of the break fraction with each differential, within pitcher,\n")
cat("after the other two are removed. Shows what a longer tunnel costs you in shape.\n\n")
pc <- rbindlist(lapply(names(GRP), function(g) {
  s <- d[grp == g]
  dm <- s[, .(brk, d_velo, d_ivb, d_hb, pitcher)]
  for (v in c("brk","d_velo","d_ivb","d_hb")) dm[, (v) := get(v) - mean(get(v)), by = pitcher]
  pr <- function(tm) { oth <- setdiff(c("d_velo","d_ivb","d_hb"), tm)
    ry <- residuals(lm(reformulate(oth, "brk"), dm))
    rx <- residuals(lm(reformulate(oth, tm), dm)); cor(ry, rx) }
  data.table(group = g, n = nrow(s), pr_velo = pr("d_velo"),
             pr_ivb = pr("d_ivb"), pr_hb = pr("d_hb")) }), fill = TRUE)
print(pc[, lapply(.SD, function(z) if (is.numeric(z)) round(z,3) else z)], row.names = FALSE)

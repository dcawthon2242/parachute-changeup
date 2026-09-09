#!/usr/bin/env Rscript

# IS A RELEASE-POINT SLOT PROXY GOOD ENOUGH TO EXTEND THE PARACHUTE BIN BACKWARDS?
#
# Reconstructing Savant's published arm angle failed at the precision this bin needs: the
# best per-pitch model reaches R2 = 0.81 and still misclassifies about a third of high-slot
# changeup seasons. That is because Savant's number uses Hawk-Eye body tracking of the
# shoulder, and release point cannot recover spine tilt.
#
# But exact reproduction is not the requirement. The requirement is a slot measure that can be
# computed identically in every season back to 2017 AND that preserves the effect. So: rebuild
# the bin inside 2023-2026 using the proxy instead of the published angle, and see whether the
# ground-ball result survives. If it does, the proxy is fit to extend. If it does not, no
# amount of extra seasons will help, because the extra seasons can only ever have the proxy.
#
# Three proxies are tested, cheapest first:
#   RAW    the fitted shoulder geometry, atan2 on release point and height
#   GBM    the boosted reconstruction of the published angle
#   RELZ   simple release height normalised by listed height, no angle at all

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(1); options(width = 200)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")
M <- readRDS(file.path(MDIR, "arm_angle_model2.rds"))

COLS <- c("pitcher","p_throws","pitch_type","arm_angle","release_pos_x","release_pos_y",
          "release_pos_z","release_extension")
d <- rbindlist(lapply(2023:2026, function(y) {
  k <- fread(sprintf("data/statcast_%d/statcast_%d_all.csv", y, y), select = COLS,
             showProgress = FALSE); k[, season := y]; k }))
d <- merge(d[pitch_type == "CH"], fread("data/pitcher_heights.csv")[, .(pitcher, ht_in)], by = "pitcher")
d <- d[is.finite(release_pos_x) & is.finite(release_pos_z) & is.finite(ht_in)]
d[, `:=`(ht_ft = ht_in/12, ax = fifelse(p_throws == "R", -release_pos_x, release_pos_x),
         thr = as.integer(p_throws == "R"))]
d[, raw := atan2(release_pos_z - M$par[1]*ht_ft, ax - M$par[2]*ht_ft)*180/pi]
ok <- stats::complete.cases(d[, M$feat, with = FALSE])
d[ok, gbmarm := predict(M$gbm, as.matrix(d[ok, M$feat, with = FALSE]))]
d[, relz_n := release_pos_z / ht_ft]

P <- d[, .(n = .N, arm_true = mean(arm_angle, na.rm = TRUE), raw = mean(raw, na.rm = TRUE),
           gbmarm = mean(gbmarm, na.rm = TRUE), relz_n = mean(relz_n, na.rm = TRUE)),
       by = .(pitcher, season)][n >= 60]

A <- fread(file.path(AST, "ext_parachute_filtered.csv"))
A <- merge(A, P[, .(pitcher, season, raw, gbmarm, relz_n)], by = c("pitcher","season"))
A <- A[is.finite(as_gap)]
cat(sprintf("%d pitcher-season changeups with the published angle, all three proxies,\nand measured active spin on both pitches\n\n", nrow(A)))

## ---- how well does each proxy rank the published angle? --------------------------
cat("=== proxy quality against the published arm angle ===\n")
for (v in c("raw","gbmarm","relz_n")) {
  s <- cor(A[[v]], A$arm, method = "spearman")
  cat(sprintf("  %-7s Spearman %.3f\n", v, s)) }

## ---- rebuild the bin with each proxy, matching the true bin's SIZE ----------------
# The published bin at arm >= 44 holds 11 seasons. Each proxy gets the threshold that
# selects the same count, so bin size is held fixed and only the membership changes.
TRUE_N <- sum(A$axis <= 10 & A$arm >= 44 & abs(A$as_gap) <= .10)
cat(sprintf("\n=== rebuilding the bin with each proxy, size held at n = %d ===\n", TRUE_N))
res <- rbindlist(lapply(c("arm","raw","gbmarm","relz_n"), function(v) {
  elig <- A$axis <= 10 & abs(A$as_gap) <= .10
  thr <- sort(A[[v]][elig], decreasing = TRUE)[min(TRUE_N, sum(elig))]
  i <- elig & A[[v]] >= thr
  tg <- t.test(A$gb[i], A$gb[!i]); tw <- t.test(A$wh[i], A$wh[!i])
  ov <- sum(i & A$axis <= 10 & A$arm >= 44 & abs(A$as_gap) <= .10)
  data.table(proxy = v, n = sum(i), overlap = ov,
             gb = diff(rev(tg$estimate)), p_gb = tg$p.value,
             wh = diff(rev(tw$estimate)), p_wh = tw$p.value,
             who = paste(sort(unique(sub(",.*","", A$player_name[i]))), collapse = " ")) }))
print(res[, .(proxy, n, shared_with_true = overlap, gb = round(gb,2), p_gb = round(p_gb,4),
              wh = round(wh,2), p_wh = round(p_wh,3))], row.names = FALSE)
cat("\nmembership under each proxy:\n")
for (i in seq_len(nrow(res))) cat(sprintf("  %-7s %s\n", res$proxy[i], res$who[i]))

## ---- the same question continuously ------------------------------------------------
cat("\n=== the axis x slot interaction, published angle vs each proxy ===\n")
for (v in c("arm","raw","gbmarm","relz_n")) {
  A[, sl := get(v)]
  for (y in c("wh","gb")) {
    co <- summary(lm(as.formula(sprintf("%s ~ axis * sl + kill + velo_sep + spin", y)), data = A))$coefficients
    cat(sprintf("  %-7s %-3s  axis:slot  beta = %+.5f   p = %.4f\n", v, y,
                co["axis:sl","Estimate"], co["axis:sl","Pr(>|t|)"])) } }

cat("\nVERDICT: compare the proxy rows to the 'arm' row. If the ground-ball effect and the\nmembership largely survive, older seasons are worth scraping; if they collapse, they are not.\n")
fwrite(A, file.path(AST, "ext_proxy_slot_feasibility.csv"))

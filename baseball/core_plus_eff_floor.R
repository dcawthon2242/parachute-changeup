#!/usr/bin/env Rscript

# Add an active-spin floor of .85 on BOTH pitches to the Core bin, and see what happens to the
# velocity-separation correlation.
#
# The Core rules already constrain the DISTANCE between the two efficiencies (|gap| <= .10) but
# put no floor on their LEVEL, so a pitcher throwing an 83 percent four-seamer and a 79 percent
# changeup satisfies the bin as comfortably as one throwing 96 and 100. Those are different
# pitches: at low active spin a large share of the rotation is gyro, the reported clock-face
# axis describes less of what the ball does, and "matched axis" starts to mean less.
#
# Both the correlation and the bin's mean residual are reported, because the floor changes who
# is in the bin and it is important to see whether the correlation moves for its own reasons or
# only because the membership changed.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
set.seed(3); options(width = 205)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

R  <- readRDS(file.path(MDIR, "whiff_tjstuff.rds"))
AS <- readRDS(file.path(MDIR, "active_spin_long.rds"))
S <- R[, .(nsw = .N, axis = mean(axis_diff), arm = mean(arm_angle), velo_sep = -mean(speed_diff),
           w = 100*mean(r4)), by = .(pitcher, player_name, season)]
S <- merge(S, AS[pitch_type == "CH", .(pitcher, season, as_ch = active_spin)], by = c("pitcher","season"))
S <- merge(S, AS[pitch_type == "FF", .(pitcher, season, as_fb = active_spin)], by = c("pitcher","season"))
S <- S[nsw >= 60]
S[, `:=`(name = trimws(paste(sub(".*,\\s*","",player_name), sub(",.*","",player_name))),
         as_gap = round(as_ch - as_fb, 4))]
S[, core := axis <= 10 & abs(as_gap) <= .10 & arm >= 44]
S[, corex := core & as_fb >= .85 & as_ch >= .85]
cat(sprintf("Core %d seasons -> Core + .85 floor %d seasons (%d dropped)\n\n",
            sum(S$core), sum(S$corex), sum(S$core) - sum(S$corex)))

cat("=== who the floor removes ===\n")
print(S[core & !corex][order(-w), .(Pitcher = name, Season = season, Swings = nsw,
      ActFF = round(as_fb,2), ActCH = round(as_ch,2), VeloSep = round(velo_sep,1),
      WhiffAbove = round(w,1))], row.names = FALSE)
cat(sprintf("\n  dropped seasons average %+.2f pp of whiff above model (%d of %d negative)\n",
            mean(S[core & !corex]$w), sum(S[core & !corex]$w < 0), sum(S$core & !S$corex)))

fmt <- function(D, lab) { ct <- cor.test(D$velo_sep, D$w)
  data.table(bin = lab, n = nrow(D), mean_resid = mean(D$w), r = ct$estimate,
             ci = sprintf("[%+.2f, %+.2f]", ct$conf.int[1], ct$conf.int[2]),
             slope = coef(lm(w ~ velo_sep, D))[2], p = ct$p.value,
             rho = cor(D$velo_sep, D$w, method = "spearman")) }
cat("\n=== velocity separation vs whiff residual ===\n")
print(rbind(fmt(S[core == TRUE],  "Core"),
            fmt(S[corex == TRUE], "Core + .85 floor"),
            fmt(S[core == FALSE], "everything else"))[, .(bin, n,
            mean_resid = round(mean_resid,2), r = round(r,3), ci, slope = round(slope,2),
            p = round(p,4), spearman = round(rho,3))], row.names = FALSE)

cat("\n=== leverage inside Core + .85 floor ===\n")
B <- S[corex == TRUE]
for (drop in list(character(0), "Dylan Cease", c("Dylan Cease","Tarik Skubal"))) {
  D <- B[!name %in% drop]; ct <- cor.test(D$velo_sep, D$w)
  cat(sprintf("  drop %-26s n=%2d  r=%+.3f  slope %+.2f  p=%.4f  spearman %+.3f\n",
      if (length(drop)) paste(drop, collapse=" + ") else "nothing", nrow(D), ct$estimate,
      coef(lm(w ~ velo_sep, D))[2], ct$p.value, cor(D$velo_sep, D$w, method="spearman"))) }

cat("\n=== interaction against the rest of the league ===\n")
for (b in c("core","corex")) { S[, bb := get(b)]
  for (wt in c(FALSE, TRUE)) {
    f <- if (wt) lm(w ~ velo_sep*bb, S, weights = nsw) else lm(w ~ velo_sep*bb, S)
    d <- summary(f)$coefficients["velo_sep:bbTRUE",]
    cat(sprintf("  %-6s %-10s %+.2f pp per mph, p = %.4f\n", b,
        if (wt) "weighted" else "unweighted", d[1], d[4])) } }

# Sweeping the floor separates "the floor is doing something" from "one threshold happened to
# land well". A real effect should move smoothly rather than spike at .85.
cat("\n=== sweeping the floor ===\n")
sweep <- rbindlist(lapply(c(0, .75, .80, .82, .85, .88, .90, .92), function(fl) {
  i <- S$core & S$as_fb >= fl & S$as_ch >= fl
  if (sum(i) < 6) return(NULL)
  D <- S[i]; ct <- cor.test(D$velo_sep, D$w)
  data.table(floor = fl, seasons = sum(i), mean_resid = mean(D$w), r = ct$estimate,
             p = ct$p.value, spearman = cor(D$velo_sep, D$w, method="spearman")) }))
print(sweep[, .(floor, seasons, mean_resid = round(mean_resid,2), r = round(r,3),
                p = round(p,4), spearman = round(spearman,3))], row.names = FALSE)

# The floor flattens the line, which raises the possibility that the Core bin's velocity
# correlation was partly an efficiency correlation wearing a disguise: if the low-efficiency
# seasons also happen to sit at low separation, they occupy the lower-left of the scatter and
# prop up a positive slope without velocity having anything to do with it.
cat("\n=== is the Core velocity correlation partly an efficiency correlation? ===\n")
C <- S[core == TRUE]; C[, as_min := pmin(as_fb, as_ch)]
for (p in list(c("as_min","w","lower of the two active spins vs whiff residual"),
               c("as_min","velo_sep","lower active spin vs velocity separation"),
               c("velo_sep","w","velocity separation vs whiff residual"))) {
  ct <- cor.test(C[[p[1]]], C[[p[2]]])
  cat(sprintf("  %-46s r = %+.3f  p = %.4f\n", p[3], ct$estimate, ct$p.value)) }
pr <- summary(lm(w ~ velo_sep + as_min, C))$coefficients
cat(sprintf("  both in one model: velo %+.2f (p=%.3f), active spin %+.1f (p=%.3f)\n",
            pr["velo_sep","Estimate"], pr["velo_sep","Pr(>|t|)"],
            pr["as_min","Estimate"], pr["as_min","Pr(>|t|)"]))

## ---- figure ---------------------------------------------------------------------------------
P <- rbind(cbind(S[core == TRUE], panel = sprintf("Core, as defined (n = %d)", sum(S$core))),
           cbind(S[corex == TRUE], panel = sprintf("Core + active spin >= .85 on both (n = %d)", sum(S$corex))))
P[, panel := factor(panel, levels = unique(panel))]
G <- S[, .(velo_sep, w)]
gg <- ggplot(P, aes(velo_sep, w)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_point(data = G, colour = "grey80", size = .7, alpha = .35) +
  geom_smooth(data = G, method = "lm", se = FALSE, colour = "grey45", linewidth = .6, formula = y ~ x) +
  geom_point(aes(size = nsw), colour = "#1d7870", alpha = .85) +
  geom_smooth(method = "lm", se = TRUE, colour = "#1d7870", fill = "#1d787033",
              linewidth = .9, formula = y ~ x) +
  geom_text_repel(aes(label = sprintf("%s '%02d", name, season %% 100)), size = 2.65,
                  seed = 3, min.segment.length = 0, max.overlaps = 20) +
  facet_wrap(~panel) +
  scale_size_continuous(range = c(1, 5.5), guide = "none") +
  labs(title = "The efficiency floor triples the bin's overperformance and weakens its velocity correlation. Those are two different findings",
       subtitle = paste0("Whiff above a shape, location and arm-angle model against velocity separation. Grey is all 1,086 changeup seasons and is flat at r = -0.01, the correct null\n",
                         "since separation is already a model feature. The Core rules cap the DISTANCE between the two active-spin figures at ten points but never require either to be\n",
                         "high, so Tanner Banks at 83 and 75 percent qualifies alongside Martin Perez at 96 and 100. Requiring .85 on both removes six seasons averaging -4.1 points,\n",
                         "including the bin's two worst, and mean overperformance goes from +0.48 to +2.33. The correlation moves the other way: +0.41 to +0.33, Spearman +0.34 to\n",
                         "+0.21, and the interaction against the league from p = .02 to p = .19. That is not a contradiction. Inside the Core bin, active-spin level is itself the\n",
                         "stronger predictor of overperformance (r = +0.50 against velocity's +0.41) and the two are almost independent of each other (r = +0.15), so both survive when\n",
                         "fitted together. Gating on one of two independent predictors truncates the outcome's range and flattens the other's slope - a mechanical consequence, not\n",
                         "evidence against velocity. What should give pause is the leverage: inside the floored bin, dropping Cease alone takes the correlation to +0.06."),
       x = "Velocity separation from the four-seamer (mph slower)",
       y = "Whiff above model (percentage points)",
       caption = "Source: Statcast 2020-2026 - four-seam anchor - point size is changeup swings - 60+ swings per season") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5), plot.subtitle = element_text(size = 8.2),
        strip.text = element_text(face = "bold", size = 9.5), panel.grid.minor = element_blank(),
        plot.caption = element_text(size = 7.5, colour = "grey40"))
ggsave(file.path(AST, "fig39_core_eff_floor.png"), gg, width = 12, height = 6.8, dpi = 150)
cat("\nwrote fig39_core_eff_floor.png\n")

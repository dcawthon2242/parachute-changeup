#!/usr/bin/env Rscript

# DOES THE WHIFF RESIDUAL STILL TRACK VELOCITY SEPARATION?
#
# It should not. speed_diff is one of the 24 features the shape-and-location model trains on,
# so a changeup 12 mph off the fastball already has a higher predicted whiff rate than one 5
# mph off, and the residual is what is left after that. Any surviving correlation is not a
# discovery about changeups - it is the model failing to absorb an input it was handed.
#
# Three levels, because they answer different questions:
#   raw       whiff rate itself against velo separation. The relationship everyone expects.
#   pitch     per-pitch residual against velo separation. The orthogonality check.
#   season    pitcher-season residual against mean velo separation. What the roster shows,
#             and where aggregation can reintroduce structure the pitch level does not have.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
options(width = 200)
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

ct <- function(x, y, lab) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  r <- cor.test(x, y); s <- suppressWarnings(cor(x, y, method = "spearman"))
  cat(sprintf("  %-42s n=%7d  r=%+.4f  [%+.3f,%+.3f]  rho=%+.3f  p=%.3g\n",
              lab, length(x), r$estimate, r$conf.int[1], r$conf.int[2], s, r$p.value))
  invisible(r)
}

## ---- per-pitch -------------------------------------------------------------------------
L  <- readRDS(file.path(MDIR, "parachute_extended_resid.rds"))
CH <- L$CH; SW <- L$SW
stopifnot(identical(CH[is_swing == TRUE]$pitcher, SW$pitcher),
          identical(CH[is_swing == TRUE]$season,  SW$season))
P <- cbind(CH[is_swing == TRUE, .(pitcher, season, speed_diff, axis_diff, arm_angle)],
           SW[, .(wh_res)])
P[, velo_sep := -speed_diff]     # mph SLOWER than the fastball; higher = more separation

RAW <- readRDS(file.path(MDIR, "parachute_extended.rds"))
RAW <- RAW[is_swing == TRUE & is.finite(speed_diff) & is.finite(whiff)]
RAW[, velo_sep := -speed_diff]

cat(sprintf("per-pitch swings: %d with residual, %d raw\n\n", nrow(P), nrow(RAW)))
cat("=== per-pitch ===\n")
ct(RAW$velo_sep, RAW$whiff,   "raw whiff vs velo separation")
ct(P$velo_sep,   P$wh_res,    "whiff RESIDUAL vs velo separation")

## ---- pitcher-season --------------------------------------------------------------------
S <- fread(file.path(AST, "ext_parachute_extended.csv"))
S <- S[is.finite(as_gap) & is.finite(arm)]
tier <- function(ax, ef, ar) S$axis <= ax & abs(S$as_gap) <= ef & S$arm >= ar
S[, wide := tier(15, .15, 42)][, core := tier(10, .10, 44)]

cat("\n=== pitcher-season (mean velo separation vs mean whiff residual) ===\n")
ct(S$velo_sep, S$wh, "all eligible changeup-seasons")
ct(S[wide == TRUE]$velo_sep, S[wide == TRUE]$wh, "inside the Wide parachute bin")
ct(S[core == TRUE]$velo_sep, S[core == TRUE]$wh, "inside the Core parachute bin")
ct(S[np >= 300]$velo_sep, S[np >= 300]$wh, "seasons with 300+ changeups (less noise)")

cat("\n=== the same slice against the other two residuals, for contrast ===\n")
ct(S$velo_sep, S$gb,    "grounder residual vs velo separation")
ct(S$velo_sep, S$rv100, "run-value residual vs velo separation")

## ---- binned curves ---------------------------------------------------------------------
bin <- function(D, v, y, n = 10) {
  D <- D[is.finite(get(v)) & is.finite(get(y))]
  D[, b := cut(get(v), quantile(get(v), seq(0, 1, length.out = n + 1)),
               include.lowest = TRUE, labels = FALSE)]
  D[, .(x = mean(get(v)), y = mean(get(y)), se = sd(get(y))/sqrt(.N), n = .N), by = b][order(b)]
}
A <- bin(RAW, "velo_sep", "whiff"); A[, `:=`(y = 100*y, se = 100*se, panel = "Raw whiff rate (%)")]
B <- bin(P,   "velo_sep", "wh_res"); B[, `:=`(y = 100*y, se = 100*se, panel = "Whiff residual above model (pp)")]
cat("\n=== per-pitch deciles of velo separation ===\n")
print(rbind(A, B)[, .(panel, decile = b, velo_sep = round(x,1), value = round(y,2),
                      se = round(se,2), n)], row.names = FALSE)

G <- rbind(A, B)
G[, panel := factor(panel, levels = c("Raw whiff rate (%)", "Whiff residual above model (pp)"))]
hl <- data.table(panel = factor("Whiff residual above model (pp)", levels = levels(G$panel)), y = 0)

gg <- ggplot(G, aes(x, y)) +
  geom_hline(data = hl, aes(yintercept = y), linetype = 2, colour = "grey45") +
  geom_ribbon(aes(ymin = y - 1.96*se, ymax = y + 1.96*se), fill = "#1d7870", alpha = .18) +
  geom_line(colour = "#1d7870", linewidth = .9) +
  geom_point(colour = "#1d7870", size = 2.1) +
  facet_wrap(~panel, scales = "free_y") +
  labs(title = "Velocity separation drives the whiff rate. It does not drive what is left over after the model sees it",
       subtitle = paste0("Every changeup swing, 2020-2026, in deciles of velocity separation from the pitcher's primary fastball. LEFT: the raw whiff rate climbs steadily with\n",
                         "separation, which is the relationship everyone already knows. RIGHT: the same swings against the residual from a shape-and-location model that trains on\n",
                         "velocity separation as one of its 24 features. The residual curve is flat, which is the correct result - it means the model absorbed the velocity effect\n",
                         "rather than leaving it on the table. Bands are 95 percent intervals on the decile mean."),
       x = "Velocity separation from the primary fastball (mph slower)", y = NULL,
       caption = "Source: Statcast 2020-2026 - 464,710 changeups - out-of-fold LightGBM residuals") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12.5), plot.subtitle = element_text(size = 8.2),
        strip.text = element_text(face = "bold", size = 10), panel.grid.minor = element_blank(),
        plot.caption = element_text(size = 7.5, colour = "grey40"))
ggsave(file.path(AST, "fig28_whiff_resid_vs_velo.png"), gg, width = 11, height = 5.6, dpi = 150)
cat("\nwrote fig28_whiff_resid_vs_velo.png\n")

#!/usr/bin/env Rscript

# WHERE THE SPIN CUE WENT.
#
# Every positive spin result in this project needed two conditions at once: the cue had to
# be axis geometry alone, and location had to be either absent from the whiff model or
# bolted on afterwards as a cubic surface. Break either condition and the effect is zero.
#
# The nested models from axis_after_tunnel.R are the cleanest evidence, because they change
# one thing at a time on a fixed set of pitchers:
#
#   M0  base shape                      residual then detrended for location post-hoc
#   M1  + path_ratio                    tunnel priced in, location still post-hoc
#   M2  + location & handedness         location moves INSIDE the model
#   M3  + path_ratio + location         both inside
#   M4  M3 + axis_diff itself           the cue given to the model directly
#
# M0 -> M1 adds the tunnel and nothing happens. M1 -> M2 moves location inside and the
# effect vanishes for every pitch type that had one. The lower panel shows why: the tunnel
# buys the model almost no AUC, location buys it twenty points.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(grid) })
MDIR <- "data/statcast_model"; AST <- file.path(MDIR, "article_assets")

LAB <- c(sliders="Sliders", sweepers="Sweepers", curves="Curveballs",
         changeups="Changeups", splitters="Splitters")
PAL <- c(Sliders="#8d99ae", Sweepers="#5c6b80", Curveballs="#3d4a5c",
         Changeups="#2a9d8f", Splitters="#14594f")
STAGE <- c(M0="M0\nbase shape", M1="M1\n+ tunnel", M2="M2\n+ location",
           M3="M3\n+ tunnel\n+ location", M4="M4\n+ axis_diff\nas a feature")

# The saved tables store each cell as "+0.343 (p=2.4e-08  n=251)".
parse_cell <- function(s) {
  m <- regmatches(s, regexec("([+-][0-9.]+)\\s*\\(p=([0-9.eE+-]+)\\s+n=([0-9]+)\\)", s))[[1]]
  if (length(m) != 4) return(list(r = NA_real_, p = NA_real_, n = NA_integer_))
  list(r = as.numeric(m[2]), p = as.numeric(m[3]), n = as.integer(m[4]))
}

rows <- rbindlist(lapply(c("breaking","offspeed"), function(g) {
  a <- readRDS(file.path(MDIR, sprintf("axis_after_tunnel_%s.rds", g)))
  M <- a$models[cue == "axis_sim"]
  tys <- setdiff(names(M), c("cue","model","what"))
  rbindlist(lapply(tys, function(ty) rbindlist(lapply(seq_len(nrow(M)), function(i) {
    v <- parse_cell(M[[ty]][i])
    data.table(group = g, type = LAB[[ty]], model = M$model[i],
               r = v$r, p = v$p, np = v$n)
  }))))
}))
rows[, `:=`(type = factor(type, levels = names(PAL)),
            model = factor(model, levels = names(STAGE)),
            sig = p < .05)]
fwrite(rows, file.path(AST, "ext_axis_collapse_trace.csv"))
cat("=== axis-match correlation across the nested models ===\n")
print(dcast(rows, type + np ~ model, value.var = "r"), row.names = FALSE)

fit <- rbindlist(lapply(c("breaking","offspeed"), function(g) {
  a <- readRDS(file.path(MDIR, sprintf("axis_after_tunnel_%s.rds", g)))
  data.table(group = fifelse(g == "breaking", "Breaking-ball whiff model",
                                              "Offspeed whiff model"),
             model = factor(names(a$auc), levels = names(STAGE)), auc = as.numeric(a$auc))
}))

BAND <- annotate("rect", xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = Inf,
                 fill = "#c0392b", alpha = .055)

# The M2 values cluster within 0.08 of each other, so the labels get pushed apart
# vertically and joined back to their points with a leader line.
spread <- function(y, gap) {
  o <- order(y); z <- y[o]
  for (i in 2:length(z)) if (z[i] - z[i-1] < gap) z[i] <- z[i-1] + gap
  z[order(o)] - mean(z) + mean(y)
}
m2 <- rows[model == "M2"][order(-r)][, ly := spread(r, 0.052)]
BASE <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 8.4), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), legend.position = "top",
        legend.key.width = unit(26,"pt"), axis.text.x = element_text(size = 8.2))

p1 <- ggplot(rows, aes(model, r, colour = type, group = type)) +
  BAND +
  geom_hline(yintercept = 0, colour = "black", linewidth = .4) +
  geom_line(linewidth = 1.05) +
  geom_point(aes(shape = sig, fill = type), size = 3.1, stroke = 1.05) +
  geom_text(data = rows[model == "M0"], aes(label = sprintf("%+.2f", r)),
            vjust = -1.3, size = 2.9, fontface = "bold", show.legend = FALSE) +
  geom_segment(data = m2, aes(x = 3.05, xend = 3.19, y = r, yend = ly),
               linewidth = .3, show.legend = FALSE) +
  geom_text(data = m2, aes(x = 3.21, y = ly, label = sprintf("%+.2f", r)),
            hjust = 0, size = 2.9, fontface = "bold", show.legend = FALSE) +
  annotate("segment", x = 2.5, xend = 2.5, y = -0.22, yend = 0.47,
           linetype = "dashed", colour = "#c0392b", linewidth = .5) +
  annotate("text", x = 2.56, y = 0.465, hjust = 0, size = 3, colour = "#c0392b",
           fontface = "bold", label = "location moves inside the model") +
  annotate("text", x = 2.44, y = 0.465, hjust = 1, size = 3, colour = "grey40",
           label = "location handled post-hoc") +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), guide = "none") +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_fill_manual(values = PAL, guide = "none") +
  scale_x_discrete(labels = STAGE) +
  coord_cartesian(ylim = c(-0.23, 0.52)) +
  labs(title = "The spin cue needed location to be mishandled, and the tunnel had nothing to do with it",
       subtitle = paste0("Between-pitcher Spearman correlation of axis match to the fastball (-axis gap) with whiff overperformance, one pitcher per point, min 200 pitches.\n",
                         "Each step changes exactly one thing in the whiff model that generates the residual. Filled circles are p < .05. M0 and M1 are the condition the retracted\n",
                         "Figure 11b published under: shape-only model, location removed afterwards with a cubic surface. Adding the tunnel metric at M1 moves nothing. Moving location\n",
                         "inside the model at M2 erases every significant effect - sliders +0.34 to -0.10, sweepers +0.40 to -0.02, splitters +0.30 to +0.05. Changeups are the exception\n",
                         "that proves the point: flat at +0.12 throughout, and never significant at any stage, because there was never anything there to remove."),
       x = NULL, y = "Correlation of axis match\nwith overperformance") + BASE

p2 <- ggplot(fit, aes(model, auc, colour = group, group = group)) +
  BAND +
  geom_line(linewidth = 1.05) + geom_point(size = 2.8) +
  geom_text(aes(label = sprintf("%.3f", auc),
                vjust = fifelse(group == "Breaking-ball whiff model", -1.35, 2.15)),
            size = 2.8, show.legend = FALSE) +
  annotate("segment", x = 2.5, xend = 2.5, y = .55, yend = .86,
           linetype = "dashed", colour = "#c0392b", linewidth = .5) +
  scale_colour_manual(values = c(`Breaking-ball whiff model` = "#5c6b80",
                                 `Offspeed whiff model` = "#2a9d8f"), name = NULL) +
  scale_x_discrete(labels = STAGE) + coord_cartesian(ylim = c(.55, .87)) +
  labs(subtitle = paste0("Why: out-of-fold AUC of the same five models. The tunnel metric is worth +0.003 and +0.006 AUC. Location is worth +0.204 and +0.162. The post-hoc cubic\n",
                         "surface at M0/M1 cannot remove what a 25-feature model can, and the handedness-shaped location value it leaves behind is what the axis cue was reading."),
       x = "Whiff model that generates the residual", y = "Out-of-fold AUC") + BASE

png(file.path(AST, "fig13_axis_collapse.png"), width = 13.5, height = 9.4, units = "in", res = 150)
grid.newpage(); pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(5.9, 3.5), "in"))))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

cat("\n=== the two conditions, stated as a table ===\n")
cond <- rows[model %in% c("M1","M2"), .(type, np,
             `location post-hoc (M1)` = sprintf("%+.3f%s", r, fifelse(p < .05, " *", "")),
             model)]
print(dcast(rows[model %in% c("M1","M2")], type + np ~ model,
            value.var = "r")[, .(type, pitchers = np, post_hoc = M1, in_model = M2,
                                 drop = round(M1 - M2, 3))], row.names = FALSE)
cat("\nwrote fig13_axis_collapse.png and ext_axis_collapse_trace.csv\n")

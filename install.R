#!/usr/bin/env Rscript
# Install the R packages used across baseball/*.R. `splines` and `grid` ship with base R.
pkgs <- c("data.table", "lightgbm", "ggplot2", "bit64", "mgcv", "ggrepel",
          "patchwork", "jsonlite", "curl", "scales", "dplyr")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
cat("R packages ready:", paste(pkgs, collapse = ", "), "\n")

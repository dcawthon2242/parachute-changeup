#!/usr/bin/env Rscript

# Turn the scraped NCAA rosters into D1PitcherHeights, keyed the way the TrackMan data is keyed.
#
# The awkward step is that TrackMan names teams as ARI_WIL and LOU_CAR while the rosters are keyed
# on school names. Inverting the code is unreliable - it is usually three letters of the school and
# three of the nickname, but Virginia Tech arrives as VIR_TEC, which takes both halves from the
# school. A nickname-based heuristic resolved fewer than a third of the codes confidently.
#
# Roster overlap settles it without any naming rule. A TrackMan team throws twenty-odd distinct
# pitchers in a season and the correct school's roster contains nearly all of them, while any other
# school contains almost none, so the match is decided by a wide margin and can be verified by
# looking at that margin. Names are compared as an unordered set of tokens, which keeps
# "Van Allen, Cade" and "Cade Van Allen" together.
#
# Output columns are the ones PitchReports.R expects - Pitcher, PitcherTeam, PitcherHeight,
# PitcherThrows - and Pitcher/PitcherTeam are written as TrackMan's own strings so the join
# downstream needs no normalisation.

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
options(width = 200)
DIR <- "data/ncaa_rosters"
TM  <- c(`2023` = "~/Downloads/pbp23tm.csv", `2024` = "~/Downloads/D1TM24.csv",
         `2025` = "~/Downloads/D1TM25.csv")

`%||%` <- function(a, b) if (is.null(a)) b else a

key <- function(x) {                       # order-insensitive name key
  t <- lapply(strsplit(tolower(gsub("[^a-z ,]", "", tolower(x))), "[ ,]+"), function(p) {
    p <- p[nzchar(p) & nchar(p) > 1]; paste(sort(p), collapse = "|") })
  unlist(t)
}

## ---- rosters -------------------------------------------------------------------------------
read_rosters <- function(f) {
  if (!file.exists(f)) return(NULL)
  rbindlist(lapply(readLines(f), function(l) {
    r <- fromJSON(l, simplifyDataFrame = FALSE)
    if (!is.null(r$error) || is.null(r$players) || !length(r$players)) return(NULL)
    rbindlist(lapply(r$players, function(p) data.table(
      school = r$school, year = r$year, name = p$name, pos = p$pos %||% "",
      height_in = if (is.null(p$height_in)) NA_real_ else as.numeric(p$height_in),
      throws = p$throws %||% "", bats = p$bats %||% "",
      src = r$source %||% "ncaa"))) }), fill = TRUE)
}
# Two sources, same schema. The NCAA pages carry throwing hand, which the Sidearm tables usually
# omit, so where both cover a school-season the NCAA row is kept.
R <- unique(rbind(read_rosters(file.path(DIR, "rosters.jsonl")),
                  read_rosters(file.path(DIR, "sidearm_rosters.jsonl")), fill = TRUE))
R[, k := key(name)]
setorder(R, school, year, k, -src)          # "ncaa" sorts after "sidearm", so -src puts it first
R <- unique(R, by = c("school","year","k"))
cat(sprintf("source mix: %s\n", paste(sprintf("%s %d", names(table(R$src)), table(R$src)),
                                      collapse = ", ")))
cat(sprintf("rosters: %d player-seasons, %d schools, %.1f%% with a height\n",
            nrow(R), uniqueN(R$school), 100*mean(is.finite(R$height_in))))

## ---- trackman pitchers ---------------------------------------------------------------------
P <- rbindlist(lapply(names(TM), function(y) {
  d <- fread(TM[[y]], select = c("Pitcher","PitcherId","PitcherTeam","PitcherThrows","Level"),
             nThread = 4)
  d <- d[Level == "D1" & nzchar(Pitcher)]
  d[, year := as.integer(y)]
  d[, .(pitches = .N, throws_tm = PitcherThrows[1], PitcherId = PitcherId[1]),
    by = .(year, PitcherTeam, Pitcher)] }))
P[, k := key(Pitcher)]
cat(sprintf("trackman: %d pitcher-team-seasons, %d team codes\n", nrow(P), uniqueN(P$PitcherTeam)))

## ---- resolve team codes by roster overlap ----------------------------------------------------
RK <- unique(R[, .(school, year, k)])
PK <- unique(P[pitches >= 25, .(PitcherTeam, year, k)])
J  <- merge(PK, RK, by = c("year","k"), allow.cartesian = TRUE)
HIT <- J[, .(hits = .N), by = .(PitcherTeam, year, school)]
TOT <- PK[, .(n = .N), by = .(PitcherTeam, year)]
M   <- merge(HIT, TOT, by = c("PitcherTeam","year"))[, share := hits/n]
setorder(M, PitcherTeam, year, -share)
BEST <- M[, .(school = school[1], share = share[1], n = n[1],
              margin = share[1] - fifelse(.N > 1, share[2], 0)), by = .(PitcherTeam, year)]
cat(sprintf("\nteam codes resolved: %d of %d team-seasons at >=50%% roster overlap\n",
            BEST[share >= .5, .N], nrow(BEST)))
print(BEST[, .(team_seasons = .N, median_overlap = round(median(share),3),
               median_margin = round(median(margin),3)),
           by = .(quality = fifelse(share >= .7, "strong (>=70%)",
                            fifelse(share >= .4, "moderate", "weak (<40%)")))], row.names = FALSE)

# One school per code, decided across seasons, so a thin year inherits the confident answer.
CODE <- BEST[share >= .4, .(school = school[which.max(share*n)], seasons = .N,
                            best_overlap = round(max(share),3)), by = PitcherTeam]
fwrite(BEST[order(share)], file.path(DIR, "team_map_overlap.csv"))
cat(sprintf("mapped %d of %d codes to a school\n", nrow(CODE), uniqueN(P$PitcherTeam)))
cat("\nweakest 12 resolved codes (check these):\n")
print(head(BEST[share >= .4][order(share)], 12), row.names = FALSE)
cat("\nunresolved codes:\n")
print(BEST[!PitcherTeam %in% CODE$PitcherTeam][order(-n)][1:15], row.names = FALSE)

## ---- attach heights -------------------------------------------------------------------------
P2 <- merge(P, CODE[, .(PitcherTeam, school)], by = "PitcherTeam", all.x = TRUE)
H  <- R[is.finite(height_in), .(height_in = median(height_in), throws_r = throws[1],
                                pos = pos[1]), by = .(school, year, k)]
P3 <- merge(P2, H, by = c("school","year","k"), all.x = TRUE)
# a pitcher listed in one season but not another still has a height; carry it within school
CARRY <- H[, .(h_any = median(height_in)), by = .(school, k)]
P3 <- merge(P3, CARRY, by = c("school","k"), all.x = TRUE)
P3[, PitcherHeight := fifelse(is.finite(height_in), height_in, h_any)]

cat(sprintf("\n=== coverage ===\n%.1f%% of pitcher-team-seasons matched to a height\n",
            100*mean(is.finite(P3$PitcherHeight))))
cat(sprintf("%.1f%% weighted by pitches thrown\n",
            100*sum(P3$pitches[is.finite(P3$PitcherHeight)])/sum(P3$pitches)))
print(P3[, .(pitcher_seasons = .N, with_height = sum(is.finite(PitcherHeight)),
             pct = round(100*mean(is.finite(PitcherHeight)),1)), by = year][order(year)],
      row.names = FALSE)

OUT <- unique(P3[is.finite(PitcherHeight),
                 .(Pitcher, PitcherTeam, PitcherHeight,
                   PitcherThrows = fifelse(nzchar(throws_r), throws_r, throws_tm),
                   PitcherId, school, year, pitches)])
fwrite(OUT, file.path(DIR, "D1PitcherHeights.csv"))
cat(sprintf("\nwrote %s: %d rows, %d distinct pitchers\n",
            file.path(DIR, "D1PitcherHeights.csv"), nrow(OUT), uniqueN(OUT$PitcherId)))
cat(sprintf("height distribution: median %.0f in, range %.0f-%.0f\n",
            median(OUT$PitcherHeight), min(OUT$PitcherHeight), max(OUT$PitcherHeight)))
print(head(OUT[order(-pitches)], 8), row.names = FALSE)

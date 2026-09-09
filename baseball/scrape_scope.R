#!/usr/bin/env Rscript

# Which of the remaining team-seasons actually have to be scraped?
#
# 705 of 917 are still outstanding and the host is rate limiting, so the cheapest win is to not
# fetch pages whose contents are already recoverable. Two reductions are worth testing.
#
# By season: a player's height does not change, so a pitcher who appears in the 2024 data but also
# appears on the same school's 2023 or 2025 roster needs no 2024 page. The 72 schools already
# scraped have all three seasons, which makes it possible to measure exactly how much 2024 is
# covered by its neighbours rather than guessing.
#
# By school: only schools that actually supply pitchers to the analysis matter. The parachute work
# runs on pitcher-seasons that throw a qualifying changeup and four-seamer, so a school contributing
# none of those is dead weight.

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
options(width = 200)
DIR <- "data/ncaa_rosters"
TM  <- c(`2023` = "~/Downloads/pbp23tm.csv", `2024` = "~/Downloads/D1TM24.csv",
         `2025` = "~/Downloads/D1TM25.csv")
`%||%` <- function(a, b) if (is.null(a)) b else a
key <- function(x) unlist(lapply(strsplit(tolower(gsub("[^a-z ,]", "", tolower(x))), "[ ,]+"),
                                 function(p) paste(sort(p[nzchar(p) & nchar(p) > 1]), collapse = "|")))

R <- rbindlist(lapply(readLines(file.path(DIR, "rosters.jsonl")), function(l) {
  r <- fromJSON(l, simplifyDataFrame = FALSE)
  if (!is.null(r$error) || !length(r$players %||% list())) return(NULL)
  rbindlist(lapply(r$players, function(p) data.table(
    school = r$school, year = r$year, name = p$name, pos = p$pos %||% "",
    height_in = if (is.null(p$height_in)) NA_real_ else as.numeric(p$height_in)))) }), fill = TRUE)
R[, k := key(name)]
FULL <- R[, uniqueN(year), by = school][V1 == 3, school]     # schools with all three seasons
cat(sprintf("%d schools scraped, %d of them have all three seasons\n", uniqueN(R$school), length(FULL)))

## ---- reduction 1: is the middle season redundant? --------------------------------------------
cat("\n=== if 2024 rosters were skipped, how many 2024 players are still recoverable? ===\n")
P24 <- R[school %in% FULL & year == 2024 & is.finite(height_in), .(school, k)]
NB  <- unique(R[school %in% FULL & year %in% c(2023, 2025) & is.finite(height_in), .(school, k)])
P24[, found := paste(school, k) %in% paste(NB$school, NB$k)]
cat(sprintf("  all 2024 players: %.1f%% appear on the same school's 2023 or 2025 roster\n",
            100*mean(P24$found)))
PIT <- R[school %in% FULL & year == 2024 & is.finite(height_in) &
         grepl("P", toupper(pos)), .(school, k)]
PIT[, found := paste(school, k) %in% paste(NB$school, NB$k)]
cat(sprintf("  pitchers only:    %.1f%% recoverable (n = %d)\n", 100*mean(PIT$found), nrow(PIT)))

# The number that decides it is not roster coverage but coverage of pitchers who actually threw,
# weighted by how much they threw - a missed pitcher with 4 pitches costs nothing.
TMD <- rbindlist(lapply(names(TM), function(y) {
  d <- fread(TM[[y]], select = c("Pitcher","PitcherId","PitcherTeam","Level"), nThread = 4)
  d <- d[Level == "D1" & nzchar(Pitcher)]
  d[, year := as.integer(y)]
  d[, .(pitches = .N), by = .(year, PitcherTeam, Pitcher, PitcherId)] }))
TMD[, k := key(Pitcher)]
saveRDS(TMD, file.path(DIR, "tm_pitchers.rds"))

MAP <- fread(file.path(DIR, "team_map_overlap.csv"))[share >= .4,
        .(school = school[which.max(share*n)]), by = PitcherTeam]
T24 <- merge(TMD[year == 2024], MAP, by = "PitcherTeam")[school %in% FULL]
T24[, in_nb := paste(school, k) %in% paste(NB$school, NB$k)]
cat(sprintf("  2024 pitchers who threw, weighted by pitches: %.1f%% recoverable from 2023/2025\n",
            100*sum(T24$pitches[T24$in_nb])/sum(T24$pitches)))

## ---- reduction 2: which schools supply analysis pitchers? ------------------------------------
cat("\n=== which schools does the parachute analysis actually draw from? ===\n")
PR <- as.data.table(readRDS("data/statcast_model/ncaa_spineff_pairs.rds"))
NEED <- unique(merge(PR[, .(PitcherId, year = season, clean)],
                     TMD[, .(PitcherId, year, PitcherTeam, pitches)],
                     by = c("PitcherId","year")))
cat(sprintf("  %d candidate pitcher-seasons sit on %d team codes\n",
            nrow(NEED), uniqueN(NEED$PitcherTeam)))
cat(sprintf("  %d passing the .90 screen sit on %d team codes\n",
            NEED[clean == TRUE, .N], NEED[clean == TRUE, uniqueN(PitcherTeam)]))
BY <- NEED[clean == TRUE, .(screened = .N), by = PitcherTeam][order(-screened)]
cat(sprintf("  distribution across teams: median %d, max %d, teams with just one: %d\n",
            median(BY$screened), max(BY$screened), sum(BY$screened == 1)))
ALL <- TMD[, .(pitches = sum(pitches)), by = PitcherTeam]
cat(sprintf("  team codes with any D1 pitches at all: %d\n", nrow(ALL)))
cat(sprintf("  codes carrying NO screened pitcher: %d (%.0f%% of all codes)\n",
            nrow(ALL[!PitcherTeam %in% BY$PitcherTeam]),
            100*mean(!ALL$PitcherTeam %in% BY$PitcherTeam)))

## ---- what the two reductions leave -----------------------------------------------------------
TEAMS <- fread(file.path(DIR, "d1_teams.csv"))
DONE  <- unique(R[, .(school, year)])
TEAMS[, done := paste(team_name, year) %in% paste(DONE$school, DONE$year)]
cat(sprintf("\n=== remaining work ===\nall seasons:        %d pages left\n", TEAMS[done == FALSE, .N]))
cat(sprintf("dropping 2024:      %d pages left\n", TEAMS[done == FALSE & year != 2024, .N]))
fwrite(TEAMS[done == FALSE & year != 2024][order(team_name, -year)],
       file.path(DIR, "d1_teams_priority.csv"))
cat(sprintf("wrote %s\n", file.path(DIR, "d1_teams_priority.csv")))

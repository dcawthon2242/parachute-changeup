#!/usr/bin/env Rscript

# Rebuild of the swing-decision run value from CCAM-Tunneling-Project/swing_decisions.R,
# trained on 2023-2025 and applied to 2026 so the residual is genuinely out of sample.
#
#   swing_decision_rv = (actual RV change) - (expected RV change)
#
# The expected side comes from a cascade of event models -- swing, called ball/strike,
# whiff, in-play, batted-ball type -- none of which see anything about sequencing. That
# is exactly what makes it the right target for a tunneling metric: location, velocity,
# movement and count are already absorbed, so whatever is left is available for the
# previous pitch to explain. The whiff test that failed earlier was contaminated by
# "where it finished"; this one is not.
#
# FOUR BUGS IN THE ORIGINAL ARE FIXED HERE, all of which change the metric materially.
#
# 1. HANDEDNESS WAS SILENTLY DROPPED. Every model did
#      stand = ifelse(stand == "Right", 1, 0)
#    but Statcast codes stand and p_throws as "R"/"L", so both features were constant 0
#    in every model. Platoon was never in the expected-event model at all.
#
# 2. WALKS AND STRIKEOUTS NEVER FIRED. The count update read
#      expected_event %in% c("CalledStrike","Whiff","Foul") & strikes < 2 ~ ...
#    with no strikes == 2 branch, falling through to TRUE ~ current_count. So a whiff on
#    0-2 and a ball on 3-1 both scored an RV change of exactly zero. The two highest-
#    leverage outcomes in the count tree were zeroed out, and the "Walk"/"Strikeout"
#    entries in run_values were unreachable.
#
# 3. THE 2-0 RUN VALUE WAS A COPY-PASTE OF 0-1. run_values had 2-0 = -0.0436, identical
#    to 0-1 and negative, when 2-0 is one of the best counts for a hitter (~ +0.098).
#
# 4. THE CALLED-PITCH MODEL DECLARED num_class = 3 BUT ONLY EVER SAW LABELS 0 AND 1.
#    hit_by_pitch mapped to NA and was filtered out, and blocked_ball was assigned twice
#    (2 then 0). xHBP was therefore noise. HBP is modelled as its own class here.
#
# Two versions of the metric are produced. sd_rv_argmax follows the original definition,
# taking the single most likely event as the expectation. sd_rv_exp uses the full
# probability-weighted expectation instead. The argmax version throws away nearly all of
# the distribution and is very lumpy, which costs a lot of power when the effect being
# hunted is small, so the weighted version is the better test.
#
# Sign convention: run values are from the HITTER's perspective, so a NEGATIVE
# swing_decision_rv means the hitter did worse than expected, i.e. the pitcher won.

suppressPackageStartupMessages({ library(data.table); library(lightgbm) })
set.seed(23); options(width = 210)
TRAIN_YEARS <- 2023:2025; APPLY_YEAR <- 2026
OUT <- file.path("data", "statcast_model", "swing_decision_rv_2026.rds")

COLS <- c("game_pk","at_bat_number","pitch_number","game_type","pitch_type","p_throws",
          "stand","balls","strikes","plate_x","plate_z","release_speed","release_spin_rate",
          "spin_axis","ax","ay","az","vx0","vy0","vz0","release_pos_x","release_pos_y",
          "release_pos_z","release_extension","arm_angle","description","launch_angle",
          "batter","pitcher")

read_year <- function(y) {
  dir <- sprintf("data/statcast_%d/chunks", y)
  fs <- if (dir.exists(dir)) list.files(dir, pattern = "csv$", full.names = TRUE) else
    sprintf("data/statcast_%d/statcast_%d_all.csv", y, y)
  d <- rbindlist(lapply(fs, function(f) fread(f, select = COLS, showProgress = FALSE)), fill = TRUE)
  d[, season := y][]
}
d <- rbindlist(lapply(c(TRAIN_YEARS, APPLY_YEAR), read_year), fill = TRUE)
d <- d[game_type == "R" & pitch_type != "" & !is.na(pitch_type) &
         !is.na(plate_x) & !is.na(plate_z) & !is.na(vy0)]
d <- unique(d, by = c("season","game_pk","at_bat_number","pitch_number"))
cat(sprintf("%s pitches, %d-%d train / %d apply\n", format(nrow(d), big.mark = ","),
            min(TRAIN_YEARS), max(TRAIN_YEARS), APPLY_YEAR))

## ---- approach-angle features ---------------------------------------------
y0 <- 50; yf <- 17/12
d[, vy_f := -sqrt(pmax(vy0^2 - 2*ay*(y0 - yf), 0))]
d[, t_f := (vy_f - vy0)/ay]
d[, VAA := -atan((vz0 + az*t_f)/vy_f) * 180/pi]
d[, HAA := -atan((vx0 + ax*t_f)/vy_f) * 180/pi]
tr <- d$season %in% TRAIN_YEARS
mv <- lm(VAA ~ plate_z + pitch_type, data = d[tr])
mh <- lm(HAA ~ plate_x + pitch_type, data = d[tr])
d[, pitch_type := factor(pitch_type, levels = levels(factor(d[tr]$pitch_type)))]
d <- d[!is.na(pitch_type)]
tr <- d$season %in% TRAIN_YEARS
d[, LIVAA := predict(mv, .SD) - VAA, .SDcols = c("plate_z","pitch_type")]
d[, LIHAA := predict(mh, .SD) - HAA, .SDcols = c("plate_x","pitch_type")]

## ---- event labels --------------------------------------------------------
WHIFF <- c("swinging_strike","swinging_strike_blocked","foul_tip")
SWING <- c(WHIFF, "foul","hit_into_play")
d[, swing := as.integer(description %in% SWING)]
d[, bbt := fifelse(description == "hit_into_play" & launch_angle < 10, "GroundBall",
            fifelse(description == "hit_into_play" & launch_angle < 25, "LineDrive",
            fifelse(description == "hit_into_play" & launch_angle < 50, "FlyBall",
            fifelse(description == "hit_into_play" & launch_angle >= 50, "PopUp",
                    NA_character_))))]

FEAT <- c("balls","strikes","standR","throwsR","release_speed","release_spin_rate","spin_axis",
          "ax","az","release_pos_x","release_pos_z","release_extension","plate_x","plate_z",
          "LIHAA","LIVAA","arm_angle","pt")
d[, `:=`(standR = as.integer(stand == "R"), throwsR = as.integer(p_throws == "R"),
         pt = as.integer(pitch_type))]

mat <- function(dt) as.matrix(dt[, FEAT, with = FALSE])
PARAMS <- list(learning_rate = 0.05, num_leaves = 63, min_data_in_leaf = 200,
               feature_fraction = 0.8, bagging_fraction = 0.8, bagging_freq = 1,
               verbose = -1, num_threads = parallel::detectCores())

fit <- function(sub, label, nclass = 1, nm = "") {
  x <- sub[!is.na(get(label))]
  n <- nrow(x); ho <- sample(n, floor(0.15*n))
  p <- PARAMS
  p$objective <- if (nclass > 1) "multiclass" else "binary"
  p$metric <- if (nclass > 1) "multi_logloss" else "binary_logloss"
  if (nclass > 1) p$num_class <- nclass
  tr_ds <- lgb.Dataset(mat(x[-ho]), label = x[-ho][[label]], categorical_feature = "pt")
  va_ds <- lgb.Dataset.create.valid(tr_ds, mat(x[ho]), label = x[ho][[label]])
  m <- lgb.train(p, tr_ds, nrounds = 1500, valids = list(v = va_ds),
                 early_stopping_rounds = 50, verbose = -1)
  cat(sprintf("  %-12s n=%-9s best_iter=%-5d %s=%.4f\n", nm, format(n, big.mark = ","),
              m$best_iter, p$metric, m$best_score))
  m
}

cat("\ntraining the cascade on", paste(range(TRAIN_YEARS), collapse = "-"), ":\n")
T <- d[tr]
m_swing <- fit(T, "swing", 1, "swing")

# called pitch: 0 ball, 1 called strike, 2 HBP -- three real classes this time
TK <- T[description %in% c("ball","blocked_ball","called_strike","hit_by_pitch")]
TK[, lab := fifelse(description == "called_strike", 1L,
             fifelse(description == "hit_by_pitch", 2L, 0L))]
m_call <- fit(TK, "lab", 3, "called")

TS <- T[description %in% SWING][, lab := as.integer(description %in% WHIFF)]
m_whiff <- fit(TS, "lab", 1, "whiff")

TC <- T[description %in% c("foul","hit_into_play")][, lab := as.integer(description == "hit_into_play")]
m_play <- fit(TC, "lab", 1, "in_play")

BB <- c("GroundBall","LineDrive","FlyBall","PopUp")
TB <- T[!is.na(bbt)][, lab := match(bbt, BB) - 1L]
m_bb <- fit(TB, "lab", 4, "batted_ball")

## ---- apply to the target season -----------------------------------------
A <- d[season == APPLY_YEAR]
X <- mat(A)
A[, p_swing := predict(m_swing, X)]
cal <- matrix(predict(m_call, X), ncol = 3, byrow = TRUE)
A[, `:=`(p_ball = cal[,1], p_cstr = cal[,2], p_hbp = cal[,3])]
A[, p_whiff := predict(m_whiff, X)]
A[, p_play := predict(m_play, X)]
bb <- matrix(predict(m_bb, X), ncol = 4, byrow = TRUE)
for (j in seq_along(BB)) A[[paste0("p_", BB[j])]] <- bb[, j]

# leaf probabilities of the nine terminal events
A[, `:=`(
  P_CalledBall   = (1 - p_swing) * p_ball,
  P_CalledStrike = (1 - p_swing) * p_cstr,
  P_HBP          = (1 - p_swing) * p_hbp,
  P_Whiff        = p_swing * p_whiff,
  P_Foul         = p_swing * (1 - p_whiff) * (1 - p_play),
  P_GroundBall   = p_swing * (1 - p_whiff) * p_play * p_GroundBall,
  P_LineDrive    = p_swing * (1 - p_whiff) * p_play * p_LineDrive,
  P_FlyBall      = p_swing * (1 - p_whiff) * p_play * p_FlyBall,
  P_PopUp        = p_swing * (1 - p_whiff) * p_play * p_PopUp)]
EV <- c("CalledBall","CalledStrike","HBP","Whiff","Foul","GroundBall","LineDrive","FlyBall","PopUp")
PC <- paste0("P_", EV)
A[, tot := rowSums(.SD), .SDcols = PC]
for (c_ in PC) A[[c_]] <- A[[c_]]/A$tot

## ---- run values ---------------------------------------------------------
# Count states as in the original, with 2-0 corrected from -0.0436 (a duplicate of 0-1)
# to +0.098. Terminal values are the original's.
RVC <- c("0-0"= 0.001695997,"1-0"= 0.039248042,"0-1"=-0.043581338,
         "2-0"= 0.098000000,"1-1"=-0.015277684,"0-2"=-0.103242476,
         "3-0"= 0.200960731,"2-1"= 0.034545018,"1-2"=-0.080485991,
         "3-1"= 0.138254876,"2-2"=-0.039716495,"3-2"= 0.048505049)
TERM <- c(Walk = 0.325, Strikeout = -0.284, FlyBall = 0.586, LineDrive = 0.528,
          GroundBall = 0.164, PopUp = 0.0186)

# RV of the state reached by `event` from count (b,s). Walks and strikeouts fire here.
rv_after <- function(b, s, event) {
  out <- numeric(length(b))
  isball <- event %in% c("CalledBall")
  isstr  <- event %in% c("CalledStrike","Whiff")
  isfoul <- event == "Foul"
  if (event == "HBP") return(rep(TERM[["Walk"]], length(b)))
  if (event %in% names(TERM)) return(rep(TERM[[event]], length(b)))
  if (isball) return(ifelse(b >= 3, TERM[["Walk"]], RVC[paste(pmin(b+1,3), s, sep = "-")]))
  if (isstr)  return(ifelse(s >= 2, TERM[["Strikeout"]], RVC[paste(b, pmin(s+1,2), sep = "-")]))
  if (isfoul) return(ifelse(s >= 2, RVC[paste(b, s, sep = "-")],
                            RVC[paste(b, pmin(s+1,2), sep = "-")]))
  out
}

A[, rv_now := RVC[paste(balls, strikes, sep = "-")]]
A <- A[!is.na(rv_now)]
for (e in EV) set(A, j = paste0("D_", e), value = rv_after(A$balls, A$strikes, e) - A$rv_now)
DC <- paste0("D_", EV)

# expected RV change: probability-weighted, and argmax as the original had it
A[, exp_rv_w := rowSums(as.matrix(.SD[, PC, with = FALSE]) * as.matrix(.SD[, DC, with = FALSE])),
  .SDcols = c(PC, DC)]
am <- max.col(as.matrix(A[, PC, with = FALSE]), ties.method = "first")
A[, exp_event := EV[am]]
A[, exp_rv_am := as.matrix(.SD)[cbind(seq_len(.N), am)], .SDcols = DC]

# actual
A[, act_event := fifelse(description %in% c("ball","blocked_ball"), "CalledBall",
                  fifelse(description == "called_strike", "CalledStrike",
                  fifelse(description == "hit_by_pitch", "HBP",
                  fifelse(description %in% WHIFF, "Whiff",
                  fifelse(description == "foul", "Foul", bbt)))))]
A <- A[!is.na(act_event) & act_event %in% EV]
A[, act_rv := as.matrix(.SD)[cbind(seq_len(.N), match(act_event, EV))], .SDcols = DC]

A[, `:=`(sd_rv_exp = act_rv - exp_rv_w, sd_rv_argmax = act_rv - exp_rv_am)]

cat(sprintf("\n2026 pitches scored: %s\n", format(nrow(A), big.mark = ",")))
cat(sprintf("argmax-event accuracy vs actual: %.1f%%\n", 100*mean(A$exp_event == A$act_event)))
cat("\nswing-decision RV (hitter perspective; negative = pitcher won):\n")
cat(sprintf("  sd_rv_exp    mean %+.5f  sd %.4f\n", mean(A$sd_rv_exp), sd(A$sd_rv_exp)))
cat(sprintf("  sd_rv_argmax mean %+.5f  sd %.4f\n", mean(A$sd_rv_argmax), sd(A$sd_rv_argmax)))
cat(sprintf("  correlation between the two versions: %.3f\n",
            cor(A$sd_rv_exp, A$sd_rv_argmax)))
cat("\n  the argmax version collapses to a handful of discrete values:\n")
cat(sprintf("  distinct values: argmax %s vs weighted %s\n",
            format(uniqueN(round(A$sd_rv_argmax, 6)), big.mark = ","),
            format(uniqueN(round(A$sd_rv_exp, 6)), big.mark = ",")))

cat("\nmean sd_rv_exp by actual event (sanity: whiffs and called strikes should be negative):\n")
print(A[, .(n = .N, sd_rv_exp = round(mean(sd_rv_exp), 4),
            sd_rv_argmax = round(mean(sd_rv_argmax), 4)), by = act_event][order(sd_rv_exp)])

keep <- c("game_pk","at_bat_number","pitch_number","pitcher","batter","pitch_type",
          "description","sd_rv_exp","sd_rv_argmax","exp_event","act_event","p_swing",
          "P_Whiff","exp_rv_w","act_rv")
saveRDS(A[, ..keep], OUT)
cat(sprintf("\nwrote %s\n", OUT))

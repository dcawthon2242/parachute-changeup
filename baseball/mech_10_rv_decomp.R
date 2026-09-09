#!/usr/bin/env Rscript

# IS THE WHIFF ADVANTAGE USEFUL? DECOMPOSE RUN VALUE.
#
# Every locked-cell estimate in this project has been a whiff residual. The pitch is "useful"
# only if that extra miss shows up in run value, or if the run-value lag is a small-sample
# contact artifact. This splits MLB changeup RV/100 into four mutually exclusive channels:
#
#   take     no swing (balls and called strikes together - description is not on this file)
#   whiff    swinging strike
#   foul     swing, no whiff, not in play
#   bip      ball in play
#
# For the in-play channel, actual RV is replaced with the RV implied by Statcast's
# estimated_woba_using_speedangle, fit on the same changeups. If actual BIP RV is worse than
# expected, the deficit is real contact quality. If expected BIP RV is fine and actual is not,
# the lag is BABIP noise.

suppressPackageStartupMessages({ library(data.table); library(bit64) })
options(width = 200)
MDIR <- "data/statcast_model"
L <- readRDS(file.path(MDIR, "locked_spec.rds"))$data; setDT(L)
L <- L[league == "MLB"]
F <- readRDS(file.path(MDIR, "parachute_ff.rds")); setDT(F)
F <- F[pitch_type == "CH"]
F[, id := as.character(pitcher)]
F <- merge(F, L[, .(id, season, bin, axis, arm, ec, ef, ff_use, si_use)],
           by = c("id","season"))
F[, ch := fifelse(!is_swing, "take",
           fifelse(whiff == 1, "whiff",
           fifelse(is_bip == TRUE, "bip", "foul")))]
F[, ch := factor(ch, levels = c("take","whiff","foul","bip"))]

# Expected BIP RV: linear map from xwOBA, fit only on changeups in play with a value.
BIP <- F[ch == "bip" & is.finite(estimated_woba_using_speedangle) & is.finite(rv)]
fit <- lm(rv ~ estimated_woba_using_speedangle, BIP)
F[, xrv := rv]
F[ch == "bip" & is.finite(estimated_woba_using_speedangle),
  xrv := as.numeric(predict(fit, .SD))]
cat(sprintf("xwOBA -> RV on changeups in play: intercept %.3f, slope %.3f, R2 %.3f, n=%s\n",
            coef(fit)[1], coef(fit)[2], summary(fit)$r.squared, format(nrow(BIP), big.mark = ",")))

decomp <- function(D, lab) {
  n <- nrow(D)
  chn <- D[, .(n = .N, share = .N / n, rv100 = 100 * mean(rv, na.rm = TRUE),
               xrv100 = 100 * mean(xrv, na.rm = TRUE)), by = ch]
  tot <- data.table(ch = factor("ALL", levels = c(levels(D$ch), "ALL")),
                    n = n, share = 1,
                    rv100 = 100 * mean(D$rv, na.rm = TRUE),
                    xrv100 = 100 * mean(D$xrv, na.rm = TRUE))
  # Contribution of each channel to the total: 100 * sum(rv_channel) / N.
  chn[, contrib := share * rv100]
  tot[, contrib := rv100]
  out <- rbind(chn, tot, fill = TRUE)
  out[, group := lab]
  out
}

A <- decomp(F[bin == TRUE], "bin")
B <- decomp(F[bin == FALSE], "out")
C <- merge(A, B, by = "ch", suffixes = c("_bin","_out"))
C[, `:=`(d_rv = rv100_bin - rv100_out, d_xrv = xrv100_bin - xrv100_out,
         d_share = 100 * (share_bin - share_out),
         d_contrib = contrib_bin - contrib_out)]

cat("\n=== RV/100 by channel, locked four-seam-primary bin vs the rest of the pool ===\n")
print(C[order(ch), .(ch, n_bin, share_bin = round(100*share_bin,1),
                     rv100_bin = round(rv100_bin,2), rv100_out = round(rv100_out,2),
                     d_rv = round(d_rv,2), d_contrib = round(d_contrib,2),
                     xrv100_bin = round(xrv100_bin,2), d_xrv = round(d_xrv,2))],
      row.names = FALSE)

cat("\n  d_contrib is the piece of the total RV/100 gap that lives in that channel.\n")
cat(sprintf("  actual total gap:  %+.2f RV/100\n", C[ch == "ALL", d_rv]))
cat(sprintf("  expected total gap (BIP replaced by xwOBA-implied RV): %+.2f RV/100\n",
            C[ch == "ALL", d_xrv]))

# Contact quality on BIP only, so the rate mix does not leak in.
cat("\n=== balls in play only ===\n")
BA <- F[ch == "bip" & bin == TRUE]
BO <- F[ch == "bip" & bin == FALSE]
cat(sprintf("  bin  n=%s  actual RV/100 %+6.2f  xwOBA %.3f  xRV/100 %+6.2f\n",
            format(nrow(BA), big.mark = ","), 100*mean(BA$rv, na.rm = TRUE),
            mean(BA$estimated_woba_using_speedangle, na.rm = TRUE),
            100*mean(BA$xrv, na.rm = TRUE)))
cat(sprintf("  out  n=%s  actual RV/100 %+6.2f  xwOBA %.3f  xRV/100 %+6.2f\n",
            format(nrow(BO), big.mark = ","), 100*mean(BO$rv, na.rm = TRUE),
            mean(BO$estimated_woba_using_speedangle, na.rm = TRUE),
            100*mean(BO$xrv, na.rm = TRUE)))
t_act <- t.test(BA$rv, BO$rv)
t_x   <- t.test(BA$xrv, BO$xrv)
cat(sprintf("  actual BIP gap p = %.3f | expected BIP gap p = %.3f\n",
            t_act$p.value, t_x$p.value))

# Pitcher-season RV so one arm cannot dominate the headline.
S <- F[, .(n = .N, rv100 = 100*mean(rv), xrv100 = 100*mean(xrv)),
       by = .(id, season, bin, player_name)]
cat("\n=== pitcher-season RV/100, bin vs out ===\n")
cat(sprintf("  bin  %d seasons  actual %+6.2f +/- %.2f  expected %+6.2f\n",
            sum(S$bin), mean(S[bin == TRUE]$rv100),
            sd(S[bin == TRUE]$rv100)/sqrt(sum(S$bin)),
            mean(S[bin == TRUE]$xrv100)))
cat(sprintf("  out  %d seasons  actual %+6.2f +/- %.2f  expected %+6.2f\n",
            sum(!S$bin), mean(S[bin == FALSE]$rv100),
            sd(S[bin == FALSE]$rv100)/sqrt(sum(!S$bin)),
            mean(S[bin == FALSE]$xrv100)))
tt <- t.test(S[bin == TRUE]$rv100, S[bin == FALSE]$rv100)
tx <- t.test(S[bin == TRUE]$xrv100, S[bin == FALSE]$xrv100)
cat(sprintf("  actual gap %+6.2f p = %.3f | expected gap %+6.2f p = %.3f\n",
            tt$estimate[1] - tt$estimate[2], tt$p.value,
            tx$estimate[1] - tx$estimate[2], tx$p.value))

# Roster for the prevalence half of the claim.
cat("\n=== locked MLB bin roster, RV and expected RV ===\n")
print(S[bin == TRUE][order(-xrv100), .(player_name, season, n, rv100 = round(rv100,2),
      xrv100 = round(xrv100,2))], row.names = FALSE)

saveRDS(list(channels = C, seasons = S, xwoba_map = coef(fit)),
        file.path(MDIR, "rv_decomp.rds"))
fwrite(C, file.path(MDIR, "article_assets/ext_rv_decomp.csv"))
cat("\nwrote rv_decomp.rds\n")

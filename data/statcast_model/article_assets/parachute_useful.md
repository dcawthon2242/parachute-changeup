# Are parachute changeups useful?

Short answer: they are rare, they beat a whiff model, and they do not beat run value. The
performance claim still has to come from the pre-registered D1 2022/2026 run, not from another
slice of the data already on disk.

## Prevalence (descriptive, already solid)

Four-seam-primary population, locked gates (active spin ≥ .85 both ways, axis ≤ 10, league
top-third slot, 40+ changeup swings):

| league | bin seasons | arms | share of pool |
|---|---|---|---|
| MLB 2020-2026 | 20 | 14 | 1.8% (of 1,132) |
| D1 2023-2025 | 14 | 13 | 1.7% (of 820) |

MLB arms: Kikuchi (3 seasons), Cease, Ray, Cabrera, Skubal (2), then Sulser, Underwood Jr., May,
Smith, Jax, Detmers, Naughton, Jones, Shuster (1). This is the recipe count. It is not a
performance result.

## Usefulness (RV decomposition, new)

4,262 locked-bin changeups against 276,701 other four-seam-primary changeups.

| channel | bin share | out share | RV/100 bin | RV/100 out | contribution to gap |
|---|---|---|---|---|---|
| take (no swing) | 51.8% | 48.6% | −3.45 | −3.33 | −0.17 |
| whiff | 17.0% | 16.3% | +11.02 | +11.25 | **+0.04** |
| foul | 13.9% | 15.9% | +3.45 | +3.41 | −0.06 |
| ball in play | 17.3% | 19.3% | −4.93 | −3.30 | **−0.22** |
| **total** | | | **−0.29** | **+0.12** | **−0.41** |

Expected total, BIP replaced by xwOBA-implied RV: **−0.23**. About half the in-play hole is
expected contact quality (xwOBA 0.348 vs 0.342), half is actual-vs-expected noise. Neither BIP
gap is significant. Season-level RV gap −0.54 (p = 0.30); expected −0.29 (p = 0.49).

What this rules out: the story that they miss bats and the RV lag is BABIP luck. Expected RV is
still negative. What it does not rule out: that the true RV effect is near zero and we cannot see
it in 20 seasons. "Useful" is not shown either way.

The whiff-residual result and the RV result can both be true. The model already prices location
and shape, so a +3 MLB residual is "they miss more than this stuff should." Raw whiff rate is only
+0.7 points, and run value never sees it.

## The D1 axis transfer (also new, also thin)

At inferred efficiency ≥ .95 on both pitches, TrackMan's movement axis and Hawk-Eye's measured
axis converge, so a D1 cut can mean the same thing as MLB's. MLB's own .95 pool has 4.5% of
seasons at axis ≤ 10°, so D1's matching cut is 7.1°, not 10°.

| cell | n | arms | whiff residual | p |
|---|---|---|---|---|
| D1 .95 + matched axis | 5 | 5 | +5.58 ± 3.77 | 0.14 |
| D1 .95 + matched axis + top-third slot | 3 | 3 | +3.26 ± 6.33 | 0.61 |
| MLB .95 + axis ≤ 10 (same floor, not a new search) | 12 | 11 | +0.55 ± 1.31 | 0.68 |

Underpowered. Directionally fine, not a proof.

## What would actually prove it

Already written in [parachute_precommit.md](parachute_precommit.md). Add D1 2022 and 2026 under
the frozen spec. Confirm if pooled p < 0.05 and the point estimate stays above +1.2; kill below
+0.8; otherwise wait for MLB 2027. That is the performance test. Prevalence does not need it.

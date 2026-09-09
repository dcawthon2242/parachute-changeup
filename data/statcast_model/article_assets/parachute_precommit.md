# Pre-commit: confirmatory rule for the matched-axis changeup

Written before any new estimate is produced. The bin gates are unchanged. The only new rule is a
population filter. D1 2022 and 2026 are not in this run; the decision rule for those seasons is
locked here so it cannot be rewritten after they are seen.

## Population

Keep a pitcher-season only if four-seam usage is at least sinker usage (`ff_use >= si_use`).
Missing usage is dropped, not imputed. The filter applies to the pool and to the comparison group,
not just to the bin.

This is a perception correction, not a new threshold. The axis gap is measured against the
four-seamer. For a sinker-primary pitcher that is the wrong reference: the hitter's expectation is
set by the sinker.

## Bin gates (unchanged)

- changeup active spin >= .85
- four-seam active spin >= .85
- spin axis gap <= 10 degrees
- arm slot at or above the league's own 67th percentile, recomputed after the population filter
- at least 40 changeup swings

## Primary estimator

Inverse-variance pool of the two league means. The number that decides is the arm-clustered
standard error, not the season-level one. College careers overlap; treating repeated seasons as
independent overstates precision. Season-level p is reported as secondary.

## Decision rule for D1 2022 and 2026 (locked now)

When those seasons are added, using this same specification and no further tuning:

- **confirm** if the pooled p is below 0.05 and the point estimate stays above +1.2
- **kill** if the point estimate falls below +0.8
- **unresolved** otherwise; wait for MLB 2027

## Honesty note

The sinker-primary exclusion was informed by the first look at MLB, where six of 27 bin seasons
were sinker-primary and dropping them moved the 0-10 cell from +2.05 to +3.00. It is locked before
new D1 seasons. It is not a clean out-of-sample rule for the MLB rerun in this document.

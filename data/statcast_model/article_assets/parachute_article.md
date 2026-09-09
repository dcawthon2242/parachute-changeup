# The changeup every model wants to fix

There is a changeup shape that pitch design treats as a defect. It spins on almost the same axis as
the pitcher's four-seamer, so on a movement plot the two pitches point the same direction. It is
slow — ten or eleven miles an hour off the fastball — and it is thrown from a high slot with less
spin than a normal changeup, so it falls. Put an overlay in front of a pitching coordinator and the
verdict comes back quickly: this is not a changeup, it is a slower fastball. It has no shape of its
own. Kill the spin and get some real arm-side run, or scrap it and throw a splitter.

Stuff models agree. Given the pitch's velocity, spin, movement, extension and release point, a
LightGBM model trained on 203,132 MLB changeup swings predicts a 30.3% whiff rate for the arms whose
changeups match their fastball's axis most tightly. That is below the league changeup average. By
every measurement the industry uses to grade a pitch in isolation, this shape is a problem to be
solved.

It runs a 41.3% whiff rate.

That eleven-point gap is the largest miss in the entire study, and it points the wrong way from
everything the shape is supposed to indicate. This article is about why that happens, why no amount
of feature engineering fixes it, and — the part that matters most for anyone actually holding one of
these pitches — the specific condition that separates the versions that turn deception into runs
from the versions that do not.

All figures below are MLB Statcast 2020–2026, four-seam-primary arms, 746 pitcher-seasons across 333
arms at a 75-swing floor unless noted.

## The shape in question

Call it the parachute. Four measurements define it against the league:

| | Velocity | Separation off FF | Spin-axis gap | Spin rate | Arm slot | Arm-angle gap |
|---|---|---|---|---|---|---|
| All changeups | 85.3 | 8.2 | 23.5° | 1797 | 38.5° | 4.1° |
| Tight axis match | 83.3 | 11.5 | **7.2°** | **1637** | **52.0°** | 4.9° |
| Seam-shifted match | 83.5 | 10.2 | 11.5° | 1652 | 41.1° | 3.6° |
| "Power changeup" | 85.3 | 8.1 | **28.6°** | 1783 | 34.8° | **7.9°** |

The bottom row is the pitch that pitch design likes. It is hard, it has a big axis gap from the
fastball, and it therefore shows genuine movement differentiation — a changeup with its own shape.
The middle rows are the pitches design wants to fix.

The power changeup is predicted for a 30.1% whiff rate. It delivers 29.2%. It is the only archetype
in the study that underperforms its model, and it does so significantly (whiff residual −1.96,
p = .0075).

## The inversion

Here is every archetype, scored on a strict temporal holdout: the model is trained on 2020–2023 and
asked about 2024–2026, which is the situation a club is actually in.

| Archetype | Seasons | Predicted whiff | Actual whiff | Miss |
|---|---|---|---|---|
| Extreme seam shift | 10 | 30.1% | 43.1% | **+13.0** |
| Tight axis match | 17 | 30.3% | 41.3% | **+11.0** |
| Seam-shifted match | 44 | 29.2% | 37.3% | +8.1 |
| Axis + separation | 56 | 29.9% | 36.2% | +6.3 |
| Traditional, low seam dev | 26 | 30.8% | 34.1% | +3.3 |
| Power changeup | 66 | 30.1% | 29.2% | −0.9 |

The ordering is almost perfectly inverted against the industry's own preference ranking. Sort the
league by how much a pitch designer would want to leave the changeup alone and you have sorted it by
how badly the model underestimates it.

Two things rule out the easy dismissals. First, this is not one lucky season for a handful of arms.
The whiff residual is a repeatable property of the pitcher: year over year it correlates at r = .577
for the axis-plus-separation group (n = 40, p = .0001) and r = .453 for the seam-shifted group
(n = 31, p = .0105). Arms who beat their model do it again the next year. Second, it is not a
sample-size artifact of the tiny cells — the broad definition holds 56 seasons and 31 arms and still
posts a +2.93 whiff residual at p = .0003.

## Why the model cannot see it

The obvious response is that the model is missing a feature. So we gave it the features.

Starting from a conventional stuff model, we added the per-pitch spin-axis gap against the fastball,
seam deviation, arm slot and arm-angle gap. Then season-mean arsenal summaries. Then interaction
products. Then, finally, an explicit archetype indicator with its thresholds frozen on the training
years. Here is how much of each group's gap those additions actually absorb:

| Archetype | Gap, stuff model | Gap, everything | Closed |
|---|---|---|---|
| Traditional, low seam dev | 3.76 | 0.73 | 81% |
| Axis + separation | 6.98 | 3.51 | 50% |
| Tight axis match | 10.73 | 6.38 | 41% |
| Seam-shifted match | 9.93 | 6.81 | 31% |
| Extreme seam shift | 16.28 | 12.41 | 24% |

The pattern is the tell. The archetype that is easiest to explain away is the traditional one, whose
edge was small to begin with. The stronger the effect, the less of it any feature set recovers. At
the extreme end, three quarters of a sixteen-point gap survives everything we can hand the model.

That is not a missing column. It is a missing frame. A stuff model grades a pitch as an object: this
much velocity, this much movement, released from here. But deception is not a property of an object,
it is a property of a *pair*. A changeup that spins on the fastball's axis is not deceptive because
of anything intrinsic to it — it is deceptive because for the first fifteen feet it is
indistinguishable from the pitch the hitter has already decided he is seeing. The only cue available
is velocity, and velocity is precisely the cue a hitter cannot resolve until his swing decision is
made.

The power changeup fails for the mirror-image reason. Its 28.6° axis gap and 7.9° arm-angle gap
mean it announces itself. The hitter gets shape information and release information early, and the
model's grade — which sees a hard pitch with real movement separation — is measuring exactly the
qualities that give the pitch away.

## Where the naive version breaks

If the story stopped there it would be a clean argument for the shape and a wrong one. The tightest
axis matches — the purest expression of everything above — miss more bats than anything else in
baseball and *still lose runs*.

Every changeup lands in one of five channels whose run value contributions sum to the season total.
Splitting the archetypes by channel shows exactly where the money goes and where it leaks. Values
are runs per 100 pitches against the rest of the population, arm-clustered:

| Archetype | Ball | Called strike | Whiff | Foul | Ball in play | Total |
|---|---|---|---|---|---|---|
| Tight axis match | +0.031 | +0.200 | +0.366 | −0.130 | **−0.898** | **−0.431** |
| Axis + separation | +0.070 | +0.022 | +0.411 | −0.045 | −0.195 | +0.263 |
| Seam-shifted match | +0.019 | −0.040 | +0.471 | −0.068 | −0.096 | +0.286 |
| Extreme seam shift | −0.120 | −0.090 | +0.712 | −0.108 | +0.385 | +0.779 |
| Power changeup | −0.016 | +0.064 | −0.392 | −0.001 | +0.226 | −0.119 |

The tight-axis group earns its whiff runs and then hands back more than twice as much on contact. At
the observed exchange rate — a whiff is worth +0.1117 runs, a ball in play −0.0380 — its
ten-point whiff edge implies +0.613 runs per 100. It delivers −0.327. It converts at negative
fifty-three percent.

The contact numbers say why. Weighted by balls in play, that group's ground-ball rate runs 8.3
percentage points below expectation (p < .01) and it concedes 7.2 runs per 100 balls in play more
than its model predicts (p < .01). A changeup that matches the fastball's axis and has no
seam-shifted movement is, physically, a
backspinning pitch thrown slowly from a high slot. When a hitter does square it, it is a fly ball.
Deception without depth gets you swings and misses and home runs in the same season.

After adjusting for count, batter handedness and the model's own prediction, that group's run value
lands at −0.645 per 100, significant at p < .05. The purist version of the parachute is a real
phenomenon and a bad pitch.

## The condition that makes it work

Hold the axis gate and the velocity separation constant, and split the group by one thing: whether
the changeup's observed movement departs from what its spin alone would produce. That departure is
seam-shifted wake, and it is the whole ballgame.

| | Seasons | Whiff residual | p | Run value / 100 | Converts at |
|---|---|---|---|---|---|
| Axis match + seam shift | 44 | **+4.00** | <.0001 | **+0.379** | 86% |
| Axis match, no seam shift | 26 | +0.88 | .4354 | +0.120 | 89% |
| Extreme seam shift | 10 | **+7.80** | <.0001 | **+0.818** | 114% |

Same axis match. Same velocity separation. The half with seam-shifted movement posts a four-point
whiff residual; the half without it posts nothing that clears significance. And the seam-shifted
half keeps its runs — its ball-in-play channel is −0.096 against the tight-axis group's −0.898, and
its contact quality is statistically indistinguishable from expectation.

The extreme cell is the proof of concept. Ten seasons, seven arms, a twelve-point whiff residual,
and xwOBA on contact 33.6 points *better* than its model predicts (p < .01). It is the only group in
the study that is simultaneously elite at missing bats and elite at surviving contact.

So the mechanism is two-part, and the two parts do different jobs. The matched axis buys the swing —
it is what makes the pitch indistinguishable out of the hand. The seam shift buys the miss and
protects the contact — it is what makes the pitch arrive somewhere other than where the spin
promised. Either alone is insufficient. Axis match without seam shift misses bats and gets hit.
Seam shift without the axis match is just a changeup.

After the full stack of adjustments — count, batter handedness, and the stuff model's own predicted
run value — the two usable definitions land at **+0.295 and +0.321 runs per 100, both p < .05**.
Note the direction of that last adjustment: subtracting the model's prediction *raises* the estimate,
because the model was betting against these pitches. Per season, at a typical 330–345 changeups,
that is +1.18 and +1.26 runs for the arm.

## It is also being thrown wrong

Separately from whether anyone should learn this pitch, the arms who already have it are giving runs
back in deployment.

The right way to ask whether a pitch is misused is not to look at its run value in a count, but at
the *gap* between it and the same pitcher's hard pitches in that count, since a changeup thrown is a
fastball not thrown. With two strikes, the league's changeup beats its own hard alternative by 0.82
runs per 100 (p < .0001). For parachute arms specifically, the gap is **1.65** (p = .004) — twice
the league edge. Yet they throw it 31.3% of the time with two strikes against 30.5% with one, when
the one-strike gap is −0.11. The pitch's advantage is concentrated almost entirely in a count where
its usage barely moves.

Location is the opposite error. Sorted by run value, the best place for this pitch is not the chase
zone:

| Location | Share | Swing % | Whiff % | RV / 100 |
|---|---|---|---|---|
| Just below the zone (0–5 in) | 12.4% | 64.1% | 42.9% | **+2.71** |
| In zone | 36.2% | 77.4% | 21.8% | +2.08 |
| Chase, below (5–12 in) | 11.9% | 38.0% | 63.9% | −0.41 |
| Buried (12+ in) | 5.6% | 8.6% | 91.9% | **−3.97** |

The buried changeup has a 91.9% whiff rate and is the second-worst outcome on the board, because
almost nobody swings at it. The band just below the zone gets a 64% swing rate at a 43% whiff rate,
and that combination is where the runs are. The instinct to bury a good changeup is costing these
arms real value.

Reallocating counts alone — no change to the pitch, no change to the location — is worth about
**+1.15 runs per season** for a parachute arm, which is roughly the same size as the pitch's entire
measured edge. There is as much value available in throwing it better as there was in having it.

## Who this is

The broad definition holds 31 arms. Tarik Skubal's 2021 and 2022 changeups are in it. So are Cole
Ragans 2024 and 2025, Dylan Cease 2021 and 2026, Jeffrey Springs across four seasons, Jesús Luzardo,
Matthew Boyd, José Suárez, Ryan Feltner, Carlos Rodón 2024, Robbie Ray 2025, Michael Lorenzen 2025.
These are not marginal arms who found a trick. Several of them are among the best pitchers in
baseball, and the pitch a model would have told them to fix is a meaningful part of why.

One important limit on where the value comes from: it is a platoon pitch. Against opposite-handed
hitters the seam-shifted group posts a +6.64 whiff residual (p < .01); against same-handed hitters,
+2.16 and not significant. Sample sizes on the same-hand side are thin, but the direction is
consistent across every archetype, and any expectation built on these numbers should be built for
the platoon matchup.

## What this does not prove

The honest boundaries, stated plainly.

The run value results are real but modest, and season-level run value is a noisy measurement — it
repeats year over year at only r = .20, against r = .67 for whiff. The whiff findings are on much
firmer ground than the run findings, and the sensible reading is that the whiff effect is
established and the run effect is the best available estimate of what it converts to.

The extreme cell that shows the strongest version of every effect is ten seasons and seven arms.
Treat it as a demonstration of the mechanism, not as a reliable effect size.

Nothing here establishes that teaching the shape to a pitcher who lacks it would produce these
results. That is a causal claim and this is observational data. Within MLB it cannot be tested at
all: only 7% of consecutive big-league seasons contain a changeup separation change large enough to
observe, and the design would need an effect three to five times larger than the one in question
before it could detect anything.

And the orthodoxy is not simply wrong. The traditional half of the matched-axis group — same axis
gap, same separation, no seam shift — shows no significant whiff edge and no run value. Pitch
design's instinct that a changeup needs something beyond a velocity difference is correct. What it
has wrong is the assumption that the something has to be a visible axis separation from the
fastball. It can instead be an invisible one: movement the seams produce that the spin does not
predict, on a pitch that looks exactly like a fastball until it doesn't.

## The prescription

If a pitcher throws a changeup that spins on his fastball's axis, the model grades it poorly, and a
coordinator has flagged it for reshaping, the sequence is:

1. **Do not fix the axis match.** It is the source of the deception, and every version of the pitch
   that has it beats its model on whiff.
2. **Check the seam deviation.** If observed movement departs from spin-predicted movement by 7.8°
   or more, the pitch is already in the group that produces runs. If it does not, that — not the
   axis — is what needs work.
3. **Keep the separation above roughly 8.7 mph and the axis gap under about 16°.** Those two
   measurements alone define the 56-season group that survives every adjustment.
4. **Throw it with two strikes.** The edge is twice the league's there and the usage does not
   reflect it.
5. **Stop burying it.** Target the band from the bottom of the zone to five inches below. The
   swing-and-miss combination there is worth more than a 92% whiff rate nobody offers at.

The pitch that looks like it needs fixing is, in its correct form, one of the few in baseball whose
value the industry's own grading tools are structurally unable to see. That is not a reason to trust
it blindly. It is a reason to stop letting a model that grades pitches one at a time make decisions
about a pitch whose entire function is relational.

#!/usr/bin/env python3
"""
Height-adjusted hitter eye points in the tunneling metric's coordinate frame
(Statcast feet; origin at the back tip of home plate, +x to the catcher's right,
+y toward the pitcher, +z up).

Builds on obm_eye_point.py, which recovered the OBM lab->Statcast registration.
This script adds the height adjustment and emits a per-batter table for 2026.

WHICH AXES GET ADJUSTED FOR HEIGHT, AND WHY
-------------------------------------------
Regressing the OBM per-athlete stance eye position on body height:

  eye height   r = +0.76      strongly height-dependent  -> ADJUST
  |lateral|    r = -0.05 (R), -0.11 (L)   no relationship -> DO NOT adjust
  depth        r = -0.27 (R), +0.14 (L)   opposite signs  -> DO NOT adjust

Eye height is modelled as a RATIO of stature rather than a linear fit. Over the
observed 63-77 in range a linear fit returns a slope near 1.0 ft per ft with a
large negative intercept, which is a narrow-range artifact that misbehaves when
extrapolated to MLB's 66-79 in. The ratio eye_z / stature is tight (SD ~3.3%)
and extrapolates safely.

LATERAL: TWO VARIANTS
---------------------
'obm'  : the OBM eye lateral directly. Clean and definitionally known, but from
         collegiate hitters.
'mlb'  : anchors on the MLB bat-tracking batter position (which tracks the feet,
         not the head) and applies the OBM head-minus-feet offset. Uses the right
         population but inherits whatever Savant means by "batter position".
The two differ by 1-2 in, inside the error bar. 'obm' is the default.
"""

import numpy as np
import pandas as pd

M_PER_FT = 0.3048
Y0 = -0.221   # lab y of the plate centreline, from the mirror-symmetry constraint
X0 = -0.736   # lab x of the back tip of home plate, from the box-centre assumption

# MLB bat-tracking implied batter (feet) lateral distance from the centreline, inches
MLB_FEET_LATERAL_IN = {"R": 36.2, "L": 35.2}

OUT_PER_BATTER = "data/obm/hitter_eye_points_2026_by_batter.csv"
OUT_SUMMARY = "data/obm/hitter_eye_points_summary.csv"


def load_obm():
    d = pd.read_csv("data/obm/stance_raw_lab_metres.csv")
    per = d.groupby(["user", "side"], as_index=False)[
        ["eye_x", "eye_y", "eye_z", "feet_x", "feet_y", "feet_z", "height_in"]].mean()
    # lab metres -> Statcast feet
    per["eye_lat"] = -(per.eye_y - Y0) / M_PER_FT
    per["eye_dep"] = (per.eye_x - X0) / M_PER_FT
    per["eye_ht"] = per.eye_z / M_PER_FT
    per["feet_lat"] = -(per.feet_y - Y0) / M_PER_FT
    per["eye_ht_ratio"] = per.eye_ht * 12 / per.height_in
    return per


def main():
    per = load_obm()

    print("=== OBM eye-height ratio (eye height / stature) ===")
    g = per.groupby("side")["eye_ht_ratio"].agg(["count", "mean", "std"])
    print(g.round(4))
    ratio = {s: g.loc[s, "mean"] for s in ("R", "L")}
    ratio_sd = {s: g.loc[s, "std"] for s in ("R", "L")}

    print("\n=== OBM height-independent axes, by side (Statcast feet) ===")
    ax = per.groupby("side")[["eye_lat", "eye_dep", "feet_lat"]].agg(["mean", "std"])
    print(ax.round(3))

    lat_obm = {s: per[per.side == s].eye_lat.mean() for s in ("R", "L")}
    dep = {s: per[per.side == s].eye_dep.mean() for s in ("R", "L")}
    # head-minus-feet lateral offset, signed toward the plate
    head_off = {s: abs(per[per.side == s].eye_lat.mean()) - abs(per[per.side == s].feet_lat.mean())
                for s in ("R", "L")}
    sgn = {"R": -1.0, "L": +1.0}   # RHH box is at negative Statcast x
    lat_mlb = {s: sgn[s] * (MLB_FEET_LATERAL_IN[s] / 12 + head_off[s]) for s in ("R", "L")}

    print("\n=== Lateral variants (Statcast feet) ===")
    for s in ("R", "L"):
        print(f"  {s}: obm {lat_obm[s]:+.3f} | head-minus-feet offset {head_off[s]:+.3f} | "
              f"mlb-anchored {lat_mlb[s]:+.3f}")

    # ---- per-batter -------------------------------------------------------
    bat = pd.read_csv("data/obm/mlb_batters_2026.csv")
    ht = pd.read_csv("data/obm/mlb_batter_heights.csv")[["batter", "name", "ht_in"]]
    b = bat.merge(ht, on="batter", how="left")
    missing = b.ht_in.isna().sum()
    if missing:
        print(f"\nWARNING: {missing} batter-side rows lack a height; filled with the side median")
        b["ht_in"] = b.groupby("stand")["ht_in"].transform(lambda s: s.fillna(s.median()))

    b["eye_x_obm"] = b.stand.map(lat_obm)
    b["eye_x_mlb"] = b.stand.map(lat_mlb)
    b["eye_y"] = b.stand.map(dep)
    b["eye_z"] = b.ht_in * b.stand.map(ratio) / 12.0
    b["eye_z_sd"] = b.ht_in * b.stand.map(ratio_sd) / 12.0
    b = b.sort_values("pitches", ascending=False)
    b.to_csv(OUT_PER_BATTER, index=False)

    print(f"\nwrote {OUT_PER_BATTER}  ({len(b)} batter-side rows)")
    print("\n=== Tallest and shortest, to show the spread that height buys ===")
    show = ["name", "stand", "ht_in", "eye_x_obm", "eye_y", "eye_z", "pitches"]
    print(b.nlargest(4, "ht_in")[show].round(3).to_string(index=False))
    print(b.nsmallest(4, "ht_in")[show].round(3).to_string(index=False))

    # ---- pitch-weighted league averages ----------------------------------
    print("\n=== PITCH-WEIGHTED 2026 MLB AVERAGE EYE POINT (Statcast feet) ===")
    rows = []
    for s in ("R", "L"):
        g = b[b.stand == s]
        w = g.pitches.values
        hbar = np.average(g.ht_in, weights=w)
        z = np.average(g.eye_z, weights=w)
        rows.append(dict(side=s, batter_sides=len(g), pitches=int(w.sum()),
                         mean_height_in=hbar,
                         eye_x_obm=lat_obm[s], eye_x_mlb=lat_mlb[s],
                         eye_y=dep[s], eye_z=z,
                         eye_z_unadjusted=per[per.side == s].eye_ht.mean()))
    summ = pd.DataFrame(rows)
    print(summ.round(3).to_string(index=False))
    summ.to_csv(OUT_SUMMARY, index=False)

    print("\n=== Effect of the height adjustment ===")
    for r in rows:
        d_in = (r["eye_z"] - r["eye_z_unadjusted"]) * 12
        print(f"  {r['side']}: OBM stature {per[per.side==r['side']].height_in.mean():.1f} in -> "
              f"MLB {r['mean_height_in']:.1f} in  |  eye height "
              f"{r['eye_z_unadjusted']:.3f} -> {r['eye_z']:.3f} ft ({d_in:+.1f} in)")

    print("\n=== Per-batter eye-height range across the 2026 population ===")
    for s in ("R", "L"):
        g = b[b.stand == s]
        print(f"  {s}: {g.eye_z.min():.2f} to {g.eye_z.max():.2f} ft "
              f"(stature {g.ht_in.min():.0f}-{g.ht_in.max():.0f} in), "
              f"within-batter uncertainty +/-{g.eye_z_sd.mean():.2f} ft")
    print(f"\nwrote {OUT_SUMMARY}")


if __name__ == "__main__":
    main()

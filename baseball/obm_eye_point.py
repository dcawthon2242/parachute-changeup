#!/usr/bin/env python3
"""
Extract the average left- and right-handed hitter EYE POSITION from the
OpenBiomechanics hitting C3D trials, expressed in the same coordinate frame and
units as this project's tunneling metric (Statcast feet, origin at the back tip
of home plate).

WHY A REGISTRATION IS NEEDED
----------------------------
The OBM hitting README documents a global frame with origin at the back tip of
home plate, +x toward the mound, +y toward the RHH box, +z up, and says distances
are in inches. Three of those are wrong or incomplete for the raw C3D:

  * The C3D POINT:UNITS parameter says METERS, not inches.
  * The lateral origin is NOT the plate centerline. Right-handed hitters stand at
    lab y = +0.5 m and left-handers at y = -0.95 m, which is not mirror-symmetric.
    Since the two batter's boxes ARE symmetric about the plate by rule, the
    midpoint of the two stances locates the true centerline -- it comes out near
    y = -0.22 m, an offset of about 9 inches.
  * The release notes say poi_x/poi_y/poi_z would give a documented-frame anchor,
    but this release's poi_metrics.csv does not contain those columns, so the
    frame has to be recovered from the data.

The axis DIRECTIONS in the README do check out empirically: lead-foot minus
rear-foot is positive in x (so +x is toward the mound), right-handed hitters sit
at positive y (so +y is the RHH box), and ankle joint centres sit 0.1-0.2 m above
z = 0 (so z = 0 is the turf).

WHAT IS RECOVERED HOW
---------------------
  lateral (y) : mirror-symmetry constraint between L and R stances. Tight.
  vertical (z): ground plane, z = 0. Tight.
  depth   (x) : weakest. Anchored by assuming the average stance sits at the
                lengthwise centre of the batter's box, which by rule is level
                with the middle of home plate, i.e. 8.5 in in front of the back
                tip. Reported with a sensitivity range.
"""

import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import ezc3d

ROOT = "data/obm/c3d/c3d"
M_PER_FT = 0.3048
HEAD_FRONT = ["LFHD", "RFHD"]
HEAD_ALL = ["LFHD", "RFHD", "LBHD", "RBHD"]
FEET = ["LANK", "RANK"]

# Home plate is 17 in deep; its lengthwise middle is 8.5 in in front of the back
# tip, and the batter's box is centred on that same line by rule.
PLATE_MID_FROM_TIP_M = (8.5 / 12) * M_PER_FT

FNAME_RE = re.compile(r"^(\d+)_(\d+)_(\d+)_(\d+)_([LR])_(\d+)_(\d+)\.c3d$")


def parse_name(path):
    m = FNAME_RE.match(os.path.basename(path))
    if not m:
        return None
    user, session, height_in, weight_lb, side, swing, ev = m.groups()
    return dict(user=user, session=session, height_in=int(height_in),
                weight_lb=int(weight_lb), side=side, swing=int(swing),
                session_swing=f"{int(session)}_{int(swing)}", path=path)


def stance_position(path, markers, stance_end_t):
    """Median marker position over the pre-movement stance window, in lab metres."""
    c = ezc3d.c3d(path)
    labels = [s.strip() for s in c["parameters"]["POINT"]["LABELS"]["value"]]
    rate = float(c["parameters"]["POINT"]["RATE"]["value"][0])
    data = c["data"]["points"]
    n = data.shape[2]

    # Stance window: trial start up to 0.3 s before front-foot contact. If that
    # event is unknown, fall back to the first 0.5 s.
    end_f = int((stance_end_t - 0.3) * rate) if np.isfinite(stance_end_t) else int(0.5 * rate)
    end_f = max(30, min(end_f, n))

    out = {}
    for m in markers:
        if m not in labels:
            out[m] = np.full(3, np.nan)
            continue
        s = data[:3, labels.index(m), :end_f]
        with np.errstate(invalid="ignore"):
            out[m] = np.nanmedian(np.where(s == 0.0, np.nan, s), axis=1)
    return out


def main():
    files = sorted(f for f in glob.glob(os.path.join(ROOT, "*", "*.c3d"))
                   if not f.endswith("_model.c3d"))
    meta = [parse_name(f) for f in files]
    meta = [m for m in meta if m]
    print(f"parsed {len(meta)} swing trials", file=sys.stderr)

    # front-foot-contact time per swing, for the stance window
    lm = pd.read_csv("data/obm/landmarks/landmarks.csv",
                     usecols=["session_swing", "fp_10_time"]).drop_duplicates("session_swing")
    fp = dict(zip(lm.session_swing, lm.fp_10_time))

    rows = []
    for i, m in enumerate(meta):
        if i % 100 == 0:
            print(f"  {i}/{len(meta)}", file=sys.stderr)
        try:
            pos = stance_position(m["path"], HEAD_ALL + FEET,
                                  fp.get(m["session_swing"], np.nan))
        except Exception as e:
            print(f"  FAILED {m['path']}: {e}", file=sys.stderr)
            continue
        eye = np.nanmean([pos["LFHD"], pos["RFHD"]], axis=0)
        head = np.nanmean([pos[k] for k in HEAD_ALL], axis=0)
        feet = np.nanmean([pos["LANK"], pos["RANK"]], axis=0)
        rows.append(dict(**{k: m[k] for k in
                            ("user", "session_swing", "side", "height_in", "weight_lb")},
                         eye_x=eye[0], eye_y=eye[1], eye_z=eye[2],
                         head_x=head[0], head_y=head[1], head_z=head[2],
                         feet_x=feet[0], feet_y=feet[1], feet_z=feet[2]))

    d = pd.DataFrame(rows).dropna(subset=["eye_x", "eye_y", "eye_z", "feet_y"])
    d.to_csv("data/obm/stance_raw_lab_metres.csv", index=False)
    print(f"\nusable trials: {len(d)}  "
          f"(R={(d.side=='R').sum()}, L={(d.side=='L').sum()}), "
          f"athletes={d.user.nunique()}\n")

    # ---- registration ------------------------------------------------------
    # Average per athlete first so prolific hitters don't dominate the constraint.
    per = d.groupby(["user", "side"], as_index=False)[
        ["eye_x", "eye_y", "eye_z", "head_x", "head_y", "head_z",
         "feet_x", "feet_y", "feet_z"]].mean()
    bySide = per.groupby("side")

    fy = bySide["feet_y"].mean()
    y0 = (fy["R"] + fy["L"]) / 2.0            # lab y of the plate centreline
    fx = bySide["feet_x"].mean()
    x0 = fx.mean() - PLATE_MID_FROM_TIP_M     # lab x of the back tip of the plate

    print("=== RECOVERED REGISTRATION (lab metres) ===")
    print(f"  stance feet lateral : R {fy['R']:+.3f}  L {fy['L']:+.3f}")
    print(f"  -> plate centreline  y0 = {y0:+.3f} m  ({y0/M_PER_FT*12:+.1f} in off the documented origin)")
    print(f"  symmetry residual   : R {fy['R']-y0:+.3f} vs L {-(fy['L']-y0):+.3f} m "
          f"(should match; diff {abs((fy['R']-y0)+(fy['L']-y0)):.4f})")
    print(f"  stance feet depth   : R {fx['R']:+.3f}  L {fx['L']:+.3f}  mean {fx.mean():+.3f}")
    print(f"  -> plate back tip    x0 = {x0:+.3f} m")
    print(f"  vertical            : z = 0 is the turf (ankle JC sits ~0.1 m up)\n")

    def to_statcast(px, py, pz):
        """lab metres -> Statcast feet: +x catcher's right, +y toward pitcher, +z up."""
        return (-(py - y0) / M_PER_FT, (px - x0) / M_PER_FT, pz / M_PER_FT)

    print("=== AVERAGE HITTER POSITION IN TUNNEL-METRIC COORDINATES (Statcast feet) ===")
    print(f"{'point':12s} {'side':5s} {'n':>4s} {'x (lateral)':>12s} {'y (depth)':>11s} {'z (height)':>11s}")
    res = {}
    for label, cols in (("eyes", ("eye_x", "eye_y", "eye_z")),
                        ("head centre", ("head_x", "head_y", "head_z")),
                        ("feet", ("feet_x", "feet_y", "feet_z"))):
        for side in ("R", "L"):
            g = per[per.side == side]
            sx, sy, sz = to_statcast(g[cols[0]].mean(), g[cols[1]].mean(), g[cols[2]].mean())
            n = len(g)
            print(f"{label:12s} {side:5s} {n:4d} {sx:12.3f} {sy:11.3f} {sz:11.3f}")
            res[(label, side)] = (sx, sy, sz)
        print()

    # spread across athletes, for the eye point
    print("=== ATHLETE-TO-ATHLETE SPREAD OF THE EYE POINT (Statcast feet, SD) ===")
    for side in ("R", "L"):
        g = per[per.side == side]
        sx, sy, sz = to_statcast(g.eye_x.values, g.eye_y.values, g.eye_z.values)
        print(f"  {side}: x {np.mean(sx):+.3f} +/- {np.std(sx):.3f} | "
              f"y {np.mean(sy):+.3f} +/- {np.std(sy):.3f} | "
              f"z {np.mean(sz):+.3f} +/- {np.std(sz):.3f}")

    # depth-anchor sensitivity: the x origin is the weak one
    print("\n=== SENSITIVITY OF THE DEPTH COORDINATE TO THE x ANCHOR ===")
    print("  (if the average stance is not centred in the box, y shifts by the same amount)")
    for shift_in in (-6, -3, 0, 3, 6):
        sh = shift_in / 12.0
        print(f"    stance {shift_in:+3d} in from box centre -> "
              f"eye depth R {res[('eyes','R')][1] + sh:+.3f} ft, "
              f"L {res[('eyes','L')][1] + sh:+.3f} ft")

    out = []
    for (label, side), (sx, sy, sz) in res.items():
        out.append(dict(point=label, side=side, x_ft=sx, y_ft=sy, z_ft=sz))
    pd.DataFrame(out).to_csv("data/obm/hitter_reference_points_statcast_ft.csv", index=False)
    print("\nwrote data/obm/hitter_reference_points_statcast_ft.csv")


if __name__ == "__main__":
    main()

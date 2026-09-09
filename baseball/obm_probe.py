#!/usr/bin/env python3
"""Probe the OpenBiomechanics hitting C3D files: marker labels, units, frame rate,
events, and the orientation of the global frame. Read-only reconnaissance before
building the eye-position extraction."""

import glob
import os
import numpy as np
import ezc3d

ROOT = "data/obm/c3d/c3d"

files = sorted(f for f in glob.glob(os.path.join(ROOT, "*", "*.c3d"))
               if not f.endswith("_model.c3d"))
print(f"swing trials: {len(files)}")
print(f"model files : {len(glob.glob(os.path.join(ROOT, '*', '*_model.c3d')))}\n")

f = files[0]
print("=== probing:", f, "===")
c = ezc3d.c3d(f)

pts = c["parameters"]["POINT"]
labels = [s.strip() for s in pts["LABELS"]["value"]]
print("units      :", pts["UNITS"]["value"])
print("rate       :", pts["RATE"]["value"], "Hz")
print("frames     :", c["data"]["points"].shape[2])
print("n markers  :", len(labels))
print("\nall labels:")
print(labels)

head = [l for l in labels if "HD" in l.upper() or "HEAD" in l.upper()]
print("\nhead markers found:", head)

data = c["data"]["points"]  # 4 x n_markers x n_frames
def series(name):
    return data[:3, labels.index(name), :]

for m in head + ["RANK", "LANK", "CLAV", "RSHO", "LSHO"]:
    if m in labels:
        s = series(m)
        ok = np.isfinite(s).all(axis=0)
        print(f"  {m:6s} valid {ok.sum():5d}/{s.shape[1]:5d}  "
              f"mean xyz = ({np.nanmean(s[0][ok]):8.1f}, {np.nanmean(s[1][ok]):8.1f}, {np.nanmean(s[2][ok]):8.1f})")

ev = c["parameters"].get("EVENT", {})
if ev and "LABELS" in ev:
    print("\nevents:", [s.strip() for s in ev["LABELS"]["value"]])
    if "TIMES" in ev:
        print("times :", np.round(ev["TIMES"]["value"], 3).tolist())
else:
    print("\nno EVENT parameter block")

# Orientation sanity: ankles should straddle the plate laterally, head well above ground.
print("\n=== orientation check across a few L and R trials ===")
print(f"{'file':52s} {'side':5s} {'ank_x':>8s} {'ank_y':>8s} {'ank_z':>8s} {'hd_z':>8s}")
for g in files[:400:37]:
    side = os.path.basename(g).split("_")[4]
    try:
        cc = ezc3d.c3d(g)
    except Exception as e:
        print(f"{os.path.basename(g)[:50]:52s} FAILED {e}")
        continue
    ll = [s.strip() for s in cc["parameters"]["POINT"]["LABELS"]["value"]]
    dd = cc["data"]["points"]
    def m(nm):
        if nm not in ll:
            return np.full(3, np.nan)
        s = dd[:3, ll.index(nm), :]
        n = min(60, s.shape[1])
        return np.nanmean(s[:, :n], axis=1)
    ank = np.nanmean([m("RANK"), m("LANK")], axis=0)
    hd = np.nanmean([m("LFHD"), m("RFHD")], axis=0)
    print(f"{os.path.basename(g)[:50]:52s} {side:5s} "
          f"{ank[0]:8.1f} {ank[1]:8.1f} {ank[2]:8.1f} {hd[2]:8.1f}")

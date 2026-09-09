#!/usr/bin/env python3
"""Pull Savant broadcast clips of a pitcher's changeup whiffs and write GIFs.

Companion to docs/parachute-changeup-README.md. The README argues that a
changeup spinning on the four-seamer's axis is deceptive because it is a
*pair* property; the fastest way to see that is the broadcast clip, frozen at
release, next to the same pitcher's four-seamer from the same game.

Pipeline (all public, unofficial Savant endpoints; may change):
  1. data/statcast_<season>/statcast_<season>_all.csv  -> pick the pitcher's
     regular-season changeup swinging strikes.
  2. data/rubber/play_ids_<season>.csv (from rubber_01_playids.py) -> map
     (game_pk, at_bat_number, pitch_number) to Savant playId.
  3. https://baseballsavant.mlb.com/sporty-videos?playId=... -> scrape mp4 URL.
  4. ffmpeg -> GIF (and, with --pair, a side-by-side with a four-seam swinging
     strike from the same game so release frames can be compared).

Usage:
    .venv-cv/bin/python baseball/parachute_clips.py --pitcher 656302 --season 2026 --n 3 --pair
    .venv-cv/bin/python baseball/parachute_clips.py --pitcher 669373 --season 2021 --n 2

Selection: changeup swinging strikes (incl. blocked and foul tips), ordered by
velocity separation from the pitcher's season-mean four-seam velocity, largest
first, so the clips shown are the ones where the "same look, ten-plus off"
mechanism is most visible. Use --order recent for chronological-latest instead.

Output: data/statcast_model/article_assets/clips/<pitcher>_<season>/...
        plus clips_manifest.csv in that folder (one row per clip, with mp4 URL).
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from statistics import mean

import requests

REPO_ROOT = Path(__file__).resolve().parents[1]
STATCAST_DIR = REPO_ROOT / "data"
PLAYID_DIR = REPO_ROOT / "data" / "rubber"
OUT_ROOT = REPO_ROOT / "data" / "statcast_model" / "article_assets" / "clips"

SPORTY_PAGE = "https://baseballsavant.mlb.com/sporty-videos?playId={play_id}"
MP4_RE = re.compile(r"https://sporty-clips\.mlb\.com/[A-Za-z0-9_\-=]+\.mp4")
USER_AGENT = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0 Safari/537.36"
)
REFERER = "https://baseballsavant.mlb.com/"  # sporty-clips 403s without it

WHIFF = {"swinging_strike", "swinging_strike_blocked", "foul_tip", "missed_bunt"}
GIF_FPS = 15
GIF_WIDTH = 480
CLIP_SECONDS = 7.0  # Savant clips run ~6.7 s


def ffmpeg_exe() -> str:
    exe = shutil.which("ffmpeg")
    if exe:
        return exe
    try:
        import imageio_ffmpeg  # type: ignore

        return imageio_ffmpeg.get_ffmpeg_exe()
    except Exception:  # noqa: BLE001
        sys.exit("ffmpeg not found: install ffmpeg or pip install imageio-ffmpeg")


def load_pitches(season: int, pitcher: int) -> list[dict]:
    path = STATCAST_DIR / f"statcast_{season}" / f"statcast_{season}_all.csv"
    if not path.exists():
        sys.exit(f"missing {path}; run baseball/scrape_statcast_multi.R {season}")
    rows = []
    with path.open(newline="") as fh:
        for r in csv.DictReader(fh):
            if r.get("pitcher") != str(pitcher) or r.get("game_type") != "R":
                continue
            rows.append(r)
    if not rows:
        sys.exit(f"no regular-season pitches for pitcher {pitcher} in {season}")
    return rows


def load_play_ids(season: int) -> dict[tuple[str, str, str], str]:
    path = PLAYID_DIR / f"play_ids_{season}.csv"
    if not path.exists():
        sys.exit(
            f"missing {path}; run .venv-cv/bin/python baseball/rubber_01_playids.py "
            f"--season {season}"
        )
    out = {}
    with path.open(newline="") as fh:
        for r in csv.DictReader(fh):
            out[(r["game_pk"], r["at_bat_number"], r["pitch_number"])] = r["play_id"]
    return out


def key(r: dict) -> tuple[str, str, str]:
    return (r["game_pk"], r["at_bat_number"], r["pitch_number"])


def fnum(x: str) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def resolve_mp4(session: requests.Session, play_id: str, retries: int = 4) -> str:
    url = SPORTY_PAGE.format(play_id=play_id)
    for attempt in range(1, retries + 1):
        try:
            resp = session.get(url, timeout=30)
            resp.raise_for_status()
            m = MP4_RE.search(resp.text)
            return m.group(0) if m else ""
        except Exception:  # noqa: BLE001
            if attempt == retries:
                return ""
            time.sleep(min(20.0, 2.0**attempt))
    return ""


def headers_arg() -> str:
    return f"Referer: {REFERER}\r\nUser-Agent: {USER_AGENT}\r\n"


def run(cmd: list[str], timeout: int = 180) -> None:
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr[-400:])


def download_mp4(ff: str, mp4_url: str, dst: Path) -> None:
    run([ff, "-y", "-loglevel", "error", "-headers", headers_arg(), "-i", mp4_url,
         "-t", str(CLIP_SECONDS), "-c", "copy", str(dst)])


def to_gif(ff: str, src: Path, dst: Path) -> None:
    # Two-pass palette for a clean GIF at a reasonable size.
    vf = f"fps={GIF_FPS},scale={GIF_WIDTH}:-1:flags=lanczos"
    run([ff, "-y", "-loglevel", "error", "-i", str(src),
         "-vf", f"{vf},split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse",
         str(dst)])


def side_by_side(ff: str, left: Path, right: Path, dst: Path) -> None:
    vf = (
        f"[0:v]fps={GIF_FPS},scale={GIF_WIDTH}:-1:flags=lanczos[l];"
        f"[1:v]fps={GIF_FPS},scale={GIF_WIDTH}:-1:flags=lanczos[r];"
        f"[l][r]hstack=inputs=2,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse"
    )
    run([ff, "-y", "-loglevel", "error", "-i", str(left), "-i", str(right),
         "-filter_complex", vf, "-shortest", str(dst)])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pitcher", type=int, required=True, help="MLBAM id")
    ap.add_argument("--season", type=int, required=True)
    ap.add_argument("--n", type=int, default=3, help="number of changeup whiffs")
    ap.add_argument("--order", choices=["sep", "recent"], default="sep")
    ap.add_argument("--pair", action="store_true",
                    help="also fetch a four-seam swinging strike from the same game "
                         "and write a side-by-side GIF")
    ap.add_argument("--keep-mp4", action="store_true")
    args = ap.parse_args()

    ff = ffmpeg_exe()
    pitches = load_pitches(args.season, args.pitcher)
    play_ids = load_play_ids(args.season)

    ff_velo = [fnum(r["release_speed"]) for r in pitches if r["pitch_type"] == "FF"]
    ff_velo = [v for v in ff_velo if v == v]
    if len(ff_velo) < 50:
        sys.exit("fewer than 50 four-seamers: the parachute spec anchors on FF only")
    fb_mean = mean(ff_velo)

    ch_whiffs = [r for r in pitches
                 if r["pitch_type"] == "CH" and r["description"] in WHIFF
                 and key(r) in play_ids]
    if not ch_whiffs:
        sys.exit("no changeup whiffs with a playId; is play_ids_<season>.csv complete?")
    for r in ch_whiffs:
        r["_sep"] = fb_mean - fnum(r["release_speed"])
    if args.order == "sep":
        ch_whiffs.sort(key=lambda r: -r["_sep"])
    else:
        ch_whiffs.sort(key=lambda r: (r["game_date"], int(r["at_bat_number"]),
                                      int(r["pitch_number"])), reverse=True)
    chosen = ch_whiffs[: args.n]

    ff_by_game: dict[str, list[dict]] = {}
    if args.pair:
        for r in pitches:
            if r["pitch_type"] == "FF" and r["description"] in WHIFF and key(r) in play_ids:
                ff_by_game.setdefault(r["game_pk"], []).append(r)

    name = chosen[0]["player_name"].replace(", ", "_").replace(" ", "")
    out_dir = OUT_ROOT / f"{args.pitcher}_{args.season}_{name}"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = out_dir / "clips_manifest.csv"

    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT, "Referer": REFERER})

    rows_out = []
    for i, r in enumerate(chosen, 1):
        pid = play_ids[key(r)]
        stem = f"{i:02d}_CH_{r['game_date']}_ab{r['at_bat_number']}_p{r['pitch_number']}"
        mp4_url = resolve_mp4(session, pid)
        rec = {"clip": stem, "pitch_type": "CH", "game_date": r["game_date"],
               "game_pk": r["game_pk"], "at_bat_number": r["at_bat_number"],
               "pitch_number": r["pitch_number"], "description": r["description"],
               "release_speed": r["release_speed"], "velo_sep": f"{r['_sep']:.1f}",
               "stand": r["stand"], "play_id": pid, "mp4_url": mp4_url, "status": ""}
        if not mp4_url:
            rec["status"] = "no_mp4"
            rows_out.append(rec)
            print(f"[{i}] {stem}: no mp4", file=sys.stderr)
            continue
        mp4 = out_dir / f"{stem}.mp4"
        gif = out_dir / f"{stem}.gif"
        try:
            download_mp4(ff, mp4_url, mp4)
            to_gif(ff, mp4, gif)
            rec["status"] = "ok"
            print(f"[{i}] {gif.name}  sep {r['_sep']:.1f} mph  {r['description']}")
        except Exception as e:  # noqa: BLE001
            rec["status"] = f"ffmpeg: {e}"
            rows_out.append(rec)
            print(f"[{i}] {stem}: {e}", file=sys.stderr)
            continue

        if args.pair:
            cands = ff_by_game.get(r["game_pk"], [])
            if cands:
                f = cands[0]
                fpid = play_ids[key(f)]
                fstem = f"{i:02d}_FF_{f['game_date']}_ab{f['at_bat_number']}_p{f['pitch_number']}"
                fmp4_url = resolve_mp4(session, fpid)
                if fmp4_url:
                    fmp4 = out_dir / f"{fstem}.mp4"
                    pair = out_dir / f"{i:02d}_PAIR_FF_left_CH_right.gif"
                    try:
                        download_mp4(ff, fmp4_url, fmp4)
                        side_by_side(ff, fmp4, mp4, pair)
                        rec["pair_gif"] = pair.name
                        rec["pair_ff_play_id"] = fpid
                        print(f"     paired with {fstem} ({f['release_speed']} mph) -> {pair.name}")
                    except Exception as e:  # noqa: BLE001
                        rec["pair_gif"] = f"ffmpeg: {e}"
                    if not args.keep_mp4:
                        fmp4.unlink(missing_ok=True)
            else:
                rec["pair_gif"] = "no FF whiff with playId in this game"
        if not args.keep_mp4:
            mp4.unlink(missing_ok=True)
        rows_out.append(rec)
        time.sleep(1.0)  # courtesy pause between Savant fetches

    fields = sorted({k for r in rows_out for k in r})
    with manifest.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows_out)
    print(f"\nwrote {len(rows_out)} rows to {manifest}")


if __name__ == "__main__":
    main()

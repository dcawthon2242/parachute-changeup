#!/usr/bin/env python3
"""Harvest D1 baseball rosters from school athletics sites, verified against TrackMan.

This is the fallback route for when stats.ncaa.org is rate limiting. School athletics sites have no
shared bot protection - plain HTTP, no browser, no challenge - and most run Sidearm, which puts the
roster in a predictable table where height carries an inches value in data-sort.

The catch is that each school needs its own domain, and domains are not guessable: Wake Forest is
godeacs.com, LSU is lsusports.net, Arizona State is thesundevils.com. Candidates are therefore
gathered liberally from a seed list, name patterns and Wikipedia's external links, and the accuracy
comes from verification rather than from trusting the source.

Verification is the important part. Checking only that a domain serves *a* baseball roster is not
enough - an early version happily assigned floridagators.com to Arkansas State, which would have
silently attached the wrong heights to real pitchers. Instead each fetched roster is matched against
the TrackMan pitcher lists by name overlap: the correct school shares most of a team's pitchers and
every other school shares almost none, so a wrong domain scores near zero and is discarded. That
also identifies which TrackMan team code the roster belongs to, which is needed downstream anyway.

Output matches the NCAA scraper's JSONL schema, so build_d1_heights.R consumes either or both.
"""

import argparse, csv, html, json, os, re, ssl, threading, time
import urllib.error, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor

UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/126.0 Safari/537.36")
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data", "ncaa_rosters")
OUT = os.path.join(DATA, "sidearm_rosters.jsonl")
SITES = os.path.join(DATA, "sidearm_sites.csv")
YEARS = [2025, 2024, 2023]
lock = threading.Lock()


def fetch(url, timeout=20):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "text/html"})
    with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
        return r.getcode(), r.read().decode("utf-8", "ignore")


def txt(s):
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", s))).strip()


def parse_nuxt(page):
    """Newer Sidearm sites render through Nuxt and ship the roster as a flat, index-encoded payload.

    Object values are offsets into one big array rather than literals, so {"heightFeet": 210} means
    element 210 holds the number. Dereferencing is a one-liner and gives exact feet/inches, which is
    better than scraping the rendered markup where height is only formatted text.
    """
    m = re.search(r'<script[^>]*type="application/json"[^>]*>(.*?)</script>', page, re.S)
    if not m:
        return []
    try:
        d = json.loads(m.group(1))
    except Exception:
        return []
    if not isinstance(d, list):
        return []

    def deref(v):
        return d[v] if isinstance(v, int) and 0 <= v < len(d) else v

    out = []
    for p in d:
        if not (isinstance(p, dict) and "heightFeet" in p and "lastName" in p):
            continue
        fn, ln = deref(p.get("firstName")), deref(p.get("lastName"))
        if not (isinstance(fn, str) and isinstance(ln, str)):
            continue
        ft, inch = deref(p.get("heightFeet")), deref(p.get("heightInches"))
        inches = None
        if isinstance(ft, int) and isinstance(inch, int) and 4 <= ft <= 7:
            inches = ft * 12 + inch
            if not (55 <= inches <= 90):
                inches = None
        pos = deref(p.get("positionShort"))
        out.append({"name": f"{fn} {ln}".strip(), "pos": pos if isinstance(pos, str) else "",
                    "height_in": inches, "height_raw": "", "throws": "", "bats": "",
                    "hometown": ""})
    return out


def parse_roster(page):
    """Sidearm roster table: name in the roster anchor, height in a td carrying data-sort inches."""
    out = []
    for row in re.findall(r"<tr[^>]*>(.*?)</tr>", page, re.S | re.I):
        if "height" not in row:
            continue
        nm = re.search(r'<a [^>]*href="[^"]*roster/[^"]*"[^>]*>(.*?)</a>', row, re.S | re.I)
        if not nm:
            continue
        name = txt(nm.group(1))
        if not name or len(name) < 4 or len(name) > 45:
            continue
        h = re.search(r'class="height[^"]*"[^>]*data-sort="(\d{2})"', row)
        inches = int(h.group(1)) if h else None
        if inches is None:
            m2 = re.search(r'class="height[^"]*"[^>]*>\s*([5-7])[-\'\u2019]\s*(\d{1,2})', row)
            if m2:
                inches = int(m2.group(1)) * 12 + int(m2.group(2))
        if inches is not None and not (55 <= inches <= 90):
            inches = None
        pos = re.search(r'class="rp_position_short[^"]*"[^>]*>(.*?)</td>', row, re.S)
        out.append({"name": name, "pos": txt(pos.group(1)) if pos else "",
                    "height_in": inches, "height_raw": "", "throws": "", "bats": "",
                    "hometown": ""})
    return out


def name_key(x):
    toks = [t for t in re.split(r"[^a-zA-Z]+", x.lower()) if len(t) > 1]
    return "|".join(sorted(toks))


def load_tm():
    """TrackMan pitcher name keys by (team code, year), used to verify and identify a roster."""
    tm = {}
    with open(os.path.join(DATA, "tm_pitchers.csv")) as f:
        for r in csv.DictReader(f):
            tm.setdefault((r["PitcherTeam"], int(r["year"])), set()).add(r["k"])
    return tm


def best_match(keys, tm, year):
    """Which TrackMan team does this roster look like? Returns (code, overlap, margin)."""
    scored = []
    for (code, yr), pit in tm.items():
        if yr != year or len(pit) < 6:
            continue
        ov = len(pit & keys) / len(pit)
        scored.append((ov, code))
    if not scored:
        return None, 0.0, 0.0
    scored.sort(reverse=True)
    top = scored[0]
    second = scored[1][0] if len(scored) > 1 else 0.0
    return top[1], top[0], top[0] - second


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--min-overlap", type=float, default=0.35)
    a = ap.parse_args()

    tm = load_tm()
    print(f"loaded {len(tm)} TrackMan team-seasons for verification", flush=True)

    with open(SITES) as f:
        sites = [r for r in csv.DictReader(f) if r["domain"]]
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT):
            try:
                r = json.loads(l)
                done.add((r["school"], r["year"]))
            except Exception:
                pass
    jobs = [(s["school"], s["domain"], y) for s in sites for y in YEARS
            if (s["school"], y) not in done]
    if a.limit:
        jobs = jobs[:a.limit]
    print(f"{len(sites)} domains x {len(YEARS)} years -> {len(jobs)} pages to fetch", flush=True)

    n_ok = n_rej = n_miss = 0

    def work(job):
        nonlocal n_ok, n_rej, n_miss
        school, dom, year = job
        players = []
        for path in (f"/sports/baseball/roster/{year}",
                     f"/sports/baseball/roster/{year-1}-{str(year)[2:]}"):
            try:
                code, body = fetch(f"https://{dom}{path}")
            except Exception:
                continue
            if code == 200:
                players = parse_roster(body) or parse_nuxt(body)
                if players:
                    break
        if not players:
            with lock:
                n_miss += 1
            return None
        keys = {name_key(p["name"]) for p in players}
        tcode, ov, margin = best_match(keys, tm, year)
        rec = {"school": school, "year": year, "domain": dom, "source": "sidearm",
               "trackman_team": tcode, "overlap": round(ov, 3), "margin": round(margin, 3),
               "players": players}
        with lock:
            if ov < a.min_overlap:
                n_rej += 1
                print(f"  REJECT {school:24} {year}  {dom:26} overlap {ov:.2f} "
                      f"(best guess {tcode})", flush=True)
                return None
            n_ok += 1
            print(f"  ok     {school:24} {year}  {tcode:9} overlap {ov:.2f} "
                  f"margin {margin:.2f}  {len(players)}p", flush=True)
        return rec

    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        for rec in ex.map(work, jobs):
            if rec:
                with lock, open(OUT, "a") as f:
                    f.write(json.dumps(rec) + "\n")
    print(f"\nverified {n_ok}, rejected {n_rej} (wrong school), no roster found {n_miss}")
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()

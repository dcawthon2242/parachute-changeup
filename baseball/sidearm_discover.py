#!/usr/bin/env python3
"""Find each school's athletics domain and the path that serves its baseball roster.

The NCAA route is rate limited to the point of being unusable, but school athletics sites have no
shared bot protection - a probe of six of them returned full rosters over plain HTTP with no browser
and no challenge. The one thing that route needs and the NCAA route did not is a domain per school,
which is the whole reason it was the fallback.

Guessing domains from the school name does not work reliably enough to trust (Wake Forest is
godeacs.com, LSU is lsusports.net), so candidates are gathered from Wikipedia's external links,
where the real domain appears even when the specific link is a dead 2008 article, and topped up
with pattern guesses. Nothing is trusted on the strength of where it came from: every candidate is
confirmed by actually fetching a roster path and finding heights in the response, so a wrong guess
fails loudly instead of silently producing an empty roster.

Writes data/ncaa_rosters/sidearm_sites.csv with school, domain, path and the height count seen.
"""

import argparse, csv, json, os, re, ssl, sys, threading, time
import urllib.error, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor

UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/126.0 Safari/537.36")
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data", "ncaa_rosters")
OUT = os.path.join(DATA, "sidearm_sites.csv")

PATHS = ["/sports/baseball/roster/2025", "/sports/baseball/roster/2024-25",
         "/sports/baseball/roster", "/sports/mbase/2024-25/roster",
         "/sports/baseball/roster/season/2025"]
JUNK = re.compile(r"wikipedia|wikimedia|wikidata|archive\.org|doi\.org|jstor|nytimes|espn\.com|"
                  r"twitter|facebook|instagram|youtube|google|amazon|apple\.com|si\.com|"
                  r"sports-reference|census\.gov|loc\.gov|worldcat|isbn|nih\.gov|\.pdf$")
lock = threading.Lock()


def fetch(url, timeout=20):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "text/html"})
    with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
        return r.getcode(), r.read().decode("utf-8", "ignore")


def wiki_json(params):
    u = "https://en.wikipedia.org/w/api.php?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(u, headers={"User-Agent": "roster-research/1.0"})
    return json.load(urllib.request.urlopen(req, timeout=30))


def wiki_domains(school):
    """Domains appearing in the external links of the school's athletics/baseball articles."""
    doms = {}
    try:
        srch = wiki_json({"action": "query", "list": "search", "format": "json",
                          "formatversion": "2", "srlimit": "3",
                          "srsearch": f"{school} college baseball team"})
        titles = [s["title"] for s in srch.get("query", {}).get("search", [])]
    except Exception:
        titles = []
    for t in titles[:3]:
        try:
            pg = wiki_json({"action": "query", "prop": "extlinks", "ellimit": "500",
                            "titles": t, "format": "json", "formatversion": "2"})
            links = [e["url"] for e in pg["query"]["pages"][0].get("extlinks", [])]
        except Exception:
            continue
        for u in links:
            if JUNK.search(u):
                continue
            m = re.match(r"https?://([^/]+)", u)
            if not m:
                continue
            d = m.group(1).lower().lstrip("www.")
            if d.endswith(".cstv.com"):        # legacy host, but the slug names the school
                continue
            doms[d] = doms.get(d, 0) + 1
    return [d for d, _ in sorted(doms.items(), key=lambda x: -x[1])][:6]


def guesses(school, nickname):
    s = re.sub(r"[^a-z]", "", school.lower())
    n = re.sub(r"[^a-z]", "", (nickname or "").lower())
    out = []
    for base in filter(None, [f"{s}{n}", f"go{n}", f"{n}sports", f"{s}sports", f"go{s}",
                              f"{s}athletics", s]):
        out += [f"{base}.com", f"{base}.net"]
    return out[:10]


def probe(domain):
    """Return (path, height_count) if this domain serves a baseball roster, else None."""
    for p in PATHS:
        try:
            code, body = fetch(f"https://{domain}{p}")
        except Exception:
            continue
        if code != 200:
            continue
        n = len(re.findall(r"\b([5-7])[-'\u2019]\s?(\d{1,2})\b", body))
        n += len(re.findall(r"(?i)\bheight\b", body))
        if n >= 10 and re.search(r"(?i)roster", body):
            return p, n
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--schools", default="")
    a = ap.parse_args()

    with open(os.path.join(DATA, "d1_teams.csv")) as f:
        schools = sorted({r["team_name"] for r in csv.DictReader(f)})
    if a.schools:
        want = {s.strip() for s in a.schools.split(",")}
        schools = [s for s in schools if s in want]
    done = {}
    if os.path.exists(OUT):
        with open(OUT) as f:
            done = {r["school"]: r for r in csv.DictReader(f)}
    todo = [s for s in schools if s not in done]
    if a.limit:
        todo = todo[:a.limit]
    print(f"{len(schools)} schools, {len(done)} already resolved, {len(todo)} to try", flush=True)

    results, n_ok = [], 0

    def work(school):
        nonlocal n_ok
        cands, seen = [], set()
        for d in wiki_domains(school) + guesses(school, ""):
            if d not in seen:
                seen.add(d)
                cands.append(d)
        for d in cands[:12]:
            r = probe(d)
            if r:
                with lock:
                    n_ok += 1
                    print(f"  OK   {school:26} -> {d}{r[0]}", flush=True)
                return {"school": school, "domain": d, "path": r[0], "signal": r[1]}
        with lock:
            print(f"  MISS {school:26} (tried {len(cands[:12])})", flush=True)
        return {"school": school, "domain": "", "path": "", "signal": 0}

    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        results = list(ex.map(work, todo))

    allr = list(done.values()) + results
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["school", "domain", "path", "signal"])
        w.writeheader()
        for r in allr:
            w.writerow({k: r.get(k, "") for k in ["school", "domain", "path", "signal"]})
    hit = sum(1 for r in allr if r.get("domain"))
    print(f"\nresolved {hit}/{len(allr)} schools -> {OUT}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Scrape NCAA Division I baseball rosters, including listed height.

The roster page at stats.ncaa.org/teams/{season_team_id}/roster carries Name, Class, Position,
Height, Bats, Throws and Hometown in one table, so a single request per team-season is enough.
That matters: routing through individual player pages, which is what the player_url column invites,
would be roughly 13,000 requests against a protected host instead of 917.

Two things about the protection shape this design. Akamai rejects plain HTTP clients outright, so
a real browser is required. And it rate-blocks bursts: spawning a fresh Chrome per page makes every
request look like a first-time visitor that has to re-earn clearance, which starts returning
"Access Denied" within a handful of pages. So this drives ONE long-lived Chrome over the DevTools
protocol and navigates it from page to page, which keeps the clearance cookie and is also far
faster than paying process startup 917 times.

Results append to a JSONL checkpoint keyed by season_team_id, so an interrupted run resumes and a
repeat run only fetches what is missing. A block triggers exponential backoff rather than a crash.

    python3 baseball/scrape_ncaa_rosters.py [--limit N] [--retry-failed] [--delay 1.5]
"""

import argparse, atexit, html, json, os, random, re, shutil, subprocess, sys, time
import urllib.request

import websocket

CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/126.0 Safari/537.36")
PORT, PROFILE = 9333, "/tmp/ncaa_scrape_profile"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data", "ncaa_rosters")
TEAMS, OUT = os.path.join(DATA, "d1_teams.csv"), os.path.join(DATA, "rosters.jsonl")


class Browser:
    """A single headless Chrome, driven over CDP, that keeps its cookies between navigations."""

    def __init__(self):
        os.makedirs(PROFILE, exist_ok=True)
        self.proc = subprocess.Popen(
            [CHROME, "--headless=new", "--disable-gpu", "--no-sandbox", "--mute-audio",
             "--no-first-run", "--no-default-browser-check", "--disable-dev-shm-usage",
             "--disable-background-networking", "--disable-extensions",
             f"--remote-debugging-port={PORT}", "--remote-allow-origins=*",
             f"--user-data-dir={PROFILE}",
             f"--user-agent={UA}", "about:blank"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        atexit.register(self.close)
        ws_url = None
        for _ in range(60):
            try:
                tabs = json.load(urllib.request.urlopen(
                    f"http://127.0.0.1:{PORT}/json/list", timeout=2))
                pages = [t for t in tabs if t.get("type") == "page"]
                if pages:
                    ws_url = pages[0]["webSocketDebuggerUrl"]
                    break
            except Exception:
                pass
            time.sleep(0.5)
        if not ws_url:
            raise RuntimeError("could not attach to Chrome over CDP")
        self.ws = websocket.create_connection(ws_url, timeout=60)
        self.n = 0

    def cmd(self, method, **params):
        self.n += 1
        self.ws.send(json.dumps({"id": self.n, "method": method, "params": params}))
        while True:
            msg = json.loads(self.ws.recv())
            if msg.get("id") == self.n:
                return msg.get("result", {})

    def get(self, url, wait=25):
        self.cmd("Page.navigate", url=url)
        deadline = time.time() + wait
        dom = ""
        while time.time() < deadline:
            time.sleep(0.4)
            r = self.cmd("Runtime.evaluate",
                         expression="document.readyState + '|' + document.documentElement.outerHTML",
                         returnByValue=True)
            val = r.get("result", {}).get("value") or ""
            state, _, dom = val.partition("|")
            if state == "complete" and ("Height" in dom or "Access Denied" in dom):
                return dom
        return dom

    def close(self):
        try:
            self.ws.close()
        except Exception:
            pass
        try:
            self.proc.terminate()
            self.proc.wait(timeout=10)
        except Exception:
            try:
                self.proc.kill()
            except Exception:
                pass


def strip(x):
    return html.unescape(re.sub(r"<[^>]+>", " ", x)).replace("\xa0", " ").strip()


def parse(page):
    """Pull the roster table. Column order is read from the header, not assumed."""
    heads = [strip(x) for x in re.findall(r"<th[^>]*>(.*?)</th>", page, re.S | re.I)]
    seen, cols = set(), []
    for h in heads:                       # the header block repeats in thead and tfoot
        if h in seen:
            break
        seen.add(h)
        cols.append(h)
    if "Height" not in cols:
        return None
    rows = []
    for row in re.findall(r"<tr[^>]*>(.*?)</tr>", page, re.S | re.I):
        cells = [strip(x) for x in re.findall(r"<td[^>]*>(.*?)</td>", row, re.S | re.I)]
        if len(cells) >= len(cols):
            rows.append(dict(zip(cols, cells[:len(cols)])))
    return rows


def height_to_inches(s):
    if not s:
        return None
    m = re.match(r"^\s*(\d)\s*[-'\u2019]\s*(\d{1,2})", s)      # 6-2, 6'2", 6’2
    if m:
        v = int(m.group(1)) * 12 + int(m.group(2))
        return v if 55 <= v <= 90 else None
    m = re.match(r"^\s*(\d{2})\s*$", s)                         # already inches
    return int(m.group(1)) if m and 55 <= int(m.group(1)) <= 90 else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--delay", type=float, default=4.0)
    ap.add_argument("--retry-failed", action="store_true")
    ap.add_argument("--fresh-profile", action="store_true")
    ap.add_argument("--teams", default=TEAMS,
                    help="team list to work from; defaults to all 917 team-seasons")
    a = ap.parse_args()

    if not os.path.exists(CHROME):
        sys.exit(f"Chrome not found at {CHROME}")
    if a.fresh_profile:
        shutil.rmtree(PROFILE, ignore_errors=True)

    if a.retry_failed and os.path.exists(OUT):
        keep = [l for l in open(OUT) if "error" not in json.loads(l)]
        with open(OUT, "w") as f:
            f.writelines(keep)
        print(f"dropped failed rows, kept {len(keep)} good records")

    import csv
    with open(a.teams) as f:
        teams = list(csv.DictReader(f))
    done = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            try:
                done.add(json.loads(line)["season_team_id"])
            except Exception:
                pass
    todo = [t for t in teams if t["season_team_id"] not in done]
    if a.limit:
        todo = todo[:a.limit]
    print(f"{len(teams)} team-seasons, {len(done)} already scraped, {len(todo)} to go", flush=True)
    if not todo:
        return

    br = Browser()
    ok = fail = 0
    t0 = time.time()
    blocked_streak = 0
    with open(OUT, "a") as out:
        for i, t in enumerate(todo, 1):
            stid = t["season_team_id"]
            rec = {"season_team_id": stid, "team_id": t["team_id"], "school": t["team_name"],
                   "conference": t.get("conference", ""), "year": int(float(t["year"]))}
            for attempt in range(6):
                try:
                    page = br.get(f"https://stats.ncaa.org/teams/{stid}/roster")
                except Exception as e:
                    rec["error"] = f"{type(e).__name__}: {e}"[:150]
                    break
                if "Access Denied" in page:
                    # The ban is on the IP, not the session: a fresh browser profile is refused
                    # just the same. Rotating cookies therefore buys nothing and the only thing
                    # that clears it is waiting, so back off in minutes rather than seconds.
                    rec["error"] = "akamai denied"
                    blocked_streak += 1
                    back = min(900, 60 * (2 ** attempt)) + random.uniform(0, 20)
                    print(f"    blocked on {t['team_name']} {rec['year']}, "
                          f"waiting {back/60:.1f}m", flush=True)
                    time.sleep(back)
                    continue
                rows = parse(page)
                if rows is None:
                    rec["error"] = "no roster table"
                    break
                rec.pop("error", None)
                rec["players"] = [
                    {"name": r.get("Name", ""), "cls": r.get("Class", ""),
                     "pos": r.get("Position", ""), "height_raw": r.get("Height", ""),
                     "height_in": height_to_inches(r.get("Height", "")),
                     "bats": r.get("Bats", ""), "throws": r.get("Throws", ""),
                     "hometown": r.get("Hometown", "")}
                    for r in rows if r.get("Name")]
                blocked_streak = 0
                break

            out.write(json.dumps(rec) + "\n")
            out.flush()
            if "error" in rec:
                fail += 1
            else:
                ok += 1
            if i % 10 == 0 or i == len(todo):
                rate = i / max(time.time() - t0, 1)
                eta = (len(todo) - i) / max(rate, 1e-6) / 60
                print(f"  {i}/{len(todo)}  ok={ok} fail={fail}  "
                      f"{rate*60:.0f}/min  eta {eta:.0f}m", flush=True)
            if blocked_streak >= 3:
                print("  three consecutive blocks; pausing 20 minutes", flush=True)
                time.sleep(1200)
                blocked_streak = 0
            time.sleep(a.delay + random.uniform(0, 0.6))
    br.close()
    print(f"finished: ok={ok} fail={fail}  ->  {OUT}")


if __name__ == "__main__":
    main()

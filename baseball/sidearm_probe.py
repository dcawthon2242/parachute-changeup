#!/usr/bin/env python3
"""Check whether school athletics sites will serve baseball rosters with heights.

Before building anything around the Sidearm route it is worth confirming the two assumptions it
rests on: that the roster lives at a predictable path, and that plain HTTP is enough - no headless
browser, no bot challenge. If either fails the route is not actually cheaper than the NCAA one.
"""

import re, ssl, sys, urllib.error, urllib.request

UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/126.0 Safari/537.36")
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

PATHS = ["/sports/baseball/roster/2025", "/sports/baseball/roster", "/sports/bsb/2025-26/roster"]
SITES = {"Cal St. Fullerton": "fullertontitans.com", "Arizona": "arizonawildcats.com",
         "Vanderbilt": "vucommodores.com", "LSU": "lsusports.net",
         "Wake Forest": "godeacs.com", "East Carolina": "ecupirates.com"}


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "text/html"})
    with urllib.request.urlopen(req, timeout=25, context=CTX) as r:
        return r.getcode(), r.read().decode("utf-8", "ignore")


def looks_like_roster(h):
    heights = re.findall(r">\s*([5-7])['\u2019-](\d{1,2})\"?\s*<", h)
    heights += re.findall(r"\b([5-7])-(\d{1,2})\b", h)
    return len(heights), len(re.findall(r"(?i)sidearm", h))


for school, dom in SITES.items():
    hit = False
    for p in PATHS:
        url = f"https://{dom}{p}"
        try:
            code, body = get(url)
        except urllib.error.HTTPError as e:
            print(f"  {school:20} {p:32} HTTP {e.code}")
            continue
        except Exception as e:
            print(f"  {school:20} {p:32} {type(e).__name__}")
            continue
        n_h, n_sd = looks_like_roster(body)
        print(f"  {school:20} {p:32} {code}  {len(body):>7}B  heights~{n_h:<4} sidearm={n_sd>0}")
        if code == 200 and n_h >= 10:
            hit = True
            break
    if not hit:
        print(f"  {school:20} -> no roster path matched")

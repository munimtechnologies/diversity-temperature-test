# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Download helpers shared with extinction-temperature-test v0.1.0 (same PBDB
queries, cache file names and PhanDA file). In this repository only its
functions are used, by job_diversity.py.
Made with AI agents working in MV2's name.
"""
import datetime as dt
import hashlib
import os
import shutil
import time
import urllib.request

import thermo

PBDB = "https://paleobiodb.org/data1.2/occs/list.csv"
GROUPS = ["Porifera", "Cnidaria", "Bryozoa", "Brachiopoda", "Mollusca",
          "Echinodermata", "Trilobita", "Graptolithina"]
REQUIRED = {"Brachiopoda", "Mollusca", "Echinodermata", "Trilobita"}
COMMON = "envtype=marine&idreso=genus&idqual=genus_certain&pres=regular&show=class"
PHANDA_URL = ("https://raw.githubusercontent.com/EJJudd/PhanDA/main/5_Outputs/"
              "PhanDA_GMSTandCO2_percentiles.csv")
UA = "MV2-diversity-temperature/0.1 (Munim, Inc.; mv2@munimtech.com)"


def now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url, dest, timeout=1800, tries=3):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return "cached"
    tmp = f"{dest}.part{os.getpid()}"
    last = ""
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, "wb") as f:
                shutil.copyfileobj(r, f, 1 << 20)
            os.replace(tmp, dest)
            return "downloaded"
        except Exception as e:  # noqa: BLE001 - report and retry
            last = repr(e)
            print(f"fetch failed ({k + 1}/{tries}) {url}: {last}", flush=True)
            time.sleep(30 * (k + 1))
    raise RuntimeError(last)


def describe(path, url, how):
    return {"url": url, "status": how, "bytes": os.path.getsize(path), "md5": md5(path),
            "file_mtime_utc": dt.datetime.fromtimestamp(
                os.path.getmtime(path), dt.timezone.utc).isoformat()}


def get_pbdb(data_dir, prov):
    d = os.path.join(data_dir, "pbdb")
    os.makedirs(d, exist_ok=True)
    paths = []
    for g in GROUPS:
        url = f"{PBDB}?base_name={g}&{COMMON}"
        dest = os.path.join(d, f"{g}.csv")
        try:
            how = fetch(url, dest)
            with open(dest, encoding="utf-8", errors="replace") as f:
                head = f.readline()
            if "occurrence_no" not in head:
                os.remove(dest)
                raise RuntimeError("unexpected header: " + head[:300])
            prov[g] = describe(dest, url, how)
            paths.append((g, dest))
        except Exception as e:  # noqa: BLE001
            prov[g] = {"url": url, "error": repr(e)[:500]}
        print(g, prov[g].get("status", prov[g].get("error")), flush=True)
    return paths


def get_phanda(data_dir, prov):
    d = os.path.join(data_dir, "phanda")
    os.makedirs(d, exist_ok=True)
    dest = os.path.join(d, "PhanDA_GMSTandCO2_percentiles.csv")
    how = fetch(PHANDA_URL, dest)
    prov.update(describe(dest, PHANDA_URL, how))
    with open(dest, encoding="utf-8", errors="replace") as f:
        return thermo.parse_phanda(f.read())

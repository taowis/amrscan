#!/usr/bin/env python3
"""Fetch Pathogen Detection AST and isolate exports for archived PDTs in fixed batches.

The candidate-list query uses asm_acc, an isolate field, so it returns no rows
from the AST collection. AST is therefore requested by exact target_acc.
Writes <outdir>/batches/{ast,isolates}.NNNN.tsv, response headers and
batches.json, then recounts every batch to detect source changes.
"""
import argparse, csv, hashlib, json, re, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode, quote

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "bin"))
import ast_dataset as ad

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--candidates", type=Path, required=True, help="archived isolates TSV with target_acc first")
parser.add_argument("--outdir", type=Path, required=True)
parser.add_argument("--batch-size", type=int, default=100)
args = parser.parse_args()
candidates, out, size = args.candidates, args.outdir, args.batch_size

E = ad.ENDPOINT
FIELDS = {"ast": ad.AST_FIELDS, "isolates": ad.META_FIELDS + sorted(ad.META_OPTIONAL)}
ROLE = {"ast": "ast", "isolates": "metadata"}
now = lambda: datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
key = lambda p: tuple(int(x) for x in re.fullmatch(r"PDT(\d+)\.(\d+)", p).groups())

with candidates.open(encoding="utf-8-sig", newline="") as h:
    rows = [r for r in csv.reader(h, delimiter="\t") if r and r[0].strip()]
bad = sorted({r[0] for r in rows[1:] if not re.fullmatch(r"PDT\d+\.\d+", r[0])})
if bad:
    raise SystemExit(f"candidate list has {len(bad)} unversioned or malformed PDTs, e.g. {bad[:3]}")
pdts = sorted({r[0] for r in rows[1:]}, key=key)
batches = [pdts[i:i + size] for i in range(0, len(pdts), size)]
(out / "batches" / "headers").mkdir(parents=True, exist_ok=False)


def curl(url, body, headers):
    r = subprocess.run(["curl", "-sS", "--max-time", "300", "-D", str(headers), "-o", str(body),
                        "-w", "%{http_code}", url], capture_output=True, text=True)
    return r.returncode, r.stdout, r.stderr.strip()


def count(collection, flt, tag):
    url = E + "?" + urlencode({"action": "retrieve", "collection": collection, "filter": flt, "limit": 0, "fl": "id"})
    body = out / "batches" / "headers" / f"{tag}.count.json"
    for attempt in range(1, 4):
        t = now()
        code, status, err = curl(url, body, out / "batches" / "headers" / f"{tag}.count.headers.txt")
        try:
            d = json.loads(body.read_text())
            if code == 0 and status == "200" and d.get("success") is True:
                return {"url": url, "retrieved_at": t, "http_status": status, "attempts": attempt,
                        "total": d["ngout"]["data"]["totalCount"], "sha256": sha(body)}
        except (OSError, ValueError, KeyError):
            pass
        time.sleep(5 * attempt)
    raise SystemExit(f"count failed for {tag}: curl={code} http={status} {err}")


log = {"candidate_list": str(candidates), "candidate_list_sha256": sha(candidates),
       "distinct_pdts": len(pdts), "batch_size": size, "batch_count": len(batches),
       "order": "versioned PDT numeric ascending", "fields": FIELDS, "batches": []}
for n, batch in enumerate(batches, 1):
    flt = "target_acc:(" + " OR ".join(batch) + ")"
    entry = {"batch": n, "pdts": batch, "query": flt}
    for c in ("ast", "isolates"):
        tag = f"{c}.{n:04d}"
        expected = count(c, flt, tag)
        expr = f'[display()].from({c}).usingschema(/schema/pathogen).matching(q=="{quote(flt, safe="")}")'
        url = E + "?" + urlencode({"action": "solr2txt", "q": expr, "fields": ",".join(f"{f}|{f}" for f in FIELDS[c]),
                                   "filename": f"{c}.tsv", "type": "tsv", "nolimit": "on"})
        path = out / "batches" / f"{tag}.tsv"
        for attempt in range(1, 4):
            t = now()
            code, status, err = curl(url, path, out / "batches" / "headers" / f"{tag}.headers.txt")
            ok = code == 0 and status == "200"
            if ok:
                try:
                    records = ad.read_table(path, ROLE[c])
                except Exception as exc:
                    ok, err = False, f"parse: {exc}"
            if ok and len(records) == expected["total"]:
                break
            err = err or f"rows {len(records) if ok else 'NA'} != expected {expected['total']}"
            time.sleep(5 * attempt)
        else:
            raise SystemExit(f"export failed for {tag}: curl={code} http={status} {err}")
        targets = {r["target_acc"] for r in records}
        entry[c] = {"file": path.name, "url": url, "retrieved_at": t, "http_status": status, "attempts": attempt,
                    "rows": len(records), "expected_rows": expected["total"], "count_request": expected,
                    "bytes": path.stat().st_size, "sha256": sha(path),
                    "targets_outside_batch": sorted(targets - set(batch)),
                    "batch_pdts_without_rows": sorted(set(batch) - targets, key=key)}
    log["batches"].append(entry)
    print(n, entry["ast"]["rows"], entry["isolates"]["rows"], len(entry["ast"]["batch_pdts_without_rows"]), flush=True)
    (out / "batches.json").write_text(json.dumps(log, indent=1) + "\n")

# Second pass: re-count every batch to detect source changes during retrieval.
recheck_start = now()
for entry in log["batches"]:
    for c in ("ast", "isolates"):
        again = count(c, entry["query"], f"{c}.{entry['batch']:04d}.recheck")
        entry[c]["recheck_total"] = again["total"]
        entry[c]["recheck_at"] = again["retrieved_at"]
log["recheck_started_at"], log["finished_at"] = recheck_start, now()
(out / "batches.json").write_text(json.dumps(log, indent=1) + "\n")
print("done", log["finished_at"])

#!/usr/bin/env python3
"""Assemble a bin/ast_dataset.py snapshot from verified PDT batch exports.

Usage: python3 scripts/ecoli_research_v1/assemble_snapshot.py <census_dir>
Reads <census_dir>/batches.json and batches/, writes <census_dir>/snapshot/.
AST batch bodies are concatenated byte-for-byte under the shared header;
metadata batches are copied unchanged. batch_trace.json maps every combined
AST row range back to its batch file, URL, time and SHA-256.
"""
import hashlib, json, shutil, sys, tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "bin"))
census = Path(sys.argv[1])
import ast_dataset as ad

log = json.loads((census / "batches.json").read_text())
src, final = census / "batches", census / "snapshot"
if final.exists():
    raise SystemExit(f"refusing to replace {final}")
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

problems = []
for b in log["batches"]:
    for c in ("ast", "isolates"):
        e = b[c]
        if sha(src / e["file"]) != e["sha256"]:
            problems.append(f"hash changed: {e['file']}")
        if e["rows"] != e["expected_rows"] or e.get("recheck_total") != e["expected_rows"]:
            problems.append(f"count mismatch: {e['file']} rows={e['rows']} expected={e['expected_rows']} recheck={e.get('recheck_total')}")
        if e["targets_outside_batch"]:
            problems.append(f"rows outside batch: {e['file']}")
all_pdts = [p for b in log["batches"] for p in b["pdts"]]
if len(all_pdts) != len(set(all_pdts)) or len(all_pdts) != log["distinct_pdts"]:
    problems.append("batches overlap or do not cover the candidate list")
if problems:
    raise SystemExit("refusing to assemble:\n" + "\n".join(problems))

# Build in a staging directory and rename only after every check passes.
snap = Path(tempfile.mkdtemp(prefix=".snapshot-", dir=census))
try:

    header, trace, offset = None, [], 0
    with (snap / "ast.tsv").open("wb") as out:
        for b in log["batches"]:
            data = (src / b["ast"]["file"]).read_bytes()
            head, sep, body = data.partition(b"\n")
            if header is None:
                header = head
                out.write(head + b"\n")
            elif head != header:
                raise SystemExit(f"AST header differs in {b['ast']['file']}")
            if body and not body.endswith(b"\n"):
                body += b"\n"
            out.write(body)
            trace.append({"batch": b["batch"], "file": b["ast"]["file"], "sha256": b["ast"]["sha256"],
                          "url": b["ast"]["url"], "retrieved_at": b["ast"]["retrieved_at"],
                          "rows": b["ast"]["rows"], "combined_first_row": offset + 1,
                          "combined_last_row": offset + b["ast"]["rows"]})
            offset += b["ast"]["rows"]

    ast_rows = ad.read_table(snap / "ast.tsv", "ast")
    if len(ast_rows) != offset:
        raise SystemExit(f"combined AST rows {len(ast_rows)} != batch sum {offset}")
    first = min(b["ast"]["retrieved_at"] for b in log["batches"])
    last = max(b["ast"]["retrieved_at"] for b in log["batches"])
    sources = [{"role": "ast", "source_name": "NCBI Pathogen Detection",
                "source_url": f"{ad.ENDPOINT} solr2txt export of collection ast in {len(trace)} target_acc batches; see batch_trace.json",
                "retrieved_at": last, "filename": "ast.tsv",
                "original_filename": f"ast.0001.tsv..ast.{len(trace):04d}.tsv (concatenated)",
                "sha256": sha(snap / "ast.tsv"), "parser_version": ad.VERSION, "schema_version": ad.VERSION,
                "expected_rows": offset, "retrieval_window": [first, last]}]
    for b in log["batches"]:
        e = b["isolates"]
        shutil.copyfile(src / e["file"], snap / e["file"])
        sources.append({"role": "metadata", "source_name": "NCBI Pathogen Detection", "source_url": e["url"],
                        "retrieved_at": e["retrieved_at"], "filename": e["file"], "original_filename": "isolates.tsv",
                        "sha256": sha(snap / e["file"]), "parser_version": ad.VERSION, "schema_version": ad.VERSION,
                        "expected_rows": e["rows"]})
    (snap / "batch_trace.json").write_text(json.dumps({"batches_json_sha256": sha(census / "batches.json"),
                                                       "ast_header": header.decode(), "ast_batches": trace}, indent=1) + "\n")
    (snap / "manifest.json").write_text(json.dumps({
        "schema_version": ad.VERSION, "parser_version": ad.VERSION, "acquisition": "official-table-export",
        "query": f"target_acc batches of {log['batch_size']} from {log['candidate_list']} (sha256 {log['candidate_list_sha256']})",
        "sources": sources}, indent=2, sort_keys=True) + "\n")
    snap.rename(final)
except BaseException:
    shutil.rmtree(snap, ignore_errors=True)
    raise
print("assembled", offset, "AST rows,", len(sources) - 1, "metadata files")

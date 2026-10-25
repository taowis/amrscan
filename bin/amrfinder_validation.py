#!/usr/bin/env python3
"""Link accession-verified AMRFinderPlus evidence to source AST records."""

import argparse
import csv
from collections import Counter
import json
from pathlib import Path
import sys
import tempfile
from datetime import datetime, timezone
import hashlib

# Support both direct execution and test/import by file path.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import amrscan
import genome_dataset

AST_REQUIRED = {"record_id", "isolate_id", "biosample", "assembly_accession", "organism",
                "antimicrobial", "phenotype", "source_url", "source_sha256"}
EXCLUSION_FIELDS = ["pdt_accession", "biosample_accession", "assembly_accession", "reason"]
LINK_FIELDS = ["cohort_id", "isolate_id", "pdt_accession", "biosample_accession",
               "assembly_accession", "organism", "genome_file", "genome_sha256",
               "caller_version", "database_id", "genotype_status", "genotype_record_count",
               "ast_record_count"]


def read_tsv(path, required=None):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t", strict=True)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"{path}: missing or duplicate TSV header")
        if required and not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: missing columns {sorted(required - set(reader.fieldnames))}")
        rows = list(reader)
        if any(None in row or None in row.values() for row in rows):
            raise ValueError(f"{path}: malformed TSV row width")
    return rows


def validate_exclusions(path):
    if path is None:
        return []
    rows = read_tsv(path, set(EXCLUSION_FIELDS))
    for i, row in enumerate(rows, 2):
        if any(not (row[field] or "").strip() for field in EXCLUSION_FIELDS):
            raise ValueError(f"exclusions line {i}: all fields are required")
    return rows


def assemble(manifest, genomes, ast_path, calls_dir, outdir, exclusions_path=None):
    cohort = genome_dataset.read_manifest(manifest)
    genome_dataset.validate(manifest, genomes)
    ast_rows = read_tsv(ast_path, AST_REQUIRED)
    exclusions = validate_exclusions(exclusions_path)
    if {r["pdt_accession"] for r in cohort} & {r["pdt_accession"] for r in exclusions}:
        raise ValueError("an isolate cannot be both selected and excluded")

    by_pdt = {r["pdt_accession"]: r for r in cohort}
    ast_by_pdt = {pdt: [] for pdt in by_pdt}
    seen_record_ids = set()
    for row in ast_rows:
        pdt = row["isolate_id"]
        if pdt not in by_pdt:
            continue
        selected = by_pdt[pdt]
        if row["record_id"] in seen_record_ids:
            raise ValueError(f"duplicate AST record id: {row['record_id']}")
        seen_record_ids.add(row["record_id"])
        if row["biosample"] != selected["biosample_accession"]:
            raise ValueError(f"AST BioSample mismatch for {pdt}")
        if row["organism"] != selected["organism"]:
            raise ValueError(f"AST organism mismatch for {pdt}")
        if row["assembly_accession"] not in {"NA", selected["assembly_accession"]}:
            raise ValueError(f"AST Assembly mismatch for {pdt}")
        ast_by_pdt[pdt].append(row)
    missing_ast = [pdt for pdt, rows in ast_by_pdt.items() if not rows]
    if missing_ast:
        raise ValueError(f"selected isolates have no exact-PDT AST records: {missing_ast}")

    genomes = Path(genomes)
    calls_dir = Path(calls_dir)
    links, genotype_rows, linked_ast = [], [], []
    calls_by_pdt = {}
    for selected in cohort:
        pdt = selected["pdt_accession"]
        genome_entry = next(x for x in json.loads((genomes / "provenance.json").read_text())["assemblies"]
                            if x["assembly_accession"] == selected["assembly_accession"])
        stem = calls_dir / pdt
        raw_path = Path(str(stem) + ".amrfinder.tsv")
        provenance_path = Path(str(stem) + ".provenance.json")
        harmonized_path = Path(str(stem) + ".harmonized.tsv")
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        if (provenance.get("sample_id") != pdt or provenance.get("tool") != "AMRFinderPlus"
                or provenance.get("status") != "success"):
            raise ValueError(f"missing successful real caller provenance for {pdt}")
        if provenance.get("input_sha256") != genome_entry["sha256"]:
            raise ValueError(f"caller input hash does not match verified assembly for {pdt}")
        if provenance.get("caller_version") in {None, "", "NA"}:
            raise ValueError(f"caller version is missing for {pdt}")
        database_files = provenance.get("database_files_sha256")
        if not isinstance(database_files, dict) or not database_files:
            raise ValueError(f"database fingerprint is missing for {pdt}")
        database_id = "sha256:" + hashlib.sha256(amrscan.json_text(database_files).encode()).hexdigest()
        if provenance.get("database_id") != database_id:
            raise ValueError(f"database fingerprint mismatch for {pdt}")
        if provenance.get("raw_sha256") != amrscan.sha256(raw_path):
            raise ValueError(f"caller report hash mismatch for {pdt}")
        records = read_tsv(harmonized_path, set(amrscan.SCHEMA))
        if any(row["sample_id"] != pdt or row["tool"] != "AMRFinderPlus" for row in records):
            raise ValueError(f"harmonized sample/tool mismatch for {pdt}")
        if any(row["database_id"] != provenance.get("database_id") or
               row["caller_version"] != provenance.get("caller_version") for row in records):
            raise ValueError(f"harmonized provenance mismatch for {pdt}")
        with tempfile.TemporaryDirectory(prefix="amrfinder-harmonize-check-") as tmp:
            rebuilt = Path(tmp) / "expected.tsv"
            amrscan.harmonize(raw_path, provenance_path, pdt, rebuilt)
            if amrscan.sha256(rebuilt) != amrscan.sha256(harmonized_path):
                raise ValueError(f"harmonized report does not match raw caller output for {pdt}")
        calls_by_pdt[pdt] = records
        status = "determinants_reported" if any(r["element_type"] == "AMR" for r in records) else "no_determinants_reported"
        links.append({"cohort_id": selected["cohort_id"], "isolate_id": selected["isolate_id"],
                      "pdt_accession": pdt, "biosample_accession": selected["biosample_accession"],
                      "assembly_accession": selected["assembly_accession"], "organism": selected["organism"],
                      "genome_file": genome_entry["file"], "genome_sha256": genome_entry["sha256"],
                      "caller_version": provenance["caller_version"], "database_id": provenance["database_id"],
                      "genotype_status": status, "genotype_record_count": len(records),
                      "ast_record_count": len(ast_by_pdt[pdt])})
        for row in records:
            genotype_rows.append({"cohort_id": selected["cohort_id"], "isolate_id": selected["isolate_id"],
                                  "pdt_accession": pdt, "biosample_accession": selected["biosample_accession"],
                                  "assembly_accession": selected["assembly_accession"], **row})
        for row in ast_by_pdt[pdt]:
            linked_ast.append({"cohort_id": selected["cohort_id"], "pdt_accession": pdt,
                               "verified_assembly_accession": selected["assembly_accession"], **row})

    elements = [row for rows in calls_by_pdt.values() for row in rows if row["element_type"] == "AMR"]
    ast_counts = Counter(row["phenotype"] for rows in ast_by_pdt.values() for row in rows)
    summary = {
        "selected_cohort_isolates": len(cohort),
        "valid_assemblies": len(links),
        "successfully_genotyped_isolates": len(links),
        "isolates_with_determinant_calls": sum(link["genotype_status"] == "determinants_reported" for link in links),
        "isolates_with_no_determinant_calls": sum(link["genotype_status"] == "no_determinants_reported" for link in links),
        "isolates_linked_to_ast": sum(link["ast_record_count"] > 0 for link in links),
        "excluded_isolates": len(exclusions),
        "exclusions_by_reason": dict(sorted(Counter(row["reason"] for row in exclusions).items())),
        "amr_determinant_records": len(elements),
        "gene_family_counts": dict(sorted(Counter(row["gene"] for row in elements).items())),
        "element_type_counts": dict(sorted(Counter(row["element_type"] for rows in calls_by_pdt.values()
                                                    for row in rows).items())),
        "ast_observation_counts_by_source_phenotype": dict(sorted(ast_counts.items())),
        "interpretation": "Separate source AST and caller evidence; no phenotype prediction or concordance analysis.",
    }
    outdir = Path(outdir)
    if outdir.exists() or outdir.is_symlink():
        raise ValueError(f"output already exists: {outdir}")
    outdir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".amrfinder-validation-", dir=outdir.parent))
    try:
        write_tsv(staging / "isolate_links.tsv", LINK_FIELDS, links)
        write_tsv(staging / "genotype_evidence.tsv", ["cohort_id", "isolate_id", "pdt_accession",
                  "biosample_accession", "assembly_accession", *amrscan.SCHEMA], genotype_rows)
        write_tsv(staging / "ast_evidence.tsv", ["cohort_id", "pdt_accession", "verified_assembly_accession",
                  *ast_rows[0].keys()], linked_ast)
        write_tsv(staging / "excluded_isolates.tsv", EXCLUSION_FIELDS, exclusions)
        (staging / "validation_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n",
                                                        encoding="utf-8")
        inputs = {"manifest_sha256": genome_dataset.sha256(manifest),
                  "genome_provenance_sha256": genome_dataset.sha256(genomes / "provenance.json"),
                  "ast_table_sha256": genome_dataset.sha256(ast_path),
                  "caller_reports": {pdt: amrscan.sha256(calls_dir / f"{pdt}.amrfinder.tsv") for pdt in by_pdt}}
        run = {"schema_version": "1.0", "created_at": datetime.now(timezone.utc).isoformat(),
               "inputs": inputs, "summary_sha256": hashlib.sha256(
                   json.dumps(summary, sort_keys=True, separators=(",", ":")).encode()).hexdigest()}
        (staging / "validation_provenance.json").write_text(json.dumps(run, indent=2, sort_keys=True) + "\n",
                                                             encoding="utf-8")
        if outdir.exists() or outdir.is_symlink():
            raise ValueError(f"output appeared during build: {outdir}")
        staging.rename(outdir)
    finally:
        if staging.exists():
            import shutil
            shutil.rmtree(staging)
    return summary


def write_tsv(path, fields, rows):
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    cmd = parser.add_subparsers(dest="action", required=True).add_parser("summarize")
    for name in ("manifest", "genomes", "ast", "calls", "outdir"):
        cmd.add_argument("--" + name, type=Path, required=True)
    cmd.add_argument("--exclusions", type=Path)
    args = parser.parse_args(argv)
    try:
        result = assemble(args.manifest, args.genomes, args.ast, args.calls, args.outdir, args.exclusions)
        print(json.dumps(result, indent=2, sort_keys=True))
    except (ValueError, KeyError, OSError, csv.Error, json.JSONDecodeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

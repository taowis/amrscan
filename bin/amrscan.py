#!/usr/bin/env python3
"""Run AMRFinderPlus and retain its evidence in the AMRScan v1 TSV schema."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess


SCHEMA = [
    "sample_id", "tool", "gene", "allele", "resistance_class", "drug",
    "identity", "coverage", "contig", "start", "end", "strand", "evidence_type",
    "schema_version", "caller_version", "database_id", "thresholds_json",
    "element_type", "element_subtype", "scope", "subclass",
    "raw_file", "raw_row", "raw_record_json",
]
REQUIRED = {
    "Contig id", "Start", "Stop", "Strand", "Element symbol", "Class",
    "Subclass", "Method", "Type", "Subtype", "Scope",
    "% Identity to reference", "% Coverage of reference",
}
MISSING = {"", "NA", "N/A", "."}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_text(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def sample_id(value):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError("Sample ID must start with a letter/digit and contain only letters, digits, _, . or -")
    return value


def normalized(value):
    return "NA" if value in MISSING else value


def numeric(value, field, coordinate=False):
    if value in MISSING:
        return "NA"
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Non-finite {field}: {value}")
    if coordinate:
        if not number.is_integer() or number < 1:
            raise ValueError(f"Invalid 1-based {field}: {value}")
    elif not 0 <= number <= 100:
        raise ValueError(f"{field} must be a percentage in [0,100]: {value}")
    return value


def harmonize(raw, provenance, sample, output):
    sample_id(sample)
    raw, output = Path(raw), Path(output)
    if output.resolve() in {raw.resolve(), Path(provenance).resolve()}:
        raise ValueError("Output must not overwrite raw evidence or provenance")
    metadata = json.loads(Path(provenance).read_text())
    if metadata["sample_id"] != sample or metadata["tool"] != "AMRFinderPlus":
        raise ValueError("Provenance sample/tool does not match this input")
    if metadata["status"] not in {"success", "upstream_fixture"}:
        raise ValueError("Cannot harmonize a failed caller run")
    if metadata["raw_sha256"] != sha256(raw):
        raise ValueError("Raw report checksum does not match provenance")
    # Buffer a sample so a malformed row cannot publish a partially valid report.
    rows = []
    with raw.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        headers = reader.fieldnames or []
        if len(headers) != len(set(headers)) or not REQUIRED.issubset(headers):
            raise ValueError("Missing or duplicate AMRFinderPlus columns; expected the v4 nucleotide report")
        for row_number, record in enumerate(reader, start=1):
            if None in record or None in record.values():
                raise ValueError(f"Malformed TSV record {row_number}")
            method = normalized(record["Method"])
            if method == "NA":
                raise ValueError(f"Missing evidence Method in record {row_number}")
            start = numeric(record["Start"], "start", coordinate=True)
            end = numeric(record["Stop"], "end", coordinate=True)
            if (start == "NA") != (end == "NA"):
                raise ValueError("Coordinates must both be present or both be NA")
            if start != "NA" and int(start) > int(end):
                raise ValueError("Start exceeds Stop; strand is represented separately")
            strand = normalized(record["Strand"])
            if strand not in {"+", "-", "NA"}:
                raise ValueError(f"Invalid strand: {strand}")
            element = normalized(record["Element symbol"])
            rows.append(dict(zip(SCHEMA, [
                sample, "AMRFinderPlus", element,
                element if method in {"ALLELE", "ALLELEX", "ALLELEP"} else "NA",
                normalized(record["Class"]) if record["Type"] == "AMR" else "NA",
                "NA",  # Subclass can be a class or several drugs, not an individual AST endpoint.
                numeric(record["% Identity to reference"], "identity"),
                numeric(record["% Coverage of reference"], "coverage"),
                normalized(record["Contig id"]), start, end, strand, method,
                "1.0", metadata["caller_version"], metadata["database_id"],
                json_text(metadata["thresholds"]), normalized(record["Type"]),
                normalized(record["Subtype"]), normalized(record["Scope"]),
                normalized(record["Subclass"]), raw.name, row_number, json_text(record),
            ])))
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SCHEMA, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def call(args):
    sample_id(args.sample)
    if args.threads < 1:
        raise ValueError("threads must be positive")
    if not (args.identity == -1 or 0 <= args.identity <= 1) or not 0 <= args.coverage <= 1:
        raise ValueError("identity must be -1 or [0,1]; coverage must be [0,1]")
    database = Path(args.database)
    if not database.is_dir():
        raise ValueError("AMRFinderPlus database must be an existing directory")
    db_files = {p.relative_to(database).as_posix(): sha256(p)
                for p in sorted(database.rglob("*")) if p.is_file()}
    if not db_files:
        raise ValueError("AMRFinderPlus database directory is empty")
    executable = shutil.which(args.executable)
    if not executable:
        raise ValueError(f"AMRFinderPlus executable not found: {args.executable}")
    version = subprocess.run([executable, "--version"], check=True, text=True,
                             capture_output=True).stdout.strip()
    if not version:
        raise ValueError("AMRFinderPlus did not report its software version")
    command = [executable, "--nucleotide", str(args.input), "--database", str(database),
               "--threads", str(args.threads), "--ident_min", str(args.identity),
               "--coverage_min", str(args.coverage), "--translation_table", "11",
               "--report_all_equal", "--print_node", "-o", str(args.raw)]
    if args.organism:
        command.extend(["--organism", args.organism])
    if args.plus:
        command.append("--plus")
    metadata = {
        "schema_version": "1.0", "sample_id": args.sample, "tool": "AMRFinderPlus",
        "caller_version": version, "command": command, "input_sha256": sha256(args.input),
        "database_source": args.database_source or str(database.resolve()),
        "database_files_sha256": db_files,
        "database_id": "sha256:" + hashlib.sha256(json_text(db_files).encode()).hexdigest(),
        "thresholds": {"ident_min": args.identity, "coverage_min": args.coverage,
                       "translation_table": 11, "report_all_equal": True,
                       "organism": args.organism, "plus": args.plus},
        "status": "failed",
    }
    try:
        with Path(args.log).open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        metadata["raw_sha256"] = sha256(args.raw)
        metadata["status"] = "success"
    finally:
        Path(args.provenance).write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    run = sub.add_parser("call", help="Run a real, locally installed AMRFinderPlus")
    run.add_argument("--input", required=True, type=Path)
    run.add_argument("--database", required=True, type=Path)
    run.add_argument("--database-source")
    run.add_argument("--executable", default="amrfinder")
    run.add_argument("--threads", type=int, default=4)
    run.add_argument("--identity", type=float, default=-1)
    run.add_argument("--coverage", type=float, default=0.5)
    run.add_argument("--organism")
    run.add_argument("--plus", action="store_true")
    run.add_argument("--log", required=True, type=Path)
    convert = sub.add_parser("harmonize", help="Convert a report without discarding caller evidence")
    convert.add_argument("--output", required=True, type=Path)
    for child in (run, convert):
        child.add_argument("--sample", required=True)
        child.add_argument("--raw", required=True, type=Path)
        child.add_argument("--provenance", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        if args.action == "call":
            call(args)
        else:
            harmonize(args.raw, args.provenance, args.sample, args.output)
    except (ValueError, KeyError, OSError, csv.Error, subprocess.CalledProcessError) as error:
        parser.exit(1, f"AMRScan: {error}\n")


if __name__ == "__main__":
    main()

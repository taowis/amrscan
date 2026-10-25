#!/usr/bin/env python3
"""Retrieve and audit explicitly accessioned NCBI genome assemblies."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
import sys
import tempfile
from urllib.error import URLError
from urllib.request import Request, urlopen
import zipfile

FIELDS = ["cohort_id", "isolate_id", "pdt_accession", "biosample_accession",
          "assembly_accession", "organism", "assembly_source", "ast_available",
          "selection_reason"]
ACCESSION = re.compile(r"^(?:GCF|GCA)_\d{9}\.\d+$")
BIOSAMPLE = re.compile(r"^SAM[END]\d+$")
PDT = re.compile(r"^PDT\d+\.\d+$")
MAX_ARCHIVE = 100 * 1024 * 1024
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t", strict=True)
        if reader.fieldnames != FIELDS:
            raise ValueError(f"manifest columns must be exactly: {' '.join(FIELDS)}")
        rows, seen_assemblies, seen_isolates = [], set(), set()
        for line, row in enumerate(reader, 2):
            if None in row or any(not (v or "").strip() for v in row.values()):
                raise ValueError(f"manifest line {line}: missing or extra field")
            row = {key: value.strip() for key, value in row.items()}
            if not ACCESSION.fullmatch(row["assembly_accession"]):
                raise ValueError(f"manifest line {line}: invalid Assembly accession")
            if not BIOSAMPLE.fullmatch(row["biosample_accession"]):
                raise ValueError(f"manifest line {line}: invalid BioSample accession")
            if not PDT.fullmatch(row["pdt_accession"]):
                raise ValueError(f"manifest line {line}: invalid PDT accession")
            if row["ast_available"] != "true":
                raise ValueError(f"manifest line {line}: cohort requires a linked AST observation")
            if row["organism"] != "Escherichia coli":
                raise ValueError(f"manifest line {line}: organism must be exact Escherichia coli")
            if row["assembly_accession"] in seen_assemblies:
                raise ValueError(f"duplicate assembly accession: {row['assembly_accession']}")
            if row["isolate_id"] in seen_isolates:
                raise ValueError(f"duplicate isolate mapping: {row['isolate_id']}")
            seen_assemblies.add(row["assembly_accession"])
            seen_isolates.add(row["isolate_id"])
            rows.append(row)
    if not rows:
        raise ValueError("manifest contains no validation isolates")
    return rows


def request_bytes(url, limit):
    if not url.startswith("https://"):
        raise ValueError("only HTTPS sources are accepted")
    req = Request(url, headers={"User-Agent": "AMRScan-genome-dataset/1.0"})
    with urlopen(req, timeout=120) as response:
        blocks, size = [], 0
        while True:
            block = response.read(1024 * 1024)
            if not block:
                break
            size += len(block)
            if size > limit:
                raise ValueError(f"download exceeds {limit} bytes")
            blocks.append(block)
        return b"".join(blocks), response.geturl()


def parse_fasta(data):
    text = data.decode("ascii")
    records, bases, header = 0, 0, False
    for line_number, line in enumerate(text.splitlines(), 1):
        if line.startswith(">"):
            if not line[1:].strip():
                raise ValueError(f"empty FASTA header at line {line_number}")
            records += 1
            header = True
        elif line.strip():
            if not header or not re.fullmatch(r"[A-Za-z*.-]+", line.strip()):
                raise ValueError(f"invalid FASTA sequence at line {line_number}")
            bases += len(line.strip())
    if records == 0 or bases == 0:
        raise ValueError("FASTA must contain a header and sequence")
    return records, bases


def package_fasta(archive, assembly):
    source = io.BytesIO(archive) if isinstance(archive, bytes) else archive
    with zipfile.ZipFile(source) as package:
        names = [n for n in package.namelist()
                 if n.startswith(f"ncbi_dataset/data/{assembly}/") and n.endswith("_genomic.fna")]
        if len(names) != 1:
            raise ValueError(f"NCBI package must contain exactly one FASTA for {assembly}; found {len(names)}")
        return package.read(names[0]), names[0]


def package_identity(archive, row):
    """Read an NCBI package and verify identity across current and legacy reports."""
    source = io.BytesIO(archive) if isinstance(archive, bytes) else archive
    with zipfile.ZipFile(source) as package:
        names_in_package = set(package.namelist())
        if "md5sum.txt" not in names_in_package:
            raise ValueError("NCBI package lacks md5sum.txt")
        checksum_rows = package.read("md5sum.txt").decode("ascii").splitlines()
        checked = set()
        for line in checksum_rows:
            match = re.fullmatch(r"([0-9a-fA-F]{32})[ \t]+\*?([^\t\r\n]+)", line)
            if not match:
                raise ValueError("invalid NCBI md5sum.txt row")
            expected_md5, name = match.groups()
            if name.startswith("/") or ".." in Path(name).parts or name not in names_in_package or name in checked:
                raise ValueError(f"invalid or missing member in NCBI md5sum.txt: {name}")
            if hashlib.md5(package.read(name)).hexdigest().lower() != expected_md5.lower():
                raise ValueError(f"NCBI package MD5 mismatch: {name}")
            checked.add(name)
        if not checked:
            raise ValueError("NCBI md5sum.txt is empty")
        names = [n for n in package.namelist()
                 if n.startswith(f"ncbi_dataset/data/{row['assembly_accession']}/")
                 and n.endswith("_genomic.fna")]
        if len(names) != 1:
            raise ValueError(f"NCBI package must contain exactly one FASTA for {row['assembly_accession']}; found {len(names)}")
        member = names[0]
        fasta = package.read(member)
        reports = [n for n in package.namelist() if n.endswith("assembly_data_report.jsonl")]
        if len(reports) != 1:
            raise ValueError("NCBI package lacks a unique assembly data report")
        report_lines = package.read(reports[0]).decode("utf-8").splitlines()
    report = next((json.loads(line) for line in report_lines if line.strip()), {})
    assembly_info = report.get("assemblyInfo", {})
    biosample_info = assembly_info.get("biosample", {})
    description = biosample_info.get("description", {})
    reported_accession = assembly_info.get("assemblyAccession") or report.get("currentAccession")
    biosample = (report.get("biologicalMaterial", {}).get("biosample", {}).get("accession")
                 or biosample_info.get("accession"))
    organism = (report.get("organism", {}).get("organismName")
                or description.get("organism", {}).get("organismName"))
    if reported_accession != row["assembly_accession"] or biosample != row["biosample_accession"]:
        raise ValueError(f"NCBI identity mismatch for {row['assembly_accession']}: accession={reported_accession}, BioSample={biosample}")
    if organism != row["organism"]:
        raise ValueError(f"NCBI organism mismatch for {row['assembly_accession']}: {organism}")
    if member not in checked or reports[0] not in checked:
        raise ValueError("NCBI md5sum.txt does not cover the Assembly report and genome FASTA")
    records, bases = parse_fasta(fasta)
    return fasta, member, records, bases


def timestamp(value):
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        raise ValueError("retrieval timestamp must include timezone")
    return value


def download(manifest, outdir):
    rows = read_manifest(manifest)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    provenance = outdir / "provenance.json"
    if provenance.exists():
        raise ValueError(f"refusing to replace existing provenance: {provenance}")
    entries = []
    for row in rows:
        assembly = row["assembly_accession"]
        fasta_path = outdir / f"{row['isolate_id']}_{assembly}.fna"
        if fasta_path.exists():
            raise ValueError(f"refusing to replace existing assembly file: {fasta_path}")
        url = API.format(assembly) + "/download?include_annotation_type=GENOME_FASTA,SEQUENCE_REPORT&filename=assembly.zip"
        try:
            archive, final_url = request_bytes(url, MAX_ARCHIVE)
            with tempfile.TemporaryDirectory(dir=outdir) as tmp:
                zip_path = Path(tmp) / "assembly.zip"
                zip_path.write_bytes(archive)
                fasta, member, records, bases = package_identity(zip_path, row)
                # Exclusive create plus atomic rename; never overwrite another file.
                with tempfile.NamedTemporaryFile(dir=outdir, prefix=".assembly-", delete=False) as tmp_fasta:
                    temp_path = Path(tmp_fasta.name)
                    tmp_fasta.write(fasta)
                if fasta_path.exists():
                    raise ValueError(f"assembly file appeared during download: {fasta_path}")
                temp_path.replace(fasta_path)
            entries.append({**row, "file": fasta_path.name, "sha256": sha256(fasta_path),
                            "package_sha256": hashlib.sha256(archive).hexdigest(),
                            "source_url": final_url, "retrieved_at": datetime.now(timezone.utc).isoformat(),
                            "package_member": member, "sequence_records": records, "sequence_bases": bases})
        except Exception:
            if fasta_path.exists():
                fasta_path.unlink()
            raise
    payload = {"schema_version": "1.0", "manifest_sha256": sha256(manifest), "assemblies": entries}
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=outdir, prefix=".provenance-", delete=False) as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        temp_manifest = Path(stream.name)
    if provenance.exists():
        temp_manifest.unlink()
        raise ValueError(f"refusing to replace existing provenance: {provenance}")
    temp_manifest.replace(provenance)


def read_package_manifest(path, rows):
    fields = ["assembly_accession", "package_path", "source_url", "retrieved_at"]
    result, seen = {}, set()
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t", strict=True)
        if reader.fieldnames != fields:
            raise ValueError(f"package manifest columns must be exactly: {' '.join(fields)}")
        for line, item in enumerate(reader, 2):
            if None in item or any(not (value or "").strip() for value in item.values()):
                raise ValueError(f"package manifest line {line}: missing or extra field")
            item = {key: value.strip() for key, value in item.items()}
            assembly = item["assembly_accession"]
            if assembly not in {row["assembly_accession"] for row in rows}:
                raise ValueError(f"package manifest assembly does not match cohort: {assembly}")
            if assembly in seen:
                raise ValueError(f"duplicate package manifest assembly: {assembly}")
            if not item["source_url"].startswith("https://"):
                raise ValueError("only HTTPS sources are accepted")
            timestamp(item["retrieved_at"])
            package_path = Path(item["package_path"])
            if not package_path.is_absolute():
                package_path = Path(path).resolve().parent / package_path
            if not package_path.is_file() or package_path.is_symlink():
                raise ValueError(f"package is missing or is a symlink: {package_path}")
            if package_path.stat().st_size > MAX_ARCHIVE:
                raise ValueError(f"package exceeds {MAX_ARCHIVE} bytes: {package_path}")
            item["package_path"] = package_path
            result[assembly] = item
            seen.add(assembly)
    expected = {row["assembly_accession"] for row in rows}
    if seen != expected:
        raise ValueError(f"package manifest assemblies do not match cohort; missing={sorted(expected-seen)}, extra={sorted(seen-expected)}")
    return result


def import_packages(manifest, packages_tsv, outdir):
    """Import manually retrieved NCBI packages while preserving retrieval provenance."""
    rows = read_manifest(manifest)
    packages = read_package_manifest(packages_tsv, rows)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    provenance = outdir / "provenance.json"
    if provenance.exists():
        raise ValueError(f"refusing to replace existing provenance: {provenance}")
    destinations = [outdir / f"{row['isolate_id']}_{row['assembly_accession']}.fna" for row in rows]
    if any(path.exists() for path in destinations):
        raise ValueError("refusing to replace existing assembly file")
    entries = []
    with tempfile.TemporaryDirectory(dir=outdir, prefix=".import-") as staging_dir:
        staging = Path(staging_dir)
        for row, destination in zip(rows, destinations):
            source = packages[row["assembly_accession"]]
            package_path = source["package_path"]
            fasta, member, records, bases = package_identity(package_path, row)
            staged = staging / destination.name
            staged.write_bytes(fasta)
            entries.append({**row, "file": destination.name, "sha256": sha256(staged),
                            "source_url": source["source_url"],
                            "retrieved_at": source["retrieved_at"],
                            "package_member": member, "package_sha256": sha256(package_path),
                            "sequence_records": records, "sequence_bases": bases})
        payload = {"schema_version": "1.0", "manifest_sha256": sha256(manifest),
                   "assemblies": entries}
        staged_provenance = staging / "provenance.json"
        staged_provenance.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                                     encoding="utf-8")
        moved = []
        try:
            for destination in destinations:
                if destination.exists():
                    raise ValueError(f"assembly file appeared during import: {destination}")
                (staging / destination.name).replace(destination)
                moved.append(destination)
            if provenance.exists():
                raise ValueError(f"provenance appeared during import: {provenance}")
            staged_provenance.replace(provenance)
        except Exception:
            for path in moved:
                path.unlink(missing_ok=True)
            raise
    return len(rows)


def validate(manifest, genomes):
    rows = read_manifest(manifest)
    genomes = Path(genomes)
    payload = json.loads((genomes / "provenance.json").read_text(encoding="utf-8"))
    if payload.get("manifest_sha256") != sha256(manifest):
        raise ValueError("cohort manifest hash differs from retrieval provenance")
    by_assembly = {item["assembly_accession"]: item for item in payload["assemblies"]}
    if len(by_assembly) != len(payload["assemblies"]):
        raise ValueError("duplicate assembly entries in provenance")
    if set(by_assembly) != {row["assembly_accession"] for row in rows}:
        raise ValueError("manifest and provenance assembly sets differ")
    for row in rows:
        item = by_assembly[row["assembly_accession"]]
        for field in ("biosample_accession", "organism", "isolate_id"):
            if item[field] != row[field]:
                raise ValueError(f"{field} mismatch for {row['assembly_accession']}")
        path = genomes / item["file"]
        if path.parent != genomes or not path.is_file() or sha256(path) != item["sha256"]:
            raise ValueError(f"missing file or SHA-256 mismatch for {row['assembly_accession']}")
        parse_fasta(path.read_bytes())
    return len(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    cmd = commands.add_parser("download")
    cmd.add_argument("--manifest", required=True, type=Path)
    cmd.add_argument("--outdir", required=True, type=Path)
    cmd = commands.add_parser("validate")
    cmd.add_argument("--manifest", required=True, type=Path)
    cmd.add_argument("--genomes", required=True, type=Path)
    cmd = commands.add_parser("import-packages")
    cmd.add_argument("--manifest", required=True, type=Path)
    cmd.add_argument("--packages-tsv", required=True, type=Path)
    cmd.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        if args.action == "download":
            download(args.manifest, args.outdir)
            print("download complete")
        elif args.action == "validate":
            print(f"validated assemblies: {validate(args.manifest, args.genomes)}")
        else:
            import_packages(args.manifest, args.packages_tsv, args.outdir)
            print("package imported and verified")
    except (ValueError, OSError, KeyError, json.JSONDecodeError, UnicodeError, zipfile.BadZipFile, URLError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Archive and normalize NCBI AST evidence without interpreting breakpoints."""

import argparse
from collections import Counter, defaultdict
from contextlib import contextmanager
import csv
from datetime import date, datetime, timezone
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile
from typing import Iterator
from urllib.parse import quote, urlencode
from urllib.request import urlopen

VERSION = "1.0.0"
ENDPOINT = "https://www.ncbi.nlm.nih.gov/pathogens/pathogens-srv/"
DEFAULT_QUERY = 'scientific_name:"Escherichia coli"'
POLICY = Path(__file__).resolve().parents[1] / "ingestion/ncbi_pathogen_detection/policy.json"
NA = "NA"
MISSING = {"", "na", "n/a", "null", "missing", "not collected", "not provided", "not applicable"}
AST_FIELDS = (
    "id target_acc biosample_acc scientific_name antibiotic phenotype measurement_sign "
    "mic mic_secondary disk_diffusion standard platform vendor reagent host collection_date "
    "geo_loc_name isolation_source epi_type bioproject_acc taxgroup_name"
).split()
META_FIELDS = (
    "target_acc biosample_acc asm_acc Run scientific_name collection_date geo_loc_name "
    "isolation_source host epi_type bioproject_acc species_taxid isolate_identifiers"
).split()
# Optional metadata enrichments used by downstream cohort summaries. Keeping
# them optional preserves compatibility with previously archived snapshots.
META_OPTIONAL = {"strain", "taxid", "erd_group", "bioproject_center"}
# These two optional columns support archived tables enriched by an identified source.
AST_OPTIONAL = {"method", "standard_version"}
SCHEMA = (
    "record_id isolate_id biosample assembly_accession sra_accession organism antimicrobial "
    "antimicrobial_original phenotype phenotype_original mic mic_operator mic_units mic_raw "
    "mic_secondary mic_secondary_raw measurement_sign_original disk_diffusion_raw "
    "testing_method breakpoint_standard breakpoint_version platform vendor reagent "
    "collection_date collection_year country location isolation_source host source_category "
    "source_file source_row source_url source_sha256 retrieved_at source_record_id "
    "metadata_records_json raw_record_json schema_version policy_version qc_flags"
).split()
FLAGS = set((
    "missing_isolate_id missing_genome_link missing_assembly missing_phenotype unknown_phenotype "
    "unknown_antimicrobial invalid_mic conflicting_ast missing_country missing_collection_year "
    "missing_source organism_mismatch duplicate_record multiple_ast_records missing_metadata "
    "metadata_conflict missing_testing_method missing_breakpoint_standard missing_breakpoint_version "
    "invalid_collection_date ambiguous_source ambiguous_country"
).split())
ENDPOINT_SCHEMA = ("isolate_id antimicrobial phenotype eligible exclusion_reasons record_ids_json "
                   "record_count distinct_observations").split()


def canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def digest(path: Path) -> str:
    with path.open("rb") as handle:
        result = hashlib.sha256()
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def clean(value: str) -> str:
    return NA if value.strip().casefold() in MISSING else value.strip()


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


@contextmanager
def new_directory(destination: Path) -> Iterator[Path]:
    """Publish only completed artifacts; refuse to replace any existing path."""
    if destination.exists() or destination.is_symlink():
        raise ValueError(f"output already exists: {destination}; choose a new directory")
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".ast-", dir=destination.parent))
    try:
        yield staging
        if destination.exists() or destination.is_symlink():
            raise ValueError(f"output appeared during build: {destination}")
        staging.rename(destination)
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def read_table(path: Path, role: str) -> list[dict[str, str]]:
    fields = AST_FIELDS if role == "ast" else META_FIELDS
    allowed = set(fields) | (AST_OPTIONAL if role == "ast" else META_OPTIONAL)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", strict=True)
        header = next(reader, [])
        if header:
            header[0] = header[0].removeprefix("#")
        if len(set(header)) != len(header):
            raise ValueError(f"{path.name}: duplicate columns")
        optional = AST_OPTIONAL if role == "ast" else META_OPTIONAL
        missing, extra = set(fields) - set(header), set(header) - allowed
        missing -= optional
        if missing or extra:
            raise ValueError(f"{path.name}: schema mismatch; missing={sorted(missing)}, unexpected={sorted(extra)}")
        result = []
        for number, values in enumerate(reader, 1):
            if len(values) != len(header):
                raise ValueError(f"{path.name}: data row {number}: expected {len(header)} cells, got {len(values)}")
            result.append(dict(zip(header, values)))
    return result


def export_url(collection: str, query: str, fields: list[str]) -> str:
    expression = f'[display()].from({collection}).usingschema(/schema/pathogen).matching(q=="{quote(query, safe="")}")'
    return ENDPOINT + "?" + urlencode({
        "action": "solr2txt", "q": expression,
        "fields": ",".join(f"{field}|{field}" for field in fields),
        "filename": f"{collection}.tsv", "type": "tsv", "nolimit": "on",
    })


def retrieve(url: str, output: Path, max_bytes: int) -> None:
    with urlopen(url, timeout=60) as response, output.open("xb") as handle:
        size = 0
        while block := response.read(1024 * 1024):
            size += len(block)
            if size > max_bytes:
                raise ValueError(f"download exceeds {max_bytes} bytes; narrow --query")
            handle.write(block)


def remote_count(collection: str, query: str) -> int:
    url = ENDPOINT + "?" + urlencode({"action": "retrieve", "collection": collection,
                                     "filter": query, "limit": 0, "fl": "id"})
    with urlopen(url, timeout=60) as response:
        payload = json.loads(response.read(1024 * 1024))
    count = payload.get("ngout", {}).get("data", {}).get("totalCount")
    if payload.get("success") is not True or type(count) is not int or count < 0:
        raise ValueError("NCBI count response has an unexpected schema")
    return count


def fetch_table(collection: str, query: str, fields: list[str], path: Path, role: str, max_bytes: int) -> tuple[dict, list[dict]]:
    count = remote_count(collection, query)
    url = export_url(collection, query, fields)
    retrieve(url, path, max_bytes)
    rows = read_table(path, role)
    if len(rows) != count:
        raise ValueError("NCBI export row count changed or export was truncated; retry into a fresh snapshot")
    entry = source_entry(path, role, url, now(), f"{collection}.tsv")
    entry["expected_rows"] = count
    return entry, rows


def source_entry(path: Path, role: str, url: str, retrieved_at: str, original: str) -> dict:
    return {"role": role, "source_name": "NCBI Pathogen Detection",
            "source_url": url, "retrieved_at": retrieved_at, "filename": path.name,
            "original_filename": original, "sha256": digest(path),
            "parser_version": VERSION, "schema_version": VERSION}


def download(outdir: Path, query: str = DEFAULT_QUERY, max_bytes: int = 128 * 1024 * 1024) -> None:
    if max_bytes < 1 or not query.strip():
        raise ValueError("query must be nonempty and max-bytes positive")
    with new_directory(outdir) as staging:
        entry, rows = fetch_table("ast", query, AST_FIELDS, staging / "asts.tsv", "ast", max_bytes)
        sources = [entry]
        targets = sorted({clean(r["target_acc"]) for r in rows} - {NA})
        if any(not re.fullmatch(r"PDT\d+\.\d+", target) for target in targets):
            raise ValueError("unexpected PDT accession format in AST export")
        batches = [targets[i:i + 100] for i in range(0, len(targets), 100)] or [[]]
        for number, batch in enumerate(batches, 1):
            selection = "target_acc:(" + " OR ".join(batch) + ")" if batch else "target_acc:__empty_selection__"
            raw = staging / f"isolates.{number:04d}.tsv"
            entry, _ = fetch_table("isolates", selection, META_FIELDS + sorted(META_OPTIONAL), raw,
                                   "metadata", max_bytes)
            sources.append(entry)
        write_json(staging / "manifest.json", {"schema_version": VERSION, "parser_version": VERSION,
                   "query": query, "acquisition": "official-table-export", "sources": sources})


def archive(ast: Path, metadata: list[Path], outdir: Path, source_url: str, retrieved_at: str) -> None:
    """Copy externally obtained tables, retaining declared provenance and exact bytes."""
    if not metadata or not source_url.strip():
        raise ValueError("archive requires metadata and its source description/URL")
    if datetime.fromisoformat(retrieved_at.replace("Z", "+00:00")).tzinfo is None:
        raise ValueError("retrieval timestamp must include timezone")
    with new_directory(outdir) as staging:
        sources = []
        for number, (role, source) in enumerate([("ast", ast)] + [("metadata", p) for p in metadata]):
            read_table(source, role)
            target = staging / f"{role}.{number:04d}.tsv"
            shutil.copyfile(source, target)
            sources.append(source_entry(target, role, source_url, retrieved_at, source.name))
        write_json(staging / "manifest.json", {"schema_version": VERSION, "parser_version": VERSION,
                   "acquisition": "local-archive", "sources": sources})


def load_snapshot(snapshot: Path) -> tuple[dict, list[tuple[dict, list[dict]]]]:
    manifest = json.loads((snapshot / "manifest.json").read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("manifest must be a JSON object")
    if manifest.get("schema_version") != VERSION or manifest.get("parser_version") != VERSION:
        raise ValueError("unsupported snapshot schema/parser version")
    sources = manifest.get("sources", [])
    if not isinstance(sources, list) or not sources:
        raise ValueError("manifest sources must be a nonempty list")
    loaded, names, roles = [], set(), Counter()
    required = {"role", "source_name", "source_url", "retrieved_at", "filename", "original_filename",
                "sha256", "parser_version", "schema_version"}
    for entry in sources:
        if not isinstance(entry, dict) or not required <= entry.keys():
            raise ValueError("manifest source is missing required provenance")
        if any(not isinstance(entry[k], str) or not entry[k].strip() for k in required):
            raise ValueError("manifest provenance must contain nonempty strings")
        name, role = entry["filename"], entry["role"]
        if (name in {".", ".."} or Path(name).name != name or "\\" in name or name in names
                or role not in {"ast", "metadata"}):
            raise ValueError("invalid/duplicate manifest filename or role")
        if entry["schema_version"] != VERSION or entry["parser_version"] != VERSION:
            raise ValueError("unsupported source schema/parser version")
        timestamp = datetime.fromisoformat(entry["retrieved_at"].replace("Z", "+00:00"))
        if timestamp.tzinfo is None:
            raise ValueError("retrieval timestamp must include timezone")
        path = snapshot / name
        if path.is_symlink() or digest(path) != entry["sha256"]:
            raise ValueError(f"snapshot checksum mismatch or symlink: {name}")
        names.add(name)
        roles[role] += 1
        records = read_table(path, role)
        if "expected_rows" in entry and (type(entry["expected_rows"]) is not int or entry["expected_rows"] != len(records)):
            raise ValueError("snapshot row count differs from manifest")
        loaded.append((entry, records))
    if roles["ast"] != 1 or not roles["metadata"]:
        raise ValueError("snapshot requires exactly one AST table and at least one metadata table")
    return manifest, loaded


def year_of(value: str) -> str:
    if not re.fullmatch(r"\d{4}(?:-\d{2}(?:-\d{2})?)?", value):
        return NA
    try:
        date.fromisoformat(value + {4: "-01-01", 7: "-01", 10: ""}[len(value)])
    except ValueError:
        return NA
    return value[:4]


def parse_mic(raw: str, sign: str) -> tuple[str, str, bool]:
    if clean(raw) == NA:
        return NA, NA, False
    match = re.fullmatch(r"\s*(<=|>=|==|<|>|=)?\s*(\d+(?:\.\d+)?|\.\d+)\s*", raw)
    supplied = {"==": "="}.get(sign.strip(), sign.strip())
    if not match or supplied not in {"", "NA", "=", "<", "<=", ">", ">="}:
        return NA, NA, True
    embedded, number = match.groups()
    embedded = "=" if embedded == "==" else embedded
    if embedded and supplied not in {"", "NA", embedded}:
        return NA, NA, True
    if Decimal(number) <= 0:
        return NA, NA, True
    return format(Decimal(number).normalize(), "f"), embedded or ("=" if supplied in {"", "NA"} else supplied), False


def load_policy() -> tuple[dict, dict[str, str]]:
    policy = json.loads(POLICY.read_text(encoding="utf-8"))
    drugs = {}
    for name, aliases in policy["drugs"].items():
        for alias in [name] + aliases:
            key = alias.casefold()
            if key in drugs and drugs[key] != name:
                raise ValueError(f"ambiguous drug alias: {alias}")
            drugs[key] = name
    return policy, drugs


def normalize(raw: dict, source: dict, number: int, metadata: list[dict], policy: dict, drugs: dict) -> dict:
    row = {field: NA for field in SCHEMA}
    flags = set()
    row.update(record_id=hashlib.sha256(f'{source["sha256"]}:{number}'.encode()).hexdigest(),
               source_file=source["filename"], source_row=str(number), source_url=source["source_url"],
               source_sha256=source["sha256"], retrieved_at=source["retrieved_at"],
               source_record_id=clean(raw["id"]), raw_record_json=canonical(raw),
               metadata_records_json=canonical(metadata), schema_version=VERSION, policy_version=policy["version"])
    row["isolate_id"] = clean(raw["target_acc"])
    if row["isolate_id"] == NA:
        flags.add("missing_isolate_id")
    if not metadata:
        flags.add("missing_metadata")

    def merged(field: str) -> str:
        values = {clean(r["record"].get(field, "")) for r in metadata}
        values.add(clean(raw.get(field, "")))
        values.discard(NA)
        if len(values) > 1:
            flags.add("metadata_conflict")
            # No arbitrary accession or metadata selection when sources disagree.
            return NA
        return next(iter(values)) if values else NA

    for field, target in {"biosample_acc": "biosample", "asm_acc": "assembly_accession", "Run": "sra_accession",
                          "scientific_name": "organism", "collection_date": "collection_date", "geo_loc_name": "location",
                          "isolation_source": "isolation_source", "host": "host"}.items():
        row[target] = merged(field)
    if row["assembly_accession"] == NA:
        flags.add("missing_assembly")
        if row["sra_accession"] == NA:
            flags.add("missing_genome_link")
    if not re.fullmatch(r"Escherichia coli(?:\s+.+)?", row["organism"]):
        flags.add("organism_mismatch")
    row["antimicrobial_original"] = raw["antibiotic"]
    row["antimicrobial"] = drugs.get(raw["antibiotic"].strip().casefold(), clean(raw["antibiotic"]))
    if raw["antibiotic"].strip().casefold() not in drugs:
        flags.add("unknown_antimicrobial")
    row["phenotype_original"] = raw["phenotype"]
    value = clean(raw["phenotype"])
    row["phenotype"] = NA if value == NA else policy["phenotypes"].get(value.casefold(), "UNKNOWN")
    if row["phenotype"] == NA:
        flags.add("missing_phenotype")
    elif row["phenotype"] == "UNKNOWN":
        flags.add("unknown_phenotype")
    row["mic_raw"], row["mic_secondary_raw"] = raw["mic"], raw["mic_secondary"]
    row["measurement_sign_original"] = raw["measurement_sign"]
    row["mic"], row["mic_operator"], invalid = parse_mic(raw["mic"], raw["measurement_sign"])
    # The official export has separate concentration columns for combination drugs.
    row["mic_secondary"], _, second_invalid = parse_mic(raw["mic_secondary"], "")
    if invalid or second_invalid or (clean(raw["mic_secondary"]) != NA and row["mic"] == NA):
        flags.add("invalid_mic")
    row["mic_units"] = "mg/L" if clean(raw["mic"]) != NA else NA
    row["disk_diffusion_raw"] = raw["disk_diffusion"]
    for field, target in {"method": "testing_method", "standard": "breakpoint_standard",
                          "standard_version": "breakpoint_version", "platform": "platform",
                          "vendor": "vendor", "reagent": "reagent"}.items():
        row[target] = clean(raw.get(field, ""))
    for field in ("testing_method", "breakpoint_standard", "breakpoint_version"):
        if row[field] == NA:
            flags.add("missing_" + field)
    row["collection_year"] = year_of(row["collection_date"])
    if row["collection_year"] == NA:
        flags.add("missing_collection_year")
        if row["collection_date"] != NA:
            flags.add("invalid_collection_date")
    # NCBI geo_loc_name uses the INSDC country[: locality] convention.
    if row["location"] != NA:
        country = row["location"].split(":", 1)[0].strip()
        if country in policy["country_names"]:
            row["country"] = country
        else:
            flags.add("ambiguous_country")
    if row["country"] == NA:
        flags.add("missing_country")
    categories = {policy["host_categories"].get(row["host"].casefold()),
                  policy["source_categories"].get(row["isolation_source"].casefold())} - {None}
    row["source_category"] = next(iter(categories)) if len(categories) == 1 else "unknown"
    if len(categories) > 1:
        flags.add("ambiguous_source")
    if row["host"] == row["isolation_source"] == NA:
        flags.add("missing_source")
    row["qc_flags"] = ";".join(sorted(flags))
    return row


def add_flags(row: dict, flags: set[str]) -> None:
    row["qc_flags"] = ";".join(sorted(set(filter(None, row["qc_flags"].split(";"))) | flags))


def observation(row: dict) -> str:
    raw = json.loads(row["raw_record_json"])
    raw.pop("id", None)
    return canonical(raw)


def annotate_groups(rows: list[dict]) -> list[dict]:
    groups = defaultdict(list)
    for row in rows:
        key = (row["isolate_id"] if row["isolate_id"] != NA else "missing:" + row["record_id"], row["antimicrobial"])
        groups[key].append(row)
    endpoints = []
    for key, group in sorted(groups.items()):
        counts = Counter(observation(row) for row in group)
        phenotypes = {r["phenotype"] for r in group} - {NA, "UNKNOWN", "ND"}
        for row in group:
            flags = set()
            if counts[observation(row)] > 1:
                flags.add("duplicate_record")
            if len(counts) > 1:
                flags.add("multiple_ast_records")
            if len(phenotypes) > 1:
                flags.add("conflicting_ast")
            add_flags(row, flags)
        if all("organism_mismatch" in r["qc_flags"].split(";") for r in group):
            continue
        blocked = {"missing_isolate_id", "missing_assembly", "unknown_antimicrobial", "missing_phenotype",
                   "unknown_phenotype", "invalid_mic", "conflicting_ast", "multiple_ast_records",
                   "metadata_conflict", "organism_mismatch"}
        reasons = sorted({flag for row in group for flag in row["qc_flags"].split(";")} & blocked)
        if any(r["phenotype"] not in {"R", "S", "I", "SDD", "NS", "HLAR"} for r in group):
            reasons = sorted(set(reasons) | {"undefined_phenotype"})
        endpoints.append(dict(zip(ENDPOINT_SCHEMA, [
            group[0]["isolate_id"], key[1], next(iter(phenotypes)) if len(phenotypes) == 1 else NA,
            "false" if reasons else "true", ";".join(reasons), canonical(sorted(r["record_id"] for r in group)),
            str(len(group)), str(len(counts)),
        ])))
    return endpoints


def counts(rows: list[dict], field: str) -> dict:
    return dict(sorted(Counter(row[field] for row in rows).items()))


def unique(rows: list[dict], field: str = "isolate_id") -> set[str]:
    return {r[field] for r in rows} - {NA}


def summarize_rows(rows: list[dict], endpoints: list[dict]) -> tuple[dict, list[dict], list[dict]]:
    ecoli = [r for r in rows if "organism_mismatch" not in r["qc_flags"].split(";")]
    eligible = [r for r in endpoints if r["eligible"] == "true"]
    summary = {"schema_version": VERSION, "raw_ast_records": len(rows), "ecoli_ast_records": len(ecoli),
               "ecoli_unique_isolates": len(unique(ecoli)),
               "assembly_linked_isolates": len(unique([r for r in ecoli if r["assembly_accession"] != NA])),
               "records_with_recognized_phenotype": sum(r["phenotype"] not in {NA, "UNKNOWN", "ND"} for r in ecoli),
               "eligible_isolate_drug_pairs": len(eligible),
               "qc_record_counts_all": dict(sorted(Counter(f for r in rows for f in r["qc_flags"].split(";") if f).items())),
               "ecoli_record_counts": {field: counts(ecoli, field) for field in
                                       ["phenotype", "antimicrobial", "collection_year", "country", "source_category"]},
               "endpoint_exclusion_counts_overlapping": dict(sorted(Counter(f for r in endpoints
                          for f in r["exclusion_reasons"].split(";") if f).items())),
               "flow": [{"stage": "raw_ast", "unit": "records", "count": len(rows)},
                        {"stage": "organism_excluded", "unit": "records", "count": len(rows) - len(ecoli)},
                        {"stage": "ecoli_retained", "unit": "records", "count": len(ecoli)},
                        {"stage": "ecoli_pairs", "unit": "isolate_drug_groups", "count": len(endpoints)},
                        {"stage": "ineligible_pairs", "unit": "isolate_drug_groups", "count": len(endpoints) - len(eligible)},
                        {"stage": "eligible_pairs", "unit": "isolate_drug_groups", "count": len(eligible)}]}
    by_isolate, by_drug, eligible_by_drug = defaultdict(list), defaultdict(list), defaultdict(list)
    for row in ecoli:
        by_isolate[row["isolate_id"]].append(row)
        by_drug[row["antimicrobial"]].append(row)
    for row in eligible:
        eligible_by_drug[row["antimicrobial"]].append(row)
    isolates = []
    for isolate in sorted(unique(ecoli)):
        records = by_isolate[isolate]
        item = {"isolate_id": isolate, "ast_records": str(len(records))}
        for field in ["biosample", "assembly_accession", "sra_accession", "collection_year", "country", "source_category"]:
            item[field + "_values_json"] = canonical(sorted(unique(records, field)))
        isolates.append(item)
    drugs = []
    for drug, records in sorted(by_drug.items()):
        ready = eligible_by_drug[drug]
        r_count, s_count = sum(r["phenotype"] == "R" for r in ready), sum(r["phenotype"] == "S" for r in ready)
        years = sorted(unique(records, "collection_year"))
        drugs.append({"antimicrobial": drug, "ast_records": len(records), "isolates": len(unique(records)),
                      "eligible_pairs": len(ready), "R": r_count, "S": s_count, "other": len(ready) - r_count - s_count,
                      "R_fraction_of_RS": f"{r_count / (r_count + s_count):.6f}" if r_count + s_count else NA,
                      "countries": len(unique(records, "country")), "year_min": years[0] if years else NA,
                      "year_max": years[-1] if years else NA,
                      **{name: len(unique([r for r in records if r[field] not in {NA, "unknown"}])) for name, field in
                         [("assembly_linked_isolates", "assembly_accession"), ("country_known_isolates", "country"),
                          ("year_known_isolates", "collection_year"), ("source_classified_isolates", "source_category")]}})
    return summary, isolates, drugs


ISOLATE_SCHEMA = ["isolate_id", "ast_records"] + [f + "_values_json" for f in
                  ["biosample", "assembly_accession", "sra_accession", "collection_year", "country", "source_category"]]
DRUG_SCHEMA = ("antimicrobial ast_records isolates eligible_pairs R S other R_fraction_of_RS countries year_min year_max "
               "assembly_linked_isolates country_known_isolates year_known_isolates source_classified_isolates").split()


def prepare(snapshot: Path) -> tuple[list[dict], list[dict], dict]:
    manifest, loaded = load_snapshot(snapshot)
    policy, drugs = load_policy()
    index = defaultdict(list)
    for source, raw_rows in loaded:
        if source["role"] == "metadata":
            for number, raw in enumerate(raw_rows, 1):
                key = clean(raw["target_acc"])
                if key != NA:
                    index[key].append({"source_file": source["filename"], "source_row": number,
                                       "source_sha256": source["sha256"], "source_url": source["source_url"],
                                       "retrieved_at": source["retrieved_at"], "record": raw})
    rows = []
    for source, raw_rows in loaded:
        if source["role"] == "ast":
            rows = [normalize(raw, source, number, index.get(clean(raw["target_acc"]), []), policy, drugs)
                    for number, raw in enumerate(raw_rows, 1)]
    endpoints = annotate_groups(rows)
    summary, _, _ = summarize_rows(rows, endpoints)
    summary.update(policy_version=policy["version"], policy_sha256=digest(POLICY),
                   snapshot_manifest_sha256=digest(snapshot / "manifest.json"),
                   source_files=len(manifest["sources"]),
                   raw_metadata_records=sum(len(raw) for source, raw in loaded if source["role"] == "metadata"))
    return rows, endpoints, summary


def build(snapshot: Path, outdir: Path) -> dict:
    if snapshot.resolve() in outdir.resolve().parents:
        raise ValueError("processed output must be outside the raw snapshot")
    rows, endpoints, summary = prepare(snapshot)
    _, isolates, drugs = summarize_rows(rows, endpoints)
    with new_directory(outdir) as staging:
        write_tsv(staging / "ast.normalized.tsv", SCHEMA, rows)
        write_tsv(staging / "ecoli_ast.tsv", SCHEMA,
                  [r for r in rows if "organism_mismatch" not in r["qc_flags"].split(";")])
        write_tsv(staging / "ecoli_isolates.tsv", ISOLATE_SCHEMA, isolates)
        write_tsv(staging / "ecoli_endpoints.tsv", ENDPOINT_SCHEMA, endpoints)
        write_tsv(staging / "qc_flags.tsv", ["record_id", "flag"],
                  [{"record_id": row["record_id"], "flag": flag} for row in rows
                   for flag in row["qc_flags"].split(";") if flag])
        write_tsv(staging / "antimicrobial_summary.tsv", DRUG_SCHEMA, drugs)
        write_json(staging / "cohort_summary.json", summary)
    return summary


def validate(snapshot: Path, processed: Path | None = None) -> dict:
    """Verify raw hashes/schema and, optionally, every reproducible output byte."""
    _, _, summary = prepare(snapshot)
    if processed is not None:
        with tempfile.TemporaryDirectory() as directory:
            rebuilt = Path(directory) / "processed"
            build(snapshot, rebuilt)
            expected = {p.name for p in rebuilt.iterdir()}
            if {p.name for p in processed.iterdir()} != expected:
                raise ValueError("processed output file set differs from rebuild")
            for name in sorted(expected):
                if digest(rebuilt / name) != digest(processed / name):
                    raise ValueError(f"processed output differs from rebuild: {name}")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("download", help="fetch official AST/metadata tables into a new raw directory")
    command.add_argument("--outdir", type=Path, required=True)
    command.add_argument("--query", default=DEFAULT_QUERY)
    command.add_argument("--max-bytes", type=int, default=128 * 1024 * 1024)
    command = commands.add_parser("archive", help="archive existing exports or documented fixtures")
    command.add_argument("--ast", type=Path, required=True)
    command.add_argument("--metadata", type=Path, nargs="+", required=True)
    command.add_argument("--source-url", required=True)
    command.add_argument("--retrieved-at", required=True)
    command.add_argument("--outdir", type=Path, required=True)
    for name in ["build", "validate", "summarize"]:
        command = commands.add_parser(name)
        command.add_argument("--snapshot", type=Path, required=True)
        if name == "build":
            command.add_argument("--outdir", type=Path, required=True)
        elif name == "validate":
            command.add_argument("--processed", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "download":
            download(args.outdir, args.query, args.max_bytes)
        elif args.command == "archive":
            archive(args.ast, args.metadata, args.outdir, args.source_url, args.retrieved_at)
        else:
            summary = build(args.snapshot, args.outdir) if args.command == "build" else validate(
                args.snapshot, args.processed if args.command == "validate" else None)
            print(json.dumps(summary, indent=2, sort_keys=True))
    except (ValueError, KeyError, TypeError, OSError, csv.Error) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

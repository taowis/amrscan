#!/usr/bin/env python3
"""Build deterministic, phenotype-only E. coli AST cohort census tables."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

import ast_dataset

VERSION = "1.0.0"
DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "benchmark/cohorts/ecoli_research_v1/cohort_config.json"
NA = "NA"
PDT = re.compile(r"PDT\d+\.\d+")
BIOSAMPLE = re.compile(r"SAM[END]\d+")
COUNTRY_MISSING = {"", "NA", "unknown", "UNKNOWN"}

OBSERVATION_FIELDS = ast_dataset.SCHEMA
CANDIDATE_FIELDS = (
    "cohort_id isolate_id pdt_accession biosample_accession assembly_accession organism strain taxid country "
    "location collection_date collection_year host isolation_source source_category submitter project snp_cluster "
    "snp_clusters_json ast_record_count antimicrobial_count metadata_row_count identity_status identity_conflicts metadata_conflicts "
    "source_file source_sha256 source_url retrieved_at metadata_source_hashes_json source_record_ids_json "
    "raw_metadata_records_json"
).split()
ENDPOINT_FIELDS = (
    "cohort_id isolate_id pdt_accession biosample_accession assembly_accession antimicrobial antimicrobial_original "
    "phenotype_endpoint phenotype_status n_ast_observations n_R n_S n_I n_other duplicate_status conflict_status "
    "country collection_year source_category assembly_available eligible_binary_rs eligible_overall "
    "eligible_country_analysis eligible_temporal_analysis eligible_source_analysis exclusion_reason "
    "source_record_ids_json"
).split()
EXCLUSION_FIELDS = (
    "cohort_id analysis isolate_id pdt_accession antimicrobial reason source_record_ids_json"
).split()
ANTIMICROBIAL_FIELDS = (
    "antimicrobial total_unique_isolates binary_rs_eligible_isolates n_R n_S n_I n_other n_conflict "
    "resistance_fraction n_countries country_available_n country_available_fraction earliest_year latest_year "
    "year_span year_available_n year_available_fraction n_source_categories human_n animal_n food_n "
    "environment_n other_source_n unknown_source_n assembly_linked_n assembly_linked_fraction conflict_fraction "
    "overall_ready geography_ready temporal_ready source_ready readiness_reason"
).split()
COUNTRY_FIELDS = "country isolate_count endpoint_count n_R n_S antimicrobial_count earliest_year latest_year".split()
YEAR_FIELDS = "collection_year isolate_count endpoint_count n_R n_S antimicrobial_count".split()
SOURCE_FIELDS = "source_category antimicrobial isolate_count endpoint_count n_R n_S".split()
MISSING_FIELDS = "unit field n_total n_missing fraction_missing".split()
CONFLICT_FIELDS = (
    "isolate_id pdt_accession antimicrobial phenotype_endpoint conflict_status n_ast_observations "
    "phenotypes_json source_record_ids_json"
).split()
ELIGIBLE_FIELDS = CANDIDATE_FIELDS
INVENTORY_FIELDS = "pdt_accession biosample_accession assembly_accession organism strain taxid country collection_year project snp_clusters_json source_rows source_sha256".split()


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def missing(value: object) -> bool:
    return value is None or str(value).strip() in {"", NA}


def unique_values(values) -> list[str]:
    return sorted({str(value).strip() for value in values if not missing(value)})


def consensus(values) -> tuple[str, bool]:
    vals = unique_values(values)
    return (vals[0], False) if len(vals) == 1 else (NA, len(vals) > 1)


def is_ecoli(value: str) -> bool:
    name = value.strip().casefold()
    return name == "escherichia coli" or name.startswith("escherichia coli ") or name.startswith("escherichia coli(")


def valid_assembly(value: str, config: dict) -> bool:
    return bool(re.fullmatch(config["identity_policy"]["assembly_pattern"], value or ""))


def parse_json_cell(row: dict[str, str], field: str, expected):
    try:
        value = json.loads(row[field])
    except (KeyError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid {field} JSON in AST record {row.get('record_id', NA)}") from error
    if not isinstance(value, expected):
        raise ValueError(f"{field} must contain {expected.__name__} values")
    return value


def metadata_from_rows(rows: list[dict[str, str]]) -> list[dict]:
    """Return unique metadata source rows referenced by normalized AST records."""
    found: dict[tuple[str, str, str], dict] = {}
    for row in rows:
        for item in parse_json_cell(row, "metadata_records_json", list):
            if not isinstance(item, dict) or not isinstance(item.get("record"), dict):
                raise ValueError(f"malformed metadata provenance in {row.get('record_id', NA)}")
            key = (str(item.get("source_sha256", "")), str(item.get("source_file", "")),
                   str(item.get("source_row", "")))
            found[key] = item
    return [found[key] for key in sorted(found)]


def metadata_values(records: list[dict], key: str, fallback: str | None = None) -> list[str]:
    values = []
    for item in records:
        record = item["record"]
        value = record.get(key, "")
        if missing(value) and fallback:
            value = record.get(fallback, "")
        if isinstance(value, (dict, list)):
            value = canonical(value)
        values.append(str(value))
    return values


def make_candidate(cohort_id: str, isolate_id: str, rows: list[dict[str, str]], config: dict) -> tuple[dict, list[dict]]:
    metadata = metadata_from_rows(rows)
    raw_values = [parse_json_cell(row, "raw_record_json", dict) for row in rows]
    conflicts: list[str] = []

    def field(name: str, metadata_key: str | None = None, raw_key: str | None = None) -> str:
        values = [row.get(name, NA) for row in rows]
        if metadata_key:
            values.extend(metadata_values(metadata, metadata_key))
        if raw_key:
            values.extend(raw.get(raw_key, "") for raw in raw_values)
        value, conflict = consensus(values)
        if conflict:
            conflicts.append(name)
        return value

    biosample = field("biosample", "biosample_acc", "biosample_acc")
    assembly = field("assembly_accession", "asm_acc")
    organism = field("organism", "scientific_name", "scientific_name")
    strain = field("strain", "strain")
    taxid = field("taxid", "taxid", "taxid")
    if taxid == NA:
        taxid = field("taxid", "species_taxid", "species_taxid")
    submitter = field("submitter", "bioproject_center", "submitter")
    project = field("project", "bioproject_acc", "bioproject_acc")
    clusters = unique_values(metadata_values(metadata, "erd_group"))
    snp_cluster = clusters[0] if len(clusters) == 1 else NA
    if len(clusters) > 1:
        conflicts.append("snp_cluster")
    country = field("country")
    location = field("location")
    collection_date = field("collection_date")
    collection_year = field("collection_year")
    host = field("host")
    isolation_source = field("isolation_source")
    source_category = field("source_category")
    identity_conflicts = sorted(set(conflicts) & {
        "biosample", "assembly_accession", "organism", "taxid"
    })
    identity = config["identity_policy"]
    valid_identity = ((not identity["require_versioned_pdt"] or bool(PDT.fullmatch(isolate_id)))
                      and (not identity["require_biosample_for_overall"] or bool(BIOSAMPLE.fullmatch(biosample)))
                      and (not identity["require_genbank_assembly_for_overall"] or valid_assembly(assembly, config))
                      and is_ecoli(organism) and not identity_conflicts)
    identity_status = "verified" if valid_identity else "incomplete_or_conflicting"
    ast_files = unique_values(row.get("source_file", NA) for row in rows)
    ast_hashes = unique_values(row.get("source_sha256", NA) for row in rows)
    ast_urls = unique_values(row.get("source_url", NA) for row in rows)
    retrieved = unique_values(row.get("retrieved_at", NA) for row in rows)
    metadata_hashes = unique_values(item.get("source_sha256", "") for item in metadata)
    record_ids = sorted(str(row.get("source_record_id", NA)) for row in rows)
    drugs = unique_values(row.get("antimicrobial", NA) for row in rows)
    candidate = {
        "cohort_id": cohort_id, "isolate_id": isolate_id, "pdt_accession": isolate_id,
        "biosample_accession": biosample, "assembly_accession": assembly, "organism": organism,
        "strain": strain, "taxid": taxid, "country": country, "location": location,
        "collection_date": collection_date, "collection_year": collection_year, "host": host,
        "isolation_source": isolation_source, "source_category": source_category,
        "submitter": submitter, "project": project, "snp_cluster": snp_cluster,
        "snp_clusters_json": canonical(clusters), "ast_record_count": str(len(rows)),
        "antimicrobial_count": str(len(drugs)), "metadata_row_count": str(len(metadata)),
        "identity_status": identity_status, "identity_conflicts": ";".join(identity_conflicts) or NA,
        "metadata_conflicts": ";".join(sorted(set(conflicts))) or NA,
        "source_file": ast_files[0] if len(ast_files) == 1 else canonical(ast_files),
        "source_sha256": ast_hashes[0] if len(ast_hashes) == 1 else canonical(ast_hashes),
        "source_url": ast_urls[0] if len(ast_urls) == 1 else canonical(ast_urls),
        "retrieved_at": retrieved[0] if len(retrieved) == 1 else canonical(retrieved),
        "metadata_source_hashes_json": canonical(metadata_hashes),
        "source_record_ids_json": canonical(record_ids),
        "raw_metadata_records_json": canonical(metadata),
    }
    return candidate, raw_values


def derive_endpoint(cohort_id: str, isolate_id: str, antimicrobial: str,
                    rows: list[dict[str, str]], candidate: dict | None, config: dict) -> dict:
    phenotypes = [row.get("phenotype", NA) if not missing(row.get("phenotype", NA)) else NA for row in rows]
    counts = Counter(phenotypes)
    distinct = set(phenotypes)
    binary_categories = set(config["phenotype_policy"]["binary_categories"])
    recognized = binary_categories | set(config["phenotype_policy"]["recognized_nonbinary_categories"])
    known_drug = not any("unknown_antimicrobial" in row.get("qc_flags", "").split(";") for row in rows)
    if len(distinct) == 1 and next(iter(distinct)) in binary_categories:
        endpoint = next(iter(distinct))
        status = "single_binary_observation" if len(rows) == 1 else "repeated_consistent"
        conflict = "none"
    elif binary_categories <= distinct:
        endpoint, status, conflict = "CONFLICT", "discordant_rs", "discordant_rs"
    elif len(distinct) > 1:
        endpoint, status, conflict = "NON_BINARY_OR_CONFLICT", "mixed_phenotypes", "mixed_phenotypes"
    elif distinct == {NA}:
        endpoint, status, conflict = "NA", "missing_phenotype", "none"
    else:
        endpoint = next(iter(distinct))
        status = "non_binary" if endpoint in recognized else "unknown_phenotype"
        conflict = "none"
    repeated = len(rows) > 1
    has_exact_duplicate = any("duplicate_record" in row.get("qc_flags", "").split(";") for row in rows)
    duplicate_status = ("repeated_consistent" if repeated and len(distinct) == 1 else
                        "repeated_discordant" if repeated else "single_observation")
    if has_exact_duplicate:
        duplicate_status += ";exact_duplicate_source_rows"
    ids = sorted(str(row.get("source_record_id", NA)) for row in rows)
    identity = config["identity_policy"]
    valid_pdt = not identity["require_versioned_pdt"] or bool(PDT.fullmatch(isolate_id))
    pdt = isolate_id if PDT.fullmatch(isolate_id) else NA
    biosample = candidate["biosample_accession"] if candidate else NA
    assembly = candidate["assembly_accession"] if candidate else NA
    country = candidate["country"] if candidate else NA
    year = candidate["collection_year"] if candidate else NA
    source = candidate["source_category"] if candidate else NA
    identity_ok = bool(candidate and candidate["identity_status"] == "verified")
    assembly_ok = bool(candidate and valid_assembly(assembly, config))
    binary = (valid_pdt and known_drug and endpoint in binary_categories
              and conflict == "none")
    overall = binary and identity_ok and assembly_ok
    country_ok = overall and country not in COUNTRY_MISSING
    temporal_ok = overall and bool(re.fullmatch(r"\d{4}", year or ""))
    known_sources = set(config["eligibility"]["known_source_categories"])
    source_ok = overall and source in known_sources
    reasons = []
    if not valid_pdt:
        reasons.append("missing_or_invalid_isolate_id")
    if not known_drug:
        reasons.append("unknown_antimicrobial")
    if conflict == "discordant_rs":
        reasons.append("conflicting_rs")
    elif endpoint not in binary_categories:
        reasons.append("non_binary_phenotype")
    if not identity_ok:
        reasons.append("identity_incomplete_or_conflicting")
    if not assembly_ok:
        reasons.append("missing_or_invalid_assembly")
    if not country_ok:
        reasons.append("missing_country" if overall else "overall_ineligible")
    if not temporal_ok:
        reasons.append("missing_or_invalid_collection_year" if overall else "overall_ineligible")
    if not source_ok:
        reasons.append("missing_or_unknown_source_category" if overall else "overall_ineligible")
    return {
        "cohort_id": cohort_id, "isolate_id": isolate_id, "pdt_accession": pdt,
        "biosample_accession": biosample, "assembly_accession": assembly,
        "antimicrobial": antimicrobial, "antimicrobial_original": canonical(unique_values(
            row.get("antimicrobial_original", NA) for row in rows)),
        "phenotype_endpoint": endpoint, "phenotype_status": status,
        "n_ast_observations": str(len(rows)), "n_R": str(counts["R"]), "n_S": str(counts["S"]),
        "n_I": str(counts["I"]), "n_other": str(sum(n for p, n in counts.items() if p not in {"R", "S", "I"})),
        "duplicate_status": duplicate_status, "conflict_status": conflict,
        "country": country, "collection_year": year, "source_category": source,
        "assembly_available": "true" if assembly_ok else "false",
        "eligible_binary_rs": "true" if binary else "false", "eligible_overall": "true" if overall else "false",
        "eligible_country_analysis": "true" if country_ok else "false",
        "eligible_temporal_analysis": "true" if temporal_ok else "false",
        "eligible_source_analysis": "true" if source_ok else "false",
        "exclusion_reason": ";".join(sorted(set(reasons))) if reasons else NA,
        "source_record_ids_json": canonical(ids),
    }


def append_exclusions(endpoint: dict, output: list[dict], config: dict) -> None:
    dimensions = {
        "binary_rs": endpoint["eligible_binary_rs"],
        "overall": endpoint["eligible_overall"],
        "country": endpoint["eligible_country_analysis"],
        "temporal": endpoint["eligible_temporal_analysis"],
        "source": endpoint["eligible_source_analysis"],
    }
    reasons = set()
    if not PDT.fullmatch(endpoint["isolate_id"]):
        reasons.add("missing_or_invalid_isolate_id")
    if endpoint["conflict_status"] == "discordant_rs":
        reasons.add("conflicting_rs")
    elif endpoint["phenotype_endpoint"] not in {"R", "S"}:
        reasons.add("non_binary_phenotype")
    if "unknown_antimicrobial" in endpoint["exclusion_reason"]:
        reasons.add("unknown_antimicrobial")
    binary_reasons = set(reasons)
    overall_reasons = set(binary_reasons)
    if endpoint["exclusion_reason"] != NA and "identity_incomplete_or_conflicting" in endpoint["exclusion_reason"]:
        overall_reasons.add("identity_incomplete_or_conflicting")
    if endpoint["assembly_available"] != "true":
        overall_reasons.add("missing_or_invalid_assembly")
    known_sources = set(config["eligibility"]["known_source_categories"])
    dimension_reasons = {
        "binary_rs": binary_reasons,
        "overall": overall_reasons,
        "country": overall_reasons | ({"missing_country"} if endpoint["country"] in COUNTRY_MISSING else set()),
        "temporal": overall_reasons | ({"missing_or_invalid_collection_year"}
                                          if not re.fullmatch(r"\d{4}", endpoint["collection_year"] or "") else set()),
        "source": overall_reasons | ({"missing_or_unknown_source_category"}
                                        if endpoint["source_category"] not in known_sources
                                        else set()),
    }
    for analysis, eligible in dimensions.items():
        if eligible != "true":
            relevant = sorted(dimension_reasons[analysis])
            output.append({"cohort_id": endpoint["cohort_id"], "analysis": analysis,
                           "isolate_id": endpoint["isolate_id"], "pdt_accession": endpoint["pdt_accession"],
                           "antimicrobial": endpoint["antimicrobial"],
                           "reason": ";".join(relevant) or "ineligible",
                           "source_record_ids_json": endpoint["source_record_ids_json"]})


def groups_by(rows: list[dict], keyfunc) -> dict:
    result = defaultdict(list)
    for row in rows:
        result[keyfunc(row)].append(row)
    return result


def readiness(stats: dict, config: dict) -> tuple[dict[str, str], str]:
    policy = config.get("study_drug_readiness", {})
    thresholds = [policy.get(k) for k in ("minimum_total_isolates", "minimum_resistant_isolates",
                                           "minimum_susceptible_isolates", "minimum_countries",
                                           "minimum_year_span", "minimum_source_categories",
                                           "maximum_conflict_fraction")]
    if policy.get("decision") == "HUMAN_REVIEW_REQUIRED" or any(x is None for x in thresholds):
        state = "HUMAN_REVIEW_REQUIRED"
        return {k: state for k in ("overall_ready", "geography_ready", "temporal_ready", "source_ready")}, \
            "study-drug thresholds require human scientific review"
    ready_overall = (stats["binary_rs_eligible_isolates"] >= policy["minimum_total_isolates"]
                     and stats["n_R"] >= policy["minimum_resistant_isolates"]
                     and stats["n_S"] >= policy["minimum_susceptible_isolates"]
                     and stats["conflict_fraction_value"] <= policy["maximum_conflict_fraction"])
    flags = {
        "overall_ready": ready_overall,
        "geography_ready": ready_overall and stats["n_countries"] >= policy["minimum_countries"],
        "temporal_ready": ready_overall and stats["year_span_value"] >= policy["minimum_year_span"],
        "source_ready": ready_overall and stats["n_source_categories"] >= policy["minimum_source_categories"],
    }
    return {k: "true" if v else "false" for k, v in flags.items()}, "configured thresholds"


def summarize_antimicrobials(endpoints: list[dict], config: dict) -> list[dict]:
    by_drug = groups_by(endpoints, lambda row: row["antimicrobial"])
    output = []
    for drug, rows in sorted(by_drug.items()):
        binary = [r for r in rows if r["eligible_binary_rs"] == "true"]
        rs_counts = Counter(r["phenotype_endpoint"] for r in binary)
        countries = unique_values(r["country"] for r in rows if r["country"] not in COUNTRY_MISSING)
        country_n = len({r["isolate_id"] for r in rows if r["country"] not in COUNTRY_MISSING})
        years = sorted({int(r["collection_year"]) for r in rows if re.fullmatch(r"\d{4}", r["collection_year"] or "")})
        year_n = len({r["isolate_id"] for r in rows if re.fullmatch(r"\d{4}", r["collection_year"] or "")})
        source_counts = Counter()
        for category in {r["source_category"] for r in rows}:
            source_counts[category] = len({r["isolate_id"] for r in rows
                                           if r["source_category"] == category and r["isolate_id"] != NA})
        conflict_n = sum(r["conflict_status"] != "none" for r in rows)
        assembly_n = len({r["isolate_id"] for r in rows if r["assembly_available"] == "true"})
        total_isolates = len({r["isolate_id"] for r in rows})
        stats = {
            "binary_rs_eligible_isolates": len({r["isolate_id"] for r in binary}),
            "n_R": rs_counts["R"], "n_S": rs_counts["S"], "n_countries": len(countries),
            "n_source_categories": len({r["source_category"] for r in rows
                                         if r["source_category"] in config["eligibility"]["known_source_categories"]}),
            "conflict_fraction_value": conflict_n / len(rows) if rows else 0,
            "year_span_value": years[-1] - years[0] if len(years) > 1 else 0,
        }
        states, reason = readiness(stats, config)
        output.append({
            "antimicrobial": drug, "total_unique_isolates": str(total_isolates),
            "binary_rs_eligible_isolates": str(stats["binary_rs_eligible_isolates"]),
            "n_R": str(rs_counts["R"]), "n_S": str(rs_counts["S"]),
            "n_I": str(sum(r["phenotype_endpoint"] == "I" for r in rows)),
            "n_other": str(sum(r["phenotype_endpoint"] not in {"R", "S", "I"} for r in rows)),
            "n_conflict": str(conflict_n),
            "resistance_fraction": (f"{rs_counts['R'] / (rs_counts['R'] + rs_counts['S']):.8f}"
                                    if rs_counts["R"] + rs_counts["S"] else NA),
            "n_countries": str(stats["n_countries"]), "country_available_n": str(country_n),
            "country_available_fraction": f"{country_n / total_isolates:.8f}" if total_isolates else "0",
            "earliest_year": str(years[0]) if years else NA, "latest_year": str(years[-1]) if years else NA,
            "year_span": str(stats["year_span_value"]), "year_available_n": str(year_n),
            "year_available_fraction": f"{year_n / total_isolates:.8f}" if total_isolates else "0",
            "n_source_categories": str(stats["n_source_categories"]),
            "human_n": str(source_counts["human"]), "animal_n": str(source_counts["animal"]),
            "food_n": str(source_counts["food"]), "environment_n": str(source_counts["environment"]),
            "other_source_n": str(source_counts["other"]), "unknown_source_n": str(source_counts["unknown"]),
            "assembly_linked_n": str(assembly_n),
            "assembly_linked_fraction": f"{assembly_n / total_isolates:.8f}" if total_isolates else "0",
            "conflict_fraction": f"{conflict_n / len(rows):.8f}" if rows else "0",
            **states, "readiness_reason": reason,
        })
    return output


def summarize_groups(endpoints: list[dict], field: str, schema: list[str]) -> list[dict]:
    by_value = groups_by(endpoints, lambda r: r[field] if not missing(r[field]) else NA)
    output = []
    for value, rows in sorted(by_value.items()):
        eligible = [r for r in rows if r["eligible_binary_rs"] == "true"]
        years = sorted({r["collection_year"] for r in rows if re.fullmatch(r"\d{4}", r["collection_year"] or "")})
        output.append({
            field: value, "isolate_count": str(len({r["isolate_id"] for r in rows if r["isolate_id"] != NA})),
            "endpoint_count": str(len(rows)), "n_R": str(sum(r["phenotype_endpoint"] == "R" for r in eligible)),
            "n_S": str(sum(r["phenotype_endpoint"] == "S" for r in eligible)),
            "antimicrobial_count": str(len({r["antimicrobial"] for r in rows})),
            "earliest_year": years[0] if years else NA, "latest_year": years[-1] if years else NA,
        })
    return output


def source_summary(endpoints: list[dict]) -> list[dict]:
    by_group = groups_by(endpoints, lambda r: (r["source_category"], r["antimicrobial"]))
    result = []
    for (category, drug), rows in sorted(by_group.items()):
        eligible = [r for r in rows if r["eligible_binary_rs"] == "true"]
        result.append({"source_category": category, "antimicrobial": drug,
                       "isolate_count": str(len({r["isolate_id"] for r in rows if r["isolate_id"] != NA})),
                       "endpoint_count": str(len(rows)),
                       "n_R": str(sum(r["phenotype_endpoint"] == "R" for r in eligible)),
                       "n_S": str(sum(r["phenotype_endpoint"] == "S" for r in eligible))})
    return result


def missingness(observations: list[dict], candidates: list[dict]) -> list[dict]:
    result = []
    observation_fields = ("biosample", "assembly_accession", "country", "collection_year", "host",
                          "isolation_source", "source_category", "breakpoint_standard", "phenotype", "antimicrobial")
    candidate_fields = ("biosample_accession", "assembly_accession", "country", "collection_year", "host",
                        "isolation_source", "source_category", "strain", "taxid", "submitter", "project", "snp_cluster")
    for unit, rows, fields in (("E. coli AST observations", observations, observation_fields),
                               ("E. coli candidate isolates", candidates, candidate_fields)):
        for field in fields:
            n_missing = sum(missing(row.get(field, NA)) for row in rows)
            result.append({"unit": unit, "field": field, "n_total": str(len(rows)), "n_missing": str(n_missing),
                           "fraction_missing": f"{n_missing / len(rows):.8f}" if rows else "0"})
    return result


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def validate_config(path: Path) -> dict:
    config = json.loads(path.read_text(encoding="utf-8"))
    required = {"cohort_id", "config_version", "organism", "identity_policy", "phenotype_policy",
                "duplicate_policy", "eligibility", "study_drug_readiness"}
    if not isinstance(config, dict) or not required <= config.keys():
        raise ValueError("cohort configuration is missing required keys")
    if config["organism"] != "Escherichia coli" or not config["cohort_id"]:
        raise ValueError("this builder requires an explicit Escherichia coli cohort id")
    try:
        re.compile(config["identity_policy"]["assembly_pattern"])
    except (KeyError, re.error) as error:
        raise ValueError("identity_policy must define a valid assembly_pattern") from error
    if set(config["phenotype_policy"].get("binary_categories", [])) != {"R", "S"}:
        raise ValueError("binary phenotype policy must retain exactly R and S")
    if not config["eligibility"].get("known_source_categories"):
        raise ValueError("eligibility must define known source categories")
    return config


def inventory(candidate_path: Path, outdir: Path, source_url: str, query: str, retrieved_at: str) -> dict:
    """Summarize the archived PR #4 metadata candidate export, without implying AST rows."""
    required = ["target_acc", "biosample_acc", "asm_acc", "scientific_name", "strain", "erd_group", "bioproject_acc"]
    with candidate_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", strict=True)
        header = next(reader, [])
        if header:
            header[0] = header[0].removeprefix("#")
        if len(set(header)) != len(header):
            raise ValueError("candidate export has duplicate columns")
        missing_columns = set(required) - set(header)
        if missing_columns:
            raise ValueError(f"candidate export missing columns: {sorted(missing_columns)}")
        positions = {name: header.index(name) for name in required + ["taxid"] if name in header}
        grouped = defaultdict(list)
        rows_n = 0
        for ordinal, values in enumerate(reader, 2):
            if len(values) != len(header):
                raise ValueError(f"candidate export line {ordinal}: wrong column count")
            row = {name: values[index] for name, index in positions.items()}
            if not PDT.fullmatch(row["target_acc"]):
                raise ValueError(f"candidate export line {ordinal}: invalid target_acc")
            grouped[row["target_acc"]].append(row)
            rows_n += 1
    source_hash = digest(candidate_path)
    candidates, project_counts, identity_conflicts = [], Counter(), Counter()
    missing_by_field = Counter()
    for pdt, group in sorted(grouped.items(), key=lambda item: tuple(map(int, PDT.fullmatch(item[0]).groups()))):
        by_field = {name: unique_values(row.get(name, "") for row in group)
                    for name in ("biosample_acc", "asm_acc", "scientific_name", "strain", "bioproject_acc", "taxid")}
        for key in ("biosample_acc", "asm_acc", "scientific_name", "bioproject_acc"):
            if len(by_field[key]) > 1:
                identity_conflicts[key] += 1
        for key in ("strain", "erd_group"):
            vals = unique_values(row.get(key, "") for row in group)
            if not vals:
                missing_by_field[key] += 1
        biosample = by_field["biosample_acc"][0] if len(by_field["biosample_acc"]) == 1 else NA
        assembly = by_field["asm_acc"][0] if len(by_field["asm_acc"]) == 1 else NA
        organism = by_field["scientific_name"][0] if len(by_field["scientific_name"]) == 1 else NA
        strain = by_field["strain"][0] if len(by_field["strain"]) == 1 else NA
        project = by_field["bioproject_acc"][0] if len(by_field["bioproject_acc"]) == 1 else NA
        clusters = unique_values(row.get("erd_group", "") for row in group)
        if project != NA:
            project_counts[project] += 1
        taxids = unique_values(row.get("taxid", "") for row in group)
        candidates.append({
            "pdt_accession": pdt, "biosample_accession": biosample, "assembly_accession": assembly,
            "organism": organism, "strain": strain, "taxid": taxids[0] if len(taxids) == 1 else NA,
            "country": NA, "collection_year": NA, "project": project,
            "snp_clusters_json": canonical(clusters), "source_rows": str(len(group)), "source_sha256": source_hash,
        })
    all_identity = {
        field: {"unique_values": len({c[field] for c in candidates if c[field] != NA}),
                "missing_unique_pdt": sum(c[field] == NA for c in candidates)}
        for field in ("biosample_accession", "assembly_accession", "organism", "project")
    }
    try:
        datetime_value = datetime_from_iso(retrieved_at)
    except ValueError as error:
        raise ValueError("retrieved-at must be an ISO-8601 timestamp with timezone") from error
    summary = {
        "snapshot_type": "metadata candidate inventory only; not an AST observation census",
        "source_url": source_url, "query": query, "retrieved_at": datetime_value,
        "source_filename": candidate_path.name, "source_sha256": source_hash,
        "source_row_count": rows_n, "unique_versioned_pdt_count": len(candidates),
        "candidate_rows_per_pdt": dict(sorted(Counter(c["source_rows"] for c in candidates).items())),
        "identity_field_summary": all_identity,
        "identity_conflicting_pdt_counts": dict(sorted(identity_conflicts.items())),
        "missing_unique_pdt_counts": dict(sorted(missing_by_field.items())),
        "projects": [{"bioproject_acc": project, "unique_pdt_count": count}
                     for project, count in sorted(project_counts.items(), key=lambda item: (-item[1], item[0]))],
        "ast_observation_count": None,
        "endpoint_census_status": "unavailable_without the matching full AST observation export",
    }
    outdir = Path(outdir)
    if outdir.exists() or outdir.is_symlink():
        raise ValueError(f"output already exists: {outdir}")
    outdir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".ecoli-inventory-", dir=outdir.parent))
    try:
        write_tsv(staging / "source_candidates.tsv", INVENTORY_FIELDS, candidates)
        write_tsv(staging / "project_summary.tsv", ["bioproject_acc", "unique_pdt_count"], summary["projects"])
        (staging / "candidate_inventory.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n",
                                                           encoding="utf-8")
        checksums = [f"{digest(path)}  {path.name}" for path in sorted(staging.iterdir()) if path.is_file()]
        (staging / "SHA256SUMS").write_text("\n".join(checksums) + "\n", encoding="utf-8")
        staging.rename(outdir)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return summary


def datetime_from_iso(value: str) -> str:
    from datetime import datetime
    result = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if result.tzinfo is None:
        raise ValueError("timestamp needs timezone")
    return result.isoformat()


def build(snapshot: Path, outdir: Path, config_path: Path = DEFAULT_CONFIG) -> dict:
    config = validate_config(config_path)
    # prepare() validates the complete archived raw snapshot and reuses the
    # repository's versioned AST normalization, policy, and exact joins.
    observations, _, ast_summary = ast_dataset.prepare(snapshot)
    snapshot_manifest_path = snapshot / "manifest.json"
    snapshot_manifest = json.loads(snapshot_manifest_path.read_text(encoding="utf-8"))
    source_entries = snapshot_manifest["sources"]
    record_ids = [row.get("record_id", "") for row in observations]
    source_ids = [row.get("source_record_id", "") for row in observations if not missing(row.get("source_record_id"))]
    if len(record_ids) != len(set(record_ids)):
        raise ValueError("duplicate normalized AST record_id")
    if len(source_ids) != len(set(source_ids)):
        raise ValueError("duplicate source AST record id")

    ecoli_rows = [row for row in observations if is_ecoli(row.get("organism", ""))]
    non_ecoli = [row for row in observations if not is_ecoli(row.get("organism", ""))]
    by_isolate = groups_by([r for r in ecoli_rows if not missing(r.get("isolate_id"))],
                           lambda r: r["isolate_id"])
    candidates: dict[str, dict] = {}
    for isolate_id, rows in sorted(by_isolate.items()):
            candidates[isolate_id], _ = make_candidate(config["cohort_id"], isolate_id, rows, config)

    endpoint_groups = groups_by([r for r in ecoli_rows if not missing(r.get("isolate_id"))
                                 and not missing(r.get("antimicrobial"))],
                                lambda r: (r["isolate_id"], r["antimicrobial"]))
    endpoints = [derive_endpoint(config["cohort_id"], isolate, drug, rows, candidates.get(isolate), config)
                 for (isolate, drug), rows in sorted(endpoint_groups.items())]
    exclusions: list[dict] = []
    for row in non_ecoli:
        exclusions.append({"cohort_id": config["cohort_id"], "analysis": "cohort", "isolate_id": row.get("isolate_id", NA),
                           "pdt_accession": row.get("isolate_id", NA), "antimicrobial": row.get("antimicrobial", NA),
                           "reason": "organism_mismatch", "source_record_ids_json": canonical([row.get("source_record_id", NA)])})
    for row in ecoli_rows:
        if missing(row.get("isolate_id")):
            exclusions.append({"cohort_id": config["cohort_id"], "analysis": "cohort", "isolate_id": NA,
                               "pdt_accession": NA, "antimicrobial": row.get("antimicrobial", NA),
                               "reason": "missing_isolate_identifier",
                               "source_record_ids_json": canonical([row.get("source_record_id", NA)])})
        elif missing(row.get("antimicrobial")):
            exclusions.append({"cohort_id": config["cohort_id"], "analysis": "endpoint", "isolate_id": row["isolate_id"],
                               "pdt_accession": row["isolate_id"], "antimicrobial": NA,
                               "reason": "missing_antimicrobial",
                               "source_record_ids_json": canonical([row.get("source_record_id", NA)])})
    for endpoint in endpoints:
        append_exclusions(endpoint, exclusions, config)

    endpoint_candidates = {r["isolate_id"] for r in endpoints if r["eligible_overall"] == "true"}
    eligible_isolates = [candidates[key] for key in sorted(endpoint_candidates)]
    antimicrobial_rows = summarize_antimicrobials(endpoints, config)
    country_rows = summarize_groups(endpoints, "country", COUNTRY_FIELDS)
    year_rows = summarize_groups(endpoints, "collection_year", YEAR_FIELDS)
    source_rows = source_summary(endpoints)
    missing_rows = missingness(ecoli_rows, list(candidates.values()))
    conflict_rows = [{
        "isolate_id": row["isolate_id"], "pdt_accession": row["pdt_accession"],
        "antimicrobial": row["antimicrobial"], "phenotype_endpoint": row["phenotype_endpoint"],
        "conflict_status": row["conflict_status"], "n_ast_observations": row["n_ast_observations"],
        "phenotypes_json": canonical(sorted({x.get("phenotype", NA) for x in endpoint_groups[
            (row["isolate_id"], row["antimicrobial"])]})),
        "source_record_ids_json": row["source_record_ids_json"],
    } for row in endpoints if row["conflict_status"] != "none"]

    assembly_isolates = {key for key, row in candidates.items()
                         if valid_assembly(row["assembly_accession"], config)}
    flow = {
        "cohort_id": config["cohort_id"],
        "units": {
            "raw_ast_observations": {"unit": "AST observations", "count": len(observations)},
            "ecoli_ast_observations": {"unit": "AST observations", "count": len(ecoli_rows)},
            "organism_excluded_observations": {"unit": "AST observations", "count": len(non_ecoli)},
            "unique_ecoli_isolates": {"unit": "isolates", "count": len(candidates)},
            "isolates_with_genbank_assembly": {"unit": "isolates", "count": len(assembly_isolates)},
            "isolate_drug_endpoints": {"unit": "isolate-drug endpoints", "count": len(endpoints)},
            "binary_rs_endpoints": {"unit": "endpoints", "count": sum(r["eligible_binary_rs"] == "true" for r in endpoints)},
            "discordant_rs_endpoints": {"unit": "endpoints", "count": sum(r["conflict_status"] == "discordant_rs" for r in endpoints)},
            "eligible_overall_endpoints": {"unit": "endpoints", "count": sum(r["eligible_overall"] == "true" for r in endpoints)},
            "eligible_country_endpoints": {"unit": "endpoints", "count": sum(r["eligible_country_analysis"] == "true" for r in endpoints)},
            "eligible_temporal_endpoints": {"unit": "endpoints", "count": sum(r["eligible_temporal_analysis"] == "true" for r in endpoints)},
            "eligible_source_endpoints": {"unit": "endpoints", "count": sum(r["eligible_source_analysis"] == "true" for r in endpoints)},
        },
        "reconciliation": {
            "raw_equals_ecoli_plus_organism_excluded": len(observations) == len(ecoli_rows) + len(non_ecoli),
            "endpoint_groups_equal_emitted_endpoint_rows": len(endpoints) == len(endpoint_groups),
        },
    }
    if not all(flow["reconciliation"].values()):
        raise ValueError("cohort flow counts do not reconcile")

    source_retrieved = sorted(entry["retrieved_at"] for entry in source_entries)
    generated_at = source_retrieved[-1] if source_retrieved else NA
    code_head = "unknown"
    try:
        code_head = subprocess.run(["git", "rev-parse", "HEAD"], check=True, capture_output=True,
                                   text=True, cwd=Path(__file__).resolve().parent).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        pass
    manifest = {
        "cohort_id": config["cohort_id"], "cohort_version": config["config_version"],
        "builder_version": VERSION, "source_snapshot_id": snapshot.name,
        "source_snapshot_manifest_sha256": digest(snapshot_manifest_path),
        "source_files": [{"filename": e["filename"], "role": e["role"], "source_url": e["source_url"],
                          "retrieved_at": e["retrieved_at"], "sha256": e["sha256"],
                          "parser_version": e["parser_version"], "schema_version": e["schema_version"]}
                         for e in source_entries],
        "configuration_sha256": digest(config_path), "code_commit": code_head,
        "generated_at": generated_at,
        "generated_at_semantics": "latest source retrieval time; fixed to keep rebuild outputs deterministic",
        "candidate_isolate_count": len(candidates),
        "eligible_isolate_count": len(endpoint_candidates), "endpoint_count": len(endpoints),
        "eligible_binary_rs_endpoint_count": flow["units"]["binary_rs_endpoints"]["count"],
        "ast_parser_summary": ast_summary,
    }

    outdir = Path(outdir)
    if outdir.exists() or outdir.is_symlink():
        raise ValueError(f"output already exists: {outdir}")
    outdir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".ecoli-cohort-", dir=outdir.parent))
    try:
        write_tsv(staging / "ast_observations.tsv", OBSERVATION_FIELDS, observations)
        write_tsv(staging / "candidate_isolates.tsv", CANDIDATE_FIELDS,
                  [candidates[k] for k in sorted(candidates)])
        write_tsv(staging / "isolate_drug_endpoints.tsv", ENDPOINT_FIELDS, endpoints)
        write_tsv(staging / "eligible_isolates.tsv", ELIGIBLE_FIELDS, eligible_isolates)
        write_tsv(staging / "exclusions.tsv", EXCLUSION_FIELDS,
                  sorted(exclusions, key=lambda r: (r["analysis"], r["isolate_id"], r["antimicrobial"], r["reason"])))
        write_tsv(staging / "antimicrobial_summary.tsv", ANTIMICROBIAL_FIELDS, antimicrobial_rows)
        write_tsv(staging / "country_summary.tsv", COUNTRY_FIELDS, country_rows)
        write_tsv(staging / "year_summary.tsv", YEAR_FIELDS, year_rows)
        write_tsv(staging / "source_summary.tsv", SOURCE_FIELDS, source_rows)
        write_tsv(staging / "missingness_summary.tsv", MISSING_FIELDS, missing_rows)
        write_tsv(staging / "conflict_summary.tsv", CONFLICT_FIELDS, conflict_rows)
        (staging / "cohort_flow.json").write_text(json.dumps(flow, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        (staging / "cohort_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        checksums = []
        for path in sorted(staging.iterdir()):
            if path.is_file():
                checksums.append(f"{digest(path)}  {path.name}")
        (staging / "SHA256SUMS").write_text("\n".join(checksums) + "\n", encoding="utf-8")
        if outdir.exists() or outdir.is_symlink():
            raise ValueError(f"output appeared during build: {outdir}")
        staging.rename(outdir)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return flow


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    command = sub.add_parser("build", help="build deterministic tables from an archived AST snapshot")
    command.add_argument("--snapshot", type=Path, required=True,
                         help="raw snapshot previously validated by bin/ast_dataset.py")
    command.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    command.add_argument("--outdir", type=Path, required=True)
    inventory_command = sub.add_parser("inventory", help="summarize an archived isolate metadata candidate list")
    inventory_command.add_argument("--input", type=Path, required=True)
    inventory_command.add_argument("--source-url", required=True)
    inventory_command.add_argument("--query", required=True)
    inventory_command.add_argument("--retrieved-at", required=True)
    inventory_command.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        result = (build(args.snapshot, args.outdir, args.config) if args.command == "build" else
                  inventory(args.input, args.outdir, args.source_url, args.query, args.retrieved_at))
        print(json.dumps(result, indent=2, sort_keys=True))
    except (ValueError, KeyError, TypeError, OSError, csv.Error) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

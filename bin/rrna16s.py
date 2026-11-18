#!/usr/bin/env python3
"""Prepare, verify and summarise local 16S BLAST searches. Python standard library only."""
import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import sys

PREFIX = "16S_ribosomal_RNA"
FIELDS = ["sample_id", "query_id", "original_id", "query_length", "ambiguous_fraction",
          "status", "candidate_taxa", "candidate_taxids", "candidate_accessions",
          "best_identity", "best_query_coverage", "best_alignment_length", "best_bitscore",
          "best_evalue", "reported_subjects", "candidate_subjects", "hit_limit_reached",
          "database_version", "min_identity", "min_coverage", "min_length",
          "max_ambiguous", "score_window_percent", "max_hits", "interpretation"]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def database_files(directory):
    directory = Path(directory)
    indices = sorted(p for p in directory.glob(PREFIX + ".n*")
                     if p.is_file() and not p.name.endswith((".gz", ".md5")))
    # Volume names can contain .00, .01, ... before the index suffix.
    indices += sorted(p for p in directory.glob(PREFIX + ".[0-9]*.n*") if p.is_file())
    if not indices:
        raise ValueError(f"No decompressed {PREFIX} BLAST index files in {directory}")
    taxonomy = [directory / "taxdb.btd", directory / "taxdb.bti"]
    for path in taxonomy:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing taxonomy file: {path}. Download/decompress taxdb first.")
    return sorted(set(indices + taxonomy))


def package_database(directory, version):
    directory = Path(directory)
    output = directory / "database.json"
    if output.exists():
        raise ValueError("database.json already exists; prepare a new directory for a new snapshot")
    if not str(version).strip():
        raise ValueError("A nonempty database version/snapshot label is required")
    files = database_files(directory)
    metadata = dict(schema_version=1, prefix=PREFIX, version=str(version),
                    source="https://ftp.ncbi.nlm.nih.gov/blast/db/",
                    packaged_at_utc=datetime.now(timezone.utc).isoformat(),
                    files={p.name: {"sha256": digest(p), "bytes": p.stat().st_size} for p in files})
    output.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata


def verify_database(directory, output):
    directory = Path(directory)
    metadata = json.loads((directory / "database.json").read_text(encoding="utf-8"))
    if metadata.get("schema_version") != 1 or metadata.get("prefix") != PREFIX or not metadata.get("version"):
        raise ValueError("Invalid database.json; use this script's database command")
    actual = {p.name for p in database_files(directory)}
    if actual != set(metadata.get("files", {})):
        raise ValueError("Database file list differs from database.json; use an unchanged snapshot")
    for name, info in metadata["files"].items():
        if Path(name).name != name:
            raise ValueError("Database manifest must contain plain filenames")
        path = directory / name
        if path.stat().st_size != info["bytes"] or digest(path) != info["sha256"]:
            raise ValueError(f"Database checksum mismatch: {name}")
    shutil.copyfile(directory / "database.json", output)
    return metadata


def fasta_records(path):
    path = Path(path)
    with path.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    name, chunks = None, []
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                header = line[1:].strip()
                if not header:
                    raise ValueError("Empty FASTA header")
                name, chunks = header.split()[0], []
            else:
                if name is None:
                    raise ValueError("Expected nucleotide FASTA. Supply prepared 16S sequences, not FASTQ.")
                chunks.append("".join(line.split()).upper())
                if sum(map(len, chunks)) > 5000:
                    raise ValueError("A 16S query exceeds 5000 bases. Use the WGS FASTA route for assemblies.")
    if name is not None:
        yield name, "".join(chunks)


def prepare(source, fasta, mapping):
    records = list(fasta_records(source))
    if not records:
        raise ValueError("No FASTA sequences found")
    ids, rows = set(), []
    for index, (name, sequence) in enumerate(records, 1):
        if name in ids:
            raise ValueError(f"Duplicate query identifier: {name}")
        ids.add(name)
        if not sequence or set(sequence) - set("ACGTRYSWKMBDHVN"):
            raise ValueError(f"Empty or non-IUPAC DNA sequence: {name}")
        rows.append(dict(query_id=f"q{index:06d}", original_id=name, length=len(sequence),
                         ambiguous_fraction=sum(c not in "ACGT" for c in sequence) / len(sequence)))
    with Path(fasta).open("w", encoding="utf-8") as handle:
        for row, (_, sequence) in zip(rows, records):
            handle.write(f'>{row["query_id"]}\n')
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")
    write_tsv(mapping, rows, ["query_id", "original_id", "length", "ambiguous_fraction"])


def parse_hits(path, queries):
    hits = {q: [] for q in queries}
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t", 14)
            if len(fields) != 15:
                raise ValueError(f"BLAST row {line_number}: expected 15 columns, got {len(fields)}")
            qid, acc, identity, length, qs, qe, qlen, ss, se, slen, evalue, bitscore, taxids, names, title = fields
            if qid not in queries or int(qlen) != int(queries[qid]["length"]):
                raise ValueError(f"BLAST query does not match query_map.tsv: {qid}")
            record = dict(accession=acc, identity=float(identity), alignment_length=int(length),
                          coverage=100 * (abs(int(qe) - int(qs)) + 1) / int(qlen),
                          evalue=float(evalue), bitscore=float(bitscore), taxids=taxids,
                          names=names, title=title)
            if any(not math.isfinite(record[k]) for k in ("identity", "coverage", "evalue", "bitscore")):
                raise ValueError(f"Non-finite BLAST value on row {line_number}")
            hits[qid].append(record)
    return hits


def taxa(names):
    return {s.strip() for s in names.split(";") if s.strip().lower() not in {"", "n/a", "na", "unknown", "-"}}


def summarise(args):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.sample_id):
        raise ValueError("Invalid sample ID")
    records = read_tsv(args.map)
    queries = {r["query_id"]: r for r in records}
    if not queries or len(queries) != len(records):
        raise ValueError("Empty query map or duplicate query IDs")
    db = json.loads(Path(args.database).read_text(encoding="utf-8"))
    hits = parse_hits(args.hits, queries)
    rows = []
    for qid, query in queries.items():
        reported = hits[qid]
        n_reported = len({h["accession"] for h in reported})
        truncated = n_reported >= args.max_hits
        eligible = [h for h in reported if h["identity"] >= args.min_identity and h["coverage"] >= args.min_coverage]
        best = max(eligible, key=lambda h: (h["bitscore"], h["coverage"], h["identity"])) if eligible else None
        candidates = ([h for h in eligible if h["bitscore"] >= best["bitscore"] * (1 - args.score_window / 100)]
                      if best else [])
        names = sorted(set().union(*(taxa(h["names"]) for h in candidates)))
        if int(query["length"]) < args.min_length:
            status = "short_query"
        elif float(query["ambiguous_fraction"]) > args.max_ambiguous:
            status = "ambiguous_query"
        elif not reported:
            status = "no_hits"
        elif not candidates:
            status = "below_threshold"
        elif truncated:
            status = "hit_limit_reached"
        elif any(not taxa(h["names"]) for h in candidates):
            status = "taxonomy_unavailable"
        elif len(names) > 1:
            status = "ambiguous_match"
        else:
            status = "candidate_match"
        row = dict.fromkeys(FIELDS, "")
        row.update(sample_id=args.sample_id, query_id=qid, original_id=query["original_id"],
                   query_length=int(query["length"]), ambiguous_fraction=float(query["ambiguous_fraction"]),
                   status=status, candidate_taxa=";".join(names),
                   candidate_taxids=";".join(sorted({t for h in candidates for t in h["taxids"].split(";") if t.isdigit() and t != "0"})),
                   candidate_accessions=";".join(sorted({h["accession"] for h in candidates})),
                   reported_subjects=n_reported, candidate_subjects=len(candidates), hit_limit_reached=truncated,
                   database_version=db["version"], min_identity=args.min_identity,
                   min_coverage=args.min_coverage, min_length=args.min_length,
                   max_ambiguous=args.max_ambiguous, score_window_percent=args.score_window,
                   max_hits=args.max_hits,
                   interpretation="Sequence-similarity candidates within this database; not a confirmed species identification.")
        if best:
            row.update(best_identity=best["identity"], best_query_coverage=best["coverage"],
                       best_alignment_length=best["alignment_length"], best_bitscore=best["bitscore"],
                       best_evalue=best["evalue"])
        rows.append(row)
    write_outputs(args.output, rows)
    return rows


def write_outputs(prefix, rows):
    write_tsv(str(prefix) + ".tsv", rows, FIELDS)
    Path(str(prefix) + ".json").write_text(json.dumps(rows, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def merge(output):
    rows, seen = [], set()
    paths = sorted(Path.cwd().glob("*.16s.json"))
    if not paths:
        raise ValueError("No *.16s.json summaries staged for merging")
    for path in paths:
        for row in json.loads(path.read_text(encoding="utf-8")):
            key = (row["sample_id"], row["query_id"])
            if key in seen:
                raise ValueError(f"Duplicate sample/query in summaries: {key}")
            seen.add(key)
            rows.append(row)
    write_outputs(output, sorted(rows, key=lambda r: (r["sample_id"], r["query_id"])))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    db = sub.add_parser("database", help="Freeze/checksum an already downloaded NCBI BLAST database")
    db.add_argument("--db-dir", required=True)
    db.add_argument("--version", required=True)
    check = sub.add_parser("verify")
    check.add_argument("--db-dir", required=True)
    check.add_argument("--output", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("--input", required=True)
    prep.add_argument("--fasta", required=True)
    prep.add_argument("--map", required=True)
    summary = sub.add_parser("summarise")
    for name in ("sample-id", "hits", "map", "database", "output"):
        summary.add_argument("--" + name, required=True)
    summary.add_argument("--min-identity", type=float, default=97.0)
    summary.add_argument("--min-coverage", type=float, default=90.0)
    summary.add_argument("--min-length", type=int, default=200)
    summary.add_argument("--max-ambiguous", type=float, default=0.01)
    summary.add_argument("--score-window", type=float, default=1.0)
    summary.add_argument("--max-hits", type=int, default=100)
    collect = sub.add_parser("merge")
    collect.add_argument("--output", required=True)
    args = p.parse_args()
    try:
        if args.command == "database":
            package_database(args.db_dir, args.version)
        elif args.command == "verify":
            verify_database(args.db_dir, args.output)
        elif args.command == "prepare":
            prepare(args.input, args.fasta, args.map)
        elif args.command == "summarise":
            if not (0 < args.min_identity <= 100 and 0 < args.min_coverage <= 100
                    and 0 <= args.max_ambiguous <= 1 and 0 <= args.score_window < 100
                    and args.min_length > 0 and args.max_hits >= 2):
                raise ValueError("Invalid comparison thresholds")
            summarise(args)
        else:
            merge(args.output)
    except (OSError, ValueError, KeyError, EOFError) as exc:
        sys.exit(f"ERROR: {exc}")


if __name__ == "__main__":
    main()

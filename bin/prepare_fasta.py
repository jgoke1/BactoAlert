#!/usr/bin/env python3
"""Validate a genome assembly and write an uncompressed .fasta for GHRU tools."""

import argparse
import gzip
from pathlib import Path
import sys
import tempfile
import zlib


def prepare_fasta(source: Path, destination: Path) -> tuple[int, int]:
    if source.resolve() == destination.resolve():
        raise ValueError("input and output paths must be different")

    with source.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    alphabet = set("ACGTRYSWKMBDHVN")
    identifiers = set()
    current_id = None
    current_length = 0
    total_bases = 0
    temporary = None

    try:
        with opener(source, "rt", encoding="utf-8-sig") as reader, tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="\n",
            dir=destination.parent, prefix=destination.name + ".", suffix=".tmp",
            delete=False,
        ) as writer:
            temporary = Path(writer.name)
            for line_number, raw_line in enumerate(reader, start=1):
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_id is not None and current_length == 0:
                        raise ValueError(f"contig '{current_id}' has no sequence")
                    header = line[1:].strip()
                    if not header:
                        raise ValueError(f"line {line_number}: FASTA header has no contig identifier")
                    current_id = header.split()[0]
                    if current_id in identifiers:
                        raise ValueError(f"duplicate contig identifier '{current_id}'")
                    identifiers.add(current_id)
                    current_length = 0
                    writer.write(line + "\n")
                else:
                    if current_id is None:
                        raise ValueError(
                            f"line {line_number}: expected a FASTA header starting with '>'; "
                            "use --input_format fastq for reads"
                        )
                    sequence = "".join(line.split())
                    invalid = set(sequence.upper()) - alphabet
                    if invalid:
                        raise ValueError(
                            f"line {line_number}: non-IUPAC DNA character(s) {''.join(sorted(invalid))!r}"
                        )
                    current_length += len(sequence)
                    total_bases += len(sequence)
                    writer.write(sequence + "\n")

            if current_id is None:
                raise ValueError("the FASTA contains no sequences")
            if current_length == 0:
                raise ValueError(f"contig '{current_id}' has no sequence")

        temporary.replace(destination)
        return len(identifiers), total_bases
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_fasta", type=Path)
    parser.add_argument("output_fasta", type=Path)
    args = parser.parse_args()
    try:
        contigs, bases = prepare_fasta(args.input_fasta, args.output_fasta)
    except (OSError, EOFError, UnicodeError, ValueError, zlib.error) as exc:
        raise SystemExit(f"ERROR: {args.input_fasta}: {exc}") from exc
    print(f"Prepared {contigs} contig(s), {bases} bases: {args.output_fasta}", file=sys.stderr)


if __name__ == "__main__":
    main()

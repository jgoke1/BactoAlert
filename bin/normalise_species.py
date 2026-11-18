#!/usr/bin/env python3
"""Normalise one Speciator result and propose a downstream analysis route.

Input: GHRU species-name TXT, Speciator JSON with speciesName, or a
single-row GHRU TSV with speciesName and optional Sample_id.
Output: JSON compatible with the NORMALISE_SPECIES Nextflow module.
Uses the species name reported by Speciator to select a Kleborate preset or
the general characterisation route. Confidence is retained as evidence but
does not change routing.
"""
import argparse
import csv
import io
import json
import re
from pathlib import Path

KPSC = {
    'klebsiella pneumoniae',
    'klebsiella quasipneumoniae',
    'klebsiella quasipneumoniae subsp. quasipneumoniae',
    'klebsiella quasipneumoniae subsp. similipneumoniae',
    'klebsiella variicola',
    'klebsiella variicola subsp. variicola',
    'klebsiella variicola subsp. tropica',
    'klebsiella quasivariicola',
    'klebsiella africana',
}
KOSC = {
    'klebsiella oxytoca', 'klebsiella michiganensis',
    'klebsiella grimontii', 'klebsiella pasteurii',
    'klebsiella huaxiensis', 'klebsiella spallanzanii',
}
UNRESOLVED = {
    'unknown', 'unclassified', 'unidentified', 'uncultured',
    'bacterium', 'bacteria', 'sp', 'spp', 'complex', 'group',
    'mixed', 'ambiguous', 'null', 'none', 'na',
}


def read_speciator(path, fmt, sample_id):
    text = path.read_text(encoding='utf-8-sig')
    if fmt == 'auto':
        if path.suffix.lower() == '.json' or text.lstrip().startswith(('{', '[')):
            fmt = 'json'
        elif path.suffix.lower() == '.tsv' or '\t' in text:
            fmt = 'tsv'
        else:
            fmt = 'txt'
    if fmt == 'txt':
        return {'speciesName': text.strip()}, fmt
    if fmt == 'json':
        record = json.loads(text)
        if not isinstance(record, dict):
            raise ValueError('Expected one Speciator JSON object.')
    else:
        reader = csv.DictReader(io.StringIO(text), delimiter='\t')
        headers = reader.fieldnames or []
        if len(headers) != len(set(headers)):
            raise ValueError('Duplicate TSV column names.')
        records = list(reader)
        if len(records) != 1:
            raise ValueError('Expected exactly one sample row in the TSV.')
        record = records[0]
        if None in record:
            raise ValueError('TSV row has more fields than its header.')
    if 'speciesName' not in record:
        raise ValueError('Required Speciator field speciesName is missing.')
    recorded_id = record.get('Sample_id') or record.get('sample_id')
    if recorded_id and str(recorded_id) != sample_id:
        raise ValueError(f'Sample ID mismatch: expected {sample_id!r}, found {recorded_id!r}.')
    return record, fmt


def normalise(raw, confidence=None):
    raw = raw if isinstance(raw, str) else ''
    name = ' '.join(raw.split())
    key = name.casefold()
    result = dict(
        species_raw=raw, species_normalised=name, genus='',
        species_complex='', route='general',
        reason='Speciator result does not match a targeted Kleborate route; routed to general characterisation',
        identification_source='Speciator', confidence_assessed=False,
        mapping_version='1.2',
    )
    if len(raw.splitlines()) > 1:
        result['reason'] = 'Multiple lines in species-name value'
        return result
    # Retain Speciator confidence as provenance, but route by species name.
    confidence_value = str(confidence or '').strip().casefold()
    result['speciator_confidence'] = confidence_value or None
    result['confidence_assessed'] = bool(confidence_value)
    if key in KPSC or key in KOSC:
        is_kpsc = key in KPSC
        result.update(
            species_normalised=key[0].upper() + key[1:],
            genus='Klebsiella', species_complex='KpSC' if is_kpsc else 'KoSC',
            route='kpsc' if is_kpsc else 'kosc',
            reason='Exact match to configured species-complex mapping',
        )
        return result
    parts = name.split()
    if len(parts) >= 2:
        genus = parts[0].capitalize()
        epithet = parts[1].lower()

        if genus in {'Escherichia', 'Shigella'}:
            normalised_name = ' '.join([genus, epithet, *parts[2:]])
            result.update(
                species_normalised=normalised_name,
                genus=genus,
                route='escherichia',
                reason=(
                    'Speciator genus is Escherichia or Shigella; routed to '
                    'the Kleborate Escherichia preset'
                ),
            )
            return result
    elif parts and parts[0].capitalize() in {'Escherichia', 'Shigella'}:
        genus = parts[0].capitalize()
        result.update(
            genus=genus,
            route='escherichia',
            reason=(
                'Speciator genus is Escherichia or Shigella; routed to '
                'the Kleborate Escherichia preset'
            ),
        )
        return result

    match = re.fullmatch(r'([A-Za-z]+) ([A-Za-z][A-Za-z-]+)', name)
    if match:
        genus, epithet = match.group(1).capitalize(), match.group(2).lower()
        if genus.lower() not in UNRESOLVED and epithet not in UNRESOLVED:
            result.update(
                species_normalised=f'{genus} {epithet}', genus=genus,
                route='general',
                reason=(
                    'Klebsiella name outside configured complex mapping'
                    if genus == 'Klebsiella'
                    else 'Named non-Klebsiella identification; taxonomy not independently verified'
                ),
            )
    elif parts and parts[0].isalpha():
        # Keep genus-only or unresolved calls in the general route. In
        # particular, Escherichia and Shigella were handled above.
        result['genus'] = parts[0].capitalize()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--format', choices=['auto', 'txt', 'json', 'tsv'], default='auto')
    args = parser.parse_args()
    if args.input.resolve() == args.output.resolve():
        parser.error('Input and output must be different files.')
    try:
        record, fmt = read_speciator(args.input, args.format, args.sample_id)
        raw = record.get('speciesName')
        if raw is None:
            raw = ''
        if not isinstance(raw, str):
            raise ValueError('speciesName must be a string or null.')
        result = normalise(raw.strip(), record.get('confidence'))
        # Retain key Speciator evidence so downstream reporting can distinguish
        # a taxonomic call from the evidence supporting that call.
        for field in (
            'taxId', 'speciesId', 'genusId', 'genusName',
            'superkingdomName', 'referenceId', 'mashDistance', 'pValue',
            'matchingHashes', 'confidence', 'source',
        ):
            if field in record:
                result[f'speciator_{field[0].lower()}{field[1:]}'] = record[field]
        # A filesystem-safe folder name used by the Nextflow grouping step.
        # Keep the original scientific name in species_normalised for reports.
        slug_source = result.get('species_normalised', '').casefold()
        result['species_slug'] = re.sub(r'[^a-z0-9]+', '_', slug_source).strip('_')[:80] or 'unresolved'
        result.update(sample_id=args.sample_id, input_format=fmt)
        args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    except (OSError, ValueError, csv.Error) as exc:
        parser.exit(1, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()

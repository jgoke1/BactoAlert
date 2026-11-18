#!/usr/bin/env python3
"""Merge the per-sample taxonomy JSON files staged into the current task."""
import csv
import glob
import json
from pathlib import Path

files = sorted(glob.glob('*.taxonomy.json'))
if not files:
    raise SystemExit('ERROR: no *.taxonomy.json files were staged for merging')

records = [json.loads(Path(name).read_text(encoding='utf-8')) for name in files]
sample_ids = [str(row.get('sample_id', '')) for row in records]
if any(not sample_id for sample_id in sample_ids):
    raise SystemExit('ERROR: a taxonomy JSON is missing sample_id')
if len(sample_ids) != len(set(sample_ids)):
    raise SystemExit('ERROR: duplicate sample_id values in taxonomy JSON inputs')

Path('species_normalisation_summary.json').write_text(
    json.dumps(records, indent=2) + '\n', encoding='utf-8'
)

preferred = [
    'sample_id', 'species_raw', 'species_normalised', 'species_slug',
    'genus', 'species_complex', 'route', 'reason',
    'speciator_confidence', 'speciator_taxId', 'speciator_speciesId',
    'speciator_genusName', 'speciator_referenceId',
    'speciator_mashDistance', 'speciator_matchingHashes',
    'speciator_pValue', 'speciator_source',
]
fields = preferred + sorted({key for row in records for key in row if key not in preferred})
with Path('species_normalisation_summary.tsv').open(
    'w', newline='', encoding='utf-8'
) as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(records)

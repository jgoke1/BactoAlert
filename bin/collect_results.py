#!/usr/bin/env python3
"""Preserve source artifacts and provide lossless, long-form tables for later R work."""
import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


GROUP_ANALYSES = {'genome_similarity', 'ska_snp', 'snp_phylogeny'}
COHORT_ANALYSES = {'taxonomy', 'taxonomy_json', 'rrna16s_summary'}


def load_species(manifest):
    """Read authoritative normalisation records; never infer species from filenames."""
    samples, groups = {}, {}
    for item in manifest:
        src = Path(item['path'])
        if item['analysis'] == 'taxonomy_json' and src.is_file():
            rows = json.loads(src.read_text())
            if not isinstance(rows, list):
                raise ValueError('Taxonomy summary must be a list of sample records')
            for row in rows:
                sid = str(row['sample_id'])
                value = (row.get('species_normalised', ''), row.get('species_slug', ''))
                if sid in samples and samples[sid] != value:
                    raise ValueError('Conflicting taxonomy for ' + sid)
                samples[sid] = value
                if value[1]:
                    groups[value[1]] = value[0]
        # Standalone comparisons have members tables but no cohort taxonomy summary.
        if item['analysis'] in GROUP_ANALYSES and src.is_dir():
            for member_file in src.rglob('members.tsv'):
                with member_file.open(newline='') as handle:
                    rows = list(csv.DictReader(handle, delimiter='\t'))
                names = {r.get('species', '') for r in rows} - {''}
                if len(names) == 1:
                    group_id = str(item['sample_id']).split('/')[0]
                    name = next(iter(names))
                    if group_id in groups and groups[group_id] != name:
                        raise ValueError('Conflicting species for group ' + group_id)
                    groups[group_id] = name
    return samples, groups


def manifest_row(item, path, out, samples, groups):
    analysis = item['analysis']
    identifier = str(item['sample_id'])
    group_id, sample_id, species, slug = '', '', '', ''
    scope = 'sample'
    if analysis in GROUP_ANALYSES:
        scope = 'group'
        group_id = identifier.split('/')[0]
        slug, species = group_id, groups.get(group_id, '')
    elif analysis in COHORT_ANALYSES:
        scope = 'cohort'
    else:
        sample_id = identifier
        species, slug = samples.get(sample_id, ('', ''))
    status, reason, source_status = 'completed', '', ''
    # A skipped status describes that specific step, not every output in its group.
    if path.name == 'status.json':
        record = json.loads(path.read_text())
        source_status = str(record.get('status', ''))
        reason = str(record.get('reason', '') or '')
        if source_status in ('not_performed', 'skipped') or source_status.startswith('insufficient_'):
            status = 'skipped'
            reason = reason or source_status
    return dict(sample_id=sample_id, species=species, species_slug=slug,
                group_id=group_id, scope=scope, analysis=analysis,
                file_type=''.join(path.suffixes).lstrip('.') or 'unknown',
                relative_path=path.relative_to(out).as_posix(), status=status,
                reason=reason, source_status=source_status,
                bytes=path.stat().st_size, sha256=digest(path))


def collect(manifest, out):
    out.mkdir(parents=True, exist_ok=False)
    inventory, versions, statuses, report_manifest = [], [], [], []
    samples, groups = load_species(manifest)
    with (out / 'measurements.tsv').open('w', newline='') as handle:
        fields = ['analysis', 'sample_or_group_id', 'source_file', 'record', 'column_index', 'field', 'value']
        writer = csv.DictWriter(handle, fields, delimiter='\t')
        writer.writeheader()
        for item in sorted(manifest, key=lambda x: (x['analysis'], x['sample_id'], x['path'])):
            src = Path(item['path'])
            dest = out / 'artifacts' / item['analysis'] / item['sample_id'] / src.name
            if dest.exists():
                raise ValueError(f'Artifact collision: {dest}')
            dest.parent.mkdir(parents=True, exist_ok=True)
            if src.is_dir():
                # Keep per-isolate BAM/VCF intermediates in snp_phylogeny only.
                ignore = shutil.ignore_patterns('calls') if item['analysis'] == 'snp_phylogeny' else None
                shutil.copytree(src, dest, ignore=ignore)
                files = sorted(x for x in dest.rglob('*') if x.is_file())
            else:
                shutil.copy2(src, dest)
                files = [dest]
            for path in files:
                rel = str(path.relative_to(out))
                base = dict(analysis=item['analysis'], sample_or_group_id=item['sample_id'], source_file=rel)
                entry = manifest_row(item, path, out, samples, groups)
                report_manifest.append(entry)
                inventory.append(dict(base, bytes=entry['bytes'], sha256=entry['sha256']))
                if 'version' in path.name.lower() and path.suffix == '.txt':
                    versions.append(dict(base, version_text=path.read_text().strip()))
                if path.name == 'status.json':
                    statuses.append(dict(base, **json.loads(path.read_text())))
                # Tool-native wide tables stay intact; long form preserves fields without
                # guessing clinical meaning or conflating missing values with negatives.
                if path.suffix == '.tsv' or path.name.endswith(('_kleborate.txt', '_hAMRonization.txt')):
                    with path.open(newline='') as source:
                        reader = csv.reader(source, delimiter='\t')
                        headers = next(reader, [])
                        for index, row in enumerate(reader, 1):
                            if not row:
                                continue
                            if len(row) > len(headers):
                                raise ValueError(f'More values than header columns in {path}, record {index}')
                            row += [''] * (len(headers) - len(row))
                            # Repeated/blank headers are retained verbatim. Position is
                            # the unique column identifier; dict parsing would lose data.
                            for column, (key, value) in enumerate(zip(headers, row), 1):
                                writer.writerow(dict(base, record=index,
                                    column_index=column, field=key, value=value))
                elif item['analysis'] == 'mlst' and path.suffix == '.txt':
                    for index, line in enumerate(path.read_text().splitlines(), 1):
                        values = line.split('\t')
                        if len(values) < 3:
                            raise ValueError(f'Invalid MLST row in {path}')
                        row = dict(input_file=values[0], scheme=values[1], sequence_type=values[2], alleles=';'.join(values[3:]))
                        for column, (key, value) in enumerate(row.items(), 1):
                            writer.writerow(dict(base, record=index, column_index=column, field=key, value=value))
    for filename, rows, fields in [
        ('manifest.tsv', report_manifest, ['sample_id','species','species_slug','group_id','scope','analysis','file_type','relative_path','status','reason','source_status','bytes','sha256']),
        ('artifact_manifest.tsv', inventory, ['analysis','sample_or_group_id','source_file','bytes','sha256']),
        ('versions.tsv', versions, ['analysis','sample_or_group_id','source_file','version_text'])
    ]:
        with (out / filename).open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fields, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
    (out / 'analysis_status.json').write_text(json.dumps(statuses, indent=2) + '\n')
    (out / 'collection.json').write_text(json.dumps(dict(schema_version='1.2', artifacts=len(inventory),
        analyses=sorted({x['analysis'] for x in manifest}), missing_semantics='Absent is not negative or susceptible'), indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    collect(json.loads(args.manifest.read_text()), args.output)


if __name__ == '__main__':
    main()


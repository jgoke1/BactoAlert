"""Build the offline report's input layout from this run's Nextflow emissions."""
import csv
import json
import re
import shutil
import sys
from pathlib import Path


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, records, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, list(dict.fromkeys(fields)), delimiter="\t")
        w.writeheader()
        w.writerows(records)


def safe(value):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError("Unsafe identifier: " + value)
    return value


def place(src, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        raise ValueError("Duplicate report input: " + str(dest))
    shutil.copyfile(src, dest)


def prepare(records, out):
    out.mkdir()
    tax = {}
    native = {}
    ix = {}
    for kind, sid, src in records:
        if kind == 'taxonomy':
            place(src, out/'species_normalisation/summary/species_normalisation_summary.tsv')
            for r in rows(src):
                sid = safe(r['sample_id'])
                safe(r['species_slug'])
                if sid in tax:
                    raise ValueError('Duplicate taxonomy for ' + sid)
                tax[sid] = r
    for kind, sid, src in records:
        safe(sid)
        if kind == 'checkm':
            place(src, out/'checkm_summary'/src.name)
        elif kind == '16s':
            place(src, out/'rrna16s'/src.name)
        elif 'versions' in kind:
            place(src, out/'versions'/kind/(sid+'.versions.txt'))
        elif kind == 'kleborate':
            folder = 'Escherichia_output' if tax[sid]['route']=='escherichia' else 'Klebsiella_output'
            name = sid+'_kleborate.txt'
            place(src/name, out/'kleborate'/folder/name)
        elif kind in ('mlst','amrfinder','virulencefinder','plasmidfinder'):
            if (kind,sid) in native: raise ValueError('Duplicate result for '+sid+' '+kind)
            native[kind,sid] = src
        elif kind in ('genome_similarity','ska_snp'):
            method = 'mash' if kind=='genome_similarity' else 'ska2'
            dest = out/kind/sid/method/src.name
            meta = json.loads((src/'status.json').read_text()) if (src/'status.json').exists() else {}
            for analysis, name in ([('mash','matrix.tsv'),('mash_tree','tree.nwk')] if method=='mash' else [('snp','matrix.tsv')]):
                path = dest/name
                if (src/name).is_file():
                    place(src/name,path)
                    status,reason='valid',''
                else:
                    status,reason='unavailable',meta.get('reason') or 'No '+name+' emitted'
                ix.setdefault(analysis,[]).append(dict(species_id=sid,relative_path=str(path.relative_to(out)) if status=='valid' else '',status=status,reason=reason))
    prep=out/'report_tables'
    prep.mkdir(exist_ok=True)
    for a, entries in ix.items():
        write(prep/(a+'_index.tsv'), entries, entries[0])
    for sp in sorted({r['species_slug'] for r in tax.values() if r['route']=='general'}):
        ids=sorted(i for i,r in tax.items() if r['species_slug']==sp and r['route']=='general')
        for kind,a,col in [('mlst','mlst',None),('amrfinder','amr','Element symbol'),('virulencefinder','virulence','gene'),('plasmidfinder','plasmid','plasmid')]:
            details,availability,hits=[],[],[]
            for sid in ids:
                p=native.get((kind,sid))
                status='missing' if p is None else 'readable'
                reason=''
                if p is not None:
                    try:
                        if a=='mlst':
                            with p.open(newline='') as h: lines=list(csv.reader(h,delimiter='\t'))
                            if len(lines)!=1 or len(lines[0])<3: raise ValueError('Invalid MLST row')
                            raw=[{'scheme':lines[0][1],'ST':lines[0][2],'alleles':'; '.join(lines[0][3:])}]
                        else:
                            with p.open(newline='',encoding='utf-8-sig') as h:
                                reader=csv.DictReader(h,delimiter='\t')
                                if col not in (reader.fieldnames or []): raise ValueError('Missing '+col+' column')
                                raw=list(reader)
                                if any(None in r or None in r.values() for r in raw): raise ValueError('Malformed TSV row')
                        for r in raw:
                            details.append(dict(r,sample_id=sid))
                            if col and r[col].strip().lower() not in ('','na','n/a','-'):
                                hits.append(dict(sample_id=sid,marker=r[col].strip()))
                    except (UnicodeError,ValueError,csv.Error) as e:
                        status,reason='unreadable',str(e)
                availability.append(dict(sample_id=sid,status=status,reason=reason))
            base=prep/sp
            write(base/(a+'_availability.tsv'),availability,['sample_id','status','reason'])
            write(base/(a+'_details.tsv'),details,['sample_id']+[k for r in details for k in r])
            write(base/(a+'_detections.tsv'),hits,['sample_id','marker'])


if __name__=='__main__':
    prepare([(k,s,Path(p)) for k,s,p in json.loads(Path(sys.argv[1]).read_text())],Path(sys.argv[2]))

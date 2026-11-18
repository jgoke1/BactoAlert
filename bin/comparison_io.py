#!/usr/bin/env python3
"""Metadata and alignment handling only. Mash, SKA, SNP counting and trees use external tools."""
import argparse
import csv
import hashlib
import json
import re
import shutil
from pathlib import Path


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2) + '\n')


def tsv(path, fields, rows):
    with Path(path).open('w', newline='') as h:
        w=csv.DictWriter(h, fields, delimiter='\t'); w.writeheader(); w.writerows(rows)


def fasta(path):
    records={}; key=None
    with Path(path).open() as h:
        for line in h:
            if line.startswith('>'):
                key=line[1:].split()[0]
                if key in records: raise ValueError('Duplicate sequence ID: '+key)
                records[key]=[]
            elif line.strip():
                if key is None: raise ValueError('Missing FASTA header')
                records[key].append(line.strip().upper())
    return {k:''.join(v) for k,v in records.items()}


def write_fasta(path, records):
    with Path(path).open('w') as h:
        for k,v in records.items(): h.write(f'>{k}\n{v}\n')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()


def members(group):
    return json.loads((Path(group)/'members.json').read_text())


def prepare(manifest, out):
    rows=json.loads(Path(manifest).read_text()); out=Path(out)
    (out/'genomes').mkdir(parents=True)
    ids=[r['sample_id'] for r in rows]
    if len(set(ids))!=len(ids): raise ValueError('Duplicate sample IDs')
    result=[]
    for i,row in enumerate(sorted(rows,key=lambda r:r['sample_id'])):
        sid=row['sample_id']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',sid) or sid=='sample_id':
            raise ValueError('Unsafe/reserved sample ID: '+sid)
        alias=f'G{i+1:06d}'
        seqs=fasta(row['assembly'])
        if not seqs or not all(seqs.values()): raise ValueError('Empty assembly: '+sid)
        write_fasta(out/'genomes'/f'{alias}.fasta',seqs)
        result.append(dict(row,alias=alias,assembly_sha256=sha(row['assembly'])))
    dump(out/'members.json',result)
    tsv(out/'members.tsv',['sample_id','alias','species','input_type','assembly_sha256'],
        [{k:r.get(k,'') for k in ['sample_id','alias','species','input_type','assembly_sha256']} for r in result])
    (out/'ska_input.tsv').write_text(''.join(f"{r['alias']}\tgroup/genomes/{r['alias']}.fasta\n" for r in result))


def relabel_tree(text, mapping):
    return re.sub(r"(?<=[(,])\s*('?)(G[0-9]+)\1(?=[:),;])",lambda m:mapping[m[2]],text)


def mask_gff(path,length):
    mask=bytearray(length)
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith('#'): continue
        f=line.split('\t'); start,end=int(f[3]),int(f[4])
        if not 1<=start<=end<=length: raise ValueError('GFF coordinates out of alignment bounds')
        mask[start-1:end]=b'\1'*(end-start+1)
    return mask


def filter_alignment(source, group, out, method, min_fraction, gff=None):
    out=Path(out);out.mkdir()
    rows=members(group); aliases=[r['alias'] for r in rows]
    seqs=fasta(source)
    if method.startswith('snippy'): seqs.pop('Reference',None)
    if set(seqs)!=set(aliases): raise ValueError('Alignment sample IDs do not match group manifest')
    if len({len(s) for s in seqs.values()})!=1: raise ValueError('Unequal sequence lengths')
    seqs={k:re.sub('[^ACGT]','N',v) for k,v in seqs.items()}
    length=len(next(iter(seqs.values())))
    mask=mask_gff(gff,length) if gff else bytearray(length)
    sites=[i for i in range(length) if not mask[i] and all(seqs[k][i] in 'ACGT' for k in aliases)]
    core={k:''.join(seqs[k][i] for i in sites) for k in aliases}
    variable=sum(len({core[k][i] for k in aliases})>1 for i in range(len(sites)))
    fraction=len(sites)/length if length else 0
    # For SKA this fraction is relative to the already selected split-kmer alignment,
    # NOT genome-wide coverage. Never equate it to reference coverage.
    eligible=len(aliases)>=2 and bool(sites) and (method=='ska2' or fraction>=min_fraction)
    stats=dict(method=method,n_samples=len(aliases),alignment_columns=length,callable_sites=len(sites),
        retained_fraction_of_input_alignment=fraction,variable_sites=variable,masked_columns=sum(mask),
        distance_eligible=eligible,unique_haplotypes=len(set(core.values())),
        status='ready' if eligible else 'insufficient_samples_or_callable_sites',
        denominator='shared_split_kmer_alignment_sites' if method=='ska2' else 'shared_reference_alignment_sites',
        recombination_filtered=bool(gff))
    stats['gubbins_eligible']=method=='snippy_unmasked' and eligible and len(aliases)>=4 and stats['unique_haplotypes']>=4 and variable>=10
    stats['tree_eligible']=eligible and len(aliases)>=4 and stats['unique_haplotypes']>=4 and variable>=10
    stats['tree_model']='GTR+G+ASC' if variable==len(sites) and sites else 'GTR+G'
    dump(out/'status.json',stats)
    write_fasta(out/'core.fasta',core)
    write_fasta(out/'full.fasta',seqs)
    tsv(out/'retained_columns.tsv',['input_alignment_column'],[dict(input_alignment_column=i+1) for i in sites])
    mapping={r['alias']:r['sample_id'] for r in rows}
    write_fasta(out/'core.sample_ids.fasta',{mapping[k]:v for k,v in core.items()})
    tsv(out/'callability.tsv',['sample_id','acgt_columns','input_alignment_columns'],
        [dict(sample_id=mapping[k],acgt_columns=sum(c in 'ACGT' for c in seqs[k]),input_alignment_columns=length) for k in aliases])
    shutil.copyfile(Path(group)/'members.tsv',out/'members.tsv')
    return stats


def read_matrix(path, mapping, mash=False):
    lines=[line.split('\t') if '\t' in line else line.split() for line in Path(path).read_text().splitlines() if line.strip()]
    if not lines: raise ValueError('Empty matrix')
    if len(lines[0])==1 and lines[0][0].isdigit():
        n=int(lines.pop(0)[0]); header=[r[0] for r in lines]
        if n!=len(header): raise ValueError('Invalid matrix count')
    else: header=lines.pop(0)[1:]
    if [r[0] for r in lines]!=header or set(header)!=set(mapping) or len(set(header))!=len(header):
        raise ValueError('Matrix labels do not match samples')
    values=[[float(v) if mash else int(v) for v in r[1:]] for r in lines]
    n=len(header)
    if any(len(r)!=n for r in values): raise ValueError('Matrix is not square')
    import math
    for i in range(n):
        for j in range(n):
            if not math.isfinite(values[i][j]) or values[i][j]<0 or abs(values[i][j]-values[j][i])>1e-8:
                raise ValueError('Invalid/asymmetric matrix')
        if values[i][i]!=0: raise ValueError('Nonzero diagonal')
    return [mapping[x] for x in header],values


def report(raw, group, out, method, info=None, tree=None):
    out=Path(out);out.mkdir(); mapping={r['alias']:r['sample_id'] for r in members(group)}
    ids,values=read_matrix(raw,mapping,method=='mash')
    key='mash_distance' if method=='mash' else 'snp_distance'
    stats=json.loads(Path(info).read_text()) if info else dict(method='mash',denominator='Mash sketch estimate')
    tsv(out/'matrix.tsv',['sample_id']+ids,[dict(sample_id=k,**dict(zip(ids,row))) for k,row in zip(ids,values)])
    edges=[dict(sample_a=ids[i],sample_b=ids[j],**{key:values[i][j]},comparable_sites=stats.get('callable_sites','NA'))
           for i in range(len(ids)) for j in range(i+1,len(ids))]
    tsv(out/'edges.tsv',['sample_a','sample_b',key,'comparable_sites'],edges)
    if tree: (out/'tree.nwk').write_text(relabel_tree(Path(tree).read_text(),mapping))
    stats.update(status='complete',outbreak_interpretation='Potential genomic relatedness only; epidemiological review required')
    dump(out/'status.json',stats)
    shutil.copyfile(Path(group)/'members.tsv',out/'members.tsv')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['prepare','filter','report','tree','status','reference'])
    for key in ['manifest','group','output','alignment','gff','matrix','info','tree','reference']:
        p.add_argument('--'+key)
    p.add_argument('--method',default='ska2');p.add_argument('--min-fraction',type=float,default=.9)
    p.add_argument('--reason',default='fewer_than_two_samples')
    a=p.parse_args()
    if not 0<a.min_fraction<=1: p.error('min-fraction must be in (0,1]')
    if a.mode=='prepare': prepare(a.manifest,a.output)
    elif a.mode=='filter': filter_alignment(a.alignment,a.group,a.output,a.method,a.min_fraction,a.gff)
    elif a.mode=='report': report(a.matrix,a.group,a.output,a.method,a.info,a.tree)
    elif a.mode=='tree':
        out=Path(a.output);out.mkdir();mapping={r['alias']:r['sample_id'] for r in members(a.group)}
        (out/'tree.nwk').write_text(relabel_tree(Path(a.tree).read_text(),mapping))
        dump(out/'status.json',dict(method='IQ-TREE3',status='complete'))
    elif a.mode=='reference':
        seqs=fasta(a.reference)
        if not seqs or not all(seqs.values()): raise ValueError('Reference must be nonempty FASTA')
        dump(a.output,dict(reference_sha256=sha(a.reference),reference_length=sum(map(len,seqs.values()))))
    else:
        out=Path(a.output);out.mkdir();dump(out/'status.json',dict(method=a.method,status='not_performed',reason=a.reason))
        shutil.copyfile(Path(a.group)/'members.tsv',out/'members.tsv')


if __name__=='__main__': main()

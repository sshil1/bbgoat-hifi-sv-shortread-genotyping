#!/usr/bin/env python3
"""Safely prepare a small, balanced sites-only candidate VCF from an existing HiFi SV VCF.

This is NOT an SV genotyper. No BAM, GATK, or existing VCF is changed.
Input should be the validated, already-created 37,113-call ARS1.2 autosomal SV VCF.

Example:
 python3 02_select_sv_pilot.py --input-vcf /actual/path/to/BBG_autosomal_SV.vcf.gz \
      --output-dir /path/to/NEW_SMALL_OUTPUT_DIRECTORY --per-bin 12 --seed 20261117

Output: pilot_candidates.vcf, pilot_candidates.tsv, selection_summary.txt
Do not pass a production GenomicsDB directory as output.
"""
from __future__ import annotations
import argparse, gzip, hashlib, random, sys
from collections import Counter, defaultdict
from pathlib import Path

BINS=((50,99,"50-99"),(100,499,"100-499"),(500,4999,"500-4999"),(5000,99999,"5000-99999"))

def op(path):
    return gzip.open(path,"rt",errors="replace") if str(path).endswith((".gz",".bgz")) else open(path,"rt",errors="replace")

def tags(s):
    out={}
    for x in s.split(';'):
        if '=' in x:
            k,v=x.split('=',1);out[k]=v
    return out

def genotype(sample, fmt):
    keys=fmt.split(':'); arr=sample.split(':')
    if 'GT' not in keys: return '.'
    ix=keys.index('GT')
    return arr[ix] if ix<len(arr) else '.'

def get_len(st,alt,ref,info):
    try:
        if info.get('SVLEN'):
            x=info['SVLEN'].split(',')[0]
            return abs(int(x))
    except ValueError: pass
    if not (alt.startswith('<') and alt.endswith('>')):
        return abs(len(alt)-len(ref))
    if info.get('END') and alt=='<DEL>':
        try:return abs(int(info['END'])-int(st))
        except ValueError:pass
    return 0

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--input-vcf',required=True,type=Path)
    a.add_argument('--output-dir',required=True,type=Path)
    a.add_argument('--per-bin',type=int,default=12)
    a.add_argument('--max-per-chrom',type=int,default=8)
    a.add_argument('--min-breakpoint-distance',type=int,default=2000)
    a.add_argument('--seed',type=int,default=20261117)
    args=a.parse_args()
    if args.per_bin<1 or args.max_per_chrom<1: a.error('per-bin and max-per-chrom must be positive')
    if not args.input_vcf.is_file(): a.error('Input VCF does not exist')
    if args.input_vcf.resolve()==args.output_dir.resolve(): a.error('Output must be a directory, not the VCF')
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        a.error('Output directory is non-empty. Use a fresh directory; do not overwrite outputs.')
    rng=random.Random(args.seed)
    entries=[]; metadata=[]; samples=[]; counts=Counter()
    filehash=hashlib.sha256()
    with op(args.input_vcf) as f:
        for rownum,line in enumerate(f,1):
            if line.startswith('##'):
                metadata.append(line); continue
            if line.startswith('#CHROM'):
                samples=line.rstrip('\n').split('\t')[9:];continue
            if not line.strip() or line.startswith('#'):continue
            counts['vcf_records']+=1
            c=line.rstrip('\n').split('\t')
            if len(c)<8:counts['malformed']+=1;continue
            chrom,pos,vid,ref,alt,qual,flt,inforaw=c[:8]
            if not chrom.startswith('NC_0308') or ',' in alt or '[' in alt or ']' in alt:
                counts['non_autosome_multiallelic_or_BND']+=1;continue
            try: posi=int(pos)
            except ValueError:counts['invalid_position']+=1;continue
            info=tags(inforaw)
            svtype=info.get('SVTYPE',('INS' if alt=='<INS>' else 'DEL' if alt=='<DEL>' else ''))
            if svtype not in ('DEL','INS'):
                counts['not_del_ins']+=1;continue
            if flt not in ('PASS','.'):
                counts['not_pass']+=1;continue
            svlen=get_len(posi,alt,ref,info)
            match=[b for b in BINS if b[0]<=svlen<=b[1]]
            if not match: counts['length_out_of_range']+=1;continue
            if svtype=='INS' and alt.startswith('<'):
                seq=next((info[k] for k in ('SEQ','INSSEQ','SVINSSEQ') if info.get(k)),None)
                if not seq:
                    counts['symbolic_insertion_missing_sequence']+=1;continue
                if len(seq)!=svlen or any(base not in 'ACGTNacgtn' for base in seq):
                    counts['symbolic_insertion_bad_sequence_length_or_bases']+=1;continue
            if svtype=='DEL' and alt=='<DEL>' and 'END' not in info:
                counts['symbolic_deletion_missing_END']+=1;continue
            if svtype=='INS' and not alt.startswith('<') and len(alt)<=len(ref):
                counts['insertion_allele_malformed']+=1;continue
            if svtype=='DEL' and not alt.startswith('<') and len(ref)<=len(alt):
                counts['deletion_allele_malformed']+=1;continue
            gts=[]
            if len(c)>9:
                gts=[genotype(g,c[8]) for g in c[9:]]
            informative=sum(any(x not in ('0','.') for x in gt.replace('|','/').split('/')) for gt in gts)
            # A retained ALT in an SV VCF is acceptable even if original sample GT is absent;
            # such loci receive no genotype-based control designation and must be evaluated manually.
            label=match[0][2]
            try:
                end=int(info.get('END',posi+svlen if svtype=='DEL' else posi))
            except ValueError:
                counts['invalid_END']+=1;continue
            if end<posi:end=posi
            entries.append(dict(index=counts['vcf_records'],chrom=chrom,pos=posi,end=end,
                 svtype=svtype,size=svlen,group=f'{svtype}:{label}',
                 informative=informative,gts=gts,sample_names=samples,fields=c[:8],original_id=vid))
            counts['eligible']+=1
    if not entries:
        print('ERROR: No eligible candidates. Review input VCF type, contigs, and ALT representation.',file=sys.stderr)
        sys.exit(2)
    buckets=defaultdict(list)
    for x in entries:buckets[x['group']].append(x)
    chosen=[]; chrom_used=Counter(); loci=defaultdict(list)
    for stype in ('DEL','INS'):
        for _,__,bn in BINS:
            group=f'{stype}:{bn}'
            arr=buckets[group]
            # Prefer variants with HiFi sample genotype annotations; balance by random seed within tiers.
            rng.shuffle(arr)
            arr.sort(key=lambda x: (0 if x['informative']>0 else 1))
            n=0
            for x in arr:
                if n>=args.per_bin:break
                if chrom_used[x['chrom']]>=args.max_per_chrom:continue
                if any(not (x['end']+args.min_breakpoint_distance < s or x['pos']-args.min_breakpoint_distance>e)
                       for s,e in loci[x['chrom']]):continue
                chosen.append(x);chrom_used[x['chrom']]+=1
                loci[x['chrom']].append((x['pos'],x['end']));n+=1
    chosen.sort(key=lambda x:x['index'])
    if not chosen:print('ERROR: No variants selected; relax spacing or per-chrom cap.',file=sys.stderr);sys.exit(3)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    vcf=args.output_dir/'pilot_candidates.vcf'
    tsv=args.output_dir/'pilot_candidates.tsv'
    report=args.output_dir/'selection_summary.txt'
    if any(p.exists() for p in (vcf,tsv,report)):
        print('ERROR: Refusing overwrite',file=sys.stderr);sys.exit(4)
    with vcf.open('w') as fo:
        for m in metadata:fo.write(m)
        fo.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for i,x in enumerate(chosen,1):
            c=x['fields'].copy()
            if c[2]=='.':c[2]=f'BBGSV_PILOT_{i:04d}'
            fo.write('\t'.join(c)+'\n')
    with tsv.open('w') as fo:
        fo.write('pilot_id\toriginal_id\tchrom\tpos_1based\tend_1based\tsvtype\tsvlen_bp\tsize_bin\talt_genotypes_in_input\tinput_sample_GTs\n')
        for i,x in enumerate(chosen,1):
            gtcol=','.join(f'{name}:{gt}' for name,gt in zip(x['sample_names'],x['gts'])) or 'NOT_RECORDED'
            fo.write(f'BBGSV_PILOT_{i:04d}\t{x["original_id"]}\t{x["chrom"]}\t{x["pos"]}\t{x["end"]}\t{x["svtype"]}\t{x["size"]}\t{x["group"]}\t{x["informative"]}\t{gtcol}\n')
    bins=Counter(x['group'] for x in chosen)
    with report.open('w') as fo:
        fo.write('SV pilot candidate selection (NOT population genotyping)\n')
        fo.write(f'Input: {args.input_vcf}\nSeed: {args.seed}\nRequested: {args.per_bin} per DEL/INS size bin\n')
        fo.write(f'Selected: {len(chosen)} (requested maximum {8*args.per_bin})\n')
        fo.write(f'Input VCF header sample IDs: {", ".join(samples) if samples else "NONE"}\n\n')
        fo.write('Selected per bin:\n')
        for stype in ('DEL','INS'):
            for _,__,bn in BINS:fo.write(f'{stype}:{bn}: {bins[f"{stype}:{bn}"]}\n')
        fo.write('\nInput screening counts:\n')
        for k in sorted(counts):fo.write(f'{k}: {counts[k]}\n')
        fo.write('\nIMPORTANT: Check candidate reference and allele representation before Paragraph; symbolic INS requires sequence in INFO (typically SEQ).\n')
        fo.write('CNVR 1868 (173 kb duplication) is outside the 100 kb consensus SV discovery set and is not represented by this pilot.\n')
    print(f'WROTE {len(chosen)} candidates -> {args.output_dir}')
    print('This script did not run any SV genotyper and did not change existing analysis files.')
if __name__=='__main__':main()

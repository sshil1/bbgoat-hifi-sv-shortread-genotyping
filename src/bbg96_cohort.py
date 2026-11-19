#!/usr/bin/env python3
"""BBG known-SV 96-site Paragraph cohort workflow (21 new + existing 803m).

Non-destructive: never modifies input BAM, BAI, VCF, reference or GATK recovery.
Modes: preflight, run, summarize, status.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

HOME = Path.home()
PILOT = HOME / 'BBG_SV_extension_20261117' / 'pilot_96_20261117'
ROOT = HOME / 'BBG_SV_extension_20261117' / 'cohort96_20261118'
BAM_DIR = Path('/data/bam/dedup')
JOINT = Path('/data/gatk/corrected_joint_22_20261112/final_autosomes/BBG22_corrected_29autosomes.vcf.gz')
REF = Path('/data/db/reference/goat_ARS1.fna')
SOURCE = PILOT / 'pilot_candidates.paragraph_sites.vcf'
PARAGRAPH_ENV = HOME / 'conda_envs/bbg_paragraph'
RUNNER = PARAGRAPH_ENV / 'bin/multigrmpy.py'
PYTHON = PARAGRAPH_ENV / 'bin/python'
EXPECTED_COUNT = 96
EXPECTED_SAMPLES = 22
EXISTING_SAMPLE = '803m'
PREVIOUS = [
    PILOT / 'paragraph_72_under5kb_803m/genotypes.vcf.gz',
    PILOT / 'paragraph_6_largeSV_803m/genotypes.vcf.gz',
    PILOT / 'paragraph_18_remaining_largeSV_803m/genotypes.vcf.gz',
]
EXCLUSION_FLAGS = {'BBGSV_PILOT_0015','BBGSV_PILOT_0027','BBGSV_PILOT_0050','BBGSV_PILOT_0052'}
MIN_FREE_BYTES = 350 * 1024**3


def call(command, **kwargs):
    return subprocess.run([str(x) for x in command], check=True, text=True, **kwargs)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def require(cmd):
    if shutil.which(cmd) is None:
        raise RuntimeError(f'Required command not found: {cmd}')


def parse_source():
    if not SOURCE.is_file():
        raise RuntimeError(f'96-site source VCF missing: {SOURCE}')
    hdr = None
    rows = {}
    with SOURCE.open() as f:
        for line in f:
            if line.startswith('#CHROM'):
                hdr = line.rstrip('\n').split('\t')
            elif not line.startswith('#'):
                a = line.rstrip('\n').split('\t')
                if len(a) < 8:
                    raise RuntimeError('Malformed source sites VCF')
                if a[2] in rows:
                    raise RuntimeError('Duplicate SV ID: ' + a[2])
                rows[a[2]] = tuple(a[:5])
    if len(rows) != EXPECTED_COUNT or hdr is None:
        raise RuntimeError(f'Expected 96 unique source SVs; found {len(rows)}')
    if len(hdr) != 8:
        raise RuntimeError('Source VCF must be sites-only (8 columns); stop rather than changing genotypes')
    return rows


def list_cohort():
    if not JOINT.is_file():
        raise RuntimeError('22-sample joint VCF header unavailable: ' + str(JOINT))
    p = call(['bcftools','query','-l',JOINT], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ids = [x for x in p.stdout.splitlines() if x]
    if len(ids) != EXPECTED_SAMPLES or len(set(ids)) != len(ids):
        raise RuntimeError(f'Joint VCF sample set is not exactly 22 unique samples: {ids}')
    if EXISTING_SAMPLE not in ids:
        raise RuntimeError('803m is absent from the joint VCF sample set; review cohort membership')
    for x in ids:
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', x):
            raise RuntimeError('Unsafe sample ID in joint VCF: ' + repr(x))
    return ids


def sample_ids_from_bam(path):
    result = call(['samtools','view','-H',path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ids = set()
    sq = {}
    for line in result.stdout.splitlines():
        if line.startswith('@RG\t'):
            fields = line.split('\t')
            ids.update(x[3:] for x in fields if x.startswith('SM:'))
        if line.startswith('@SQ\t'):
            d = dict(x.split(':', 1) for x in line.split('\t')[1:] if ':' in x)
            if 'SN' in d and 'LN' in d:
                sq[d['SN']] = int(d['LN'])
    return ids, sq


def reference_contigs():
    fai = Path(str(REF) + '.fai')
    if not fai.exists():
        raise RuntimeError(f'Reference FAI missing: {fai}')
    contigs = {}
    for line in fai.read_text().splitlines():
        a = line.split('\t')
        contigs[a[0]] = int(a[1])
    return contigs


def find_index(bam):
    candidates = [Path(str(bam)+'.bai'), bam.with_suffix('.bai'), Path(str(bam)+'.csi'), bam.with_suffix('.csi')]
    return next((x for x in candidates if x.is_file() and x.stat().st_size > 0), None)


def atomic_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix='.partial_', suffix='.tsv', dir=path.parent)
    try:
        with os.fdopen(fd, 'w', newline='') as f:
            w = csv.DictWriter(f, delimiter='\t', fieldnames=fields, extrasaction='ignore')
            w.writeheader()
            w.writerows(rows)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def preflight():
    for name in ('samtools','bcftools','timeout','nice','ionice'):
        require(name)
    if not (PYTHON.is_file() and RUNNER.is_file() and REF.is_file()):
        raise RuntimeError('Paragraph interpreter/script or reference missing')
    expected = parse_source()
    svtype_counts = {'DEL':0, 'INS':0}
    with (PILOT/'pilot_candidates.tsv').open() as meta_f:
        for r in csv.DictReader(meta_f, delimiter='\t'):
            if r['pilot_id'] in expected and r['svtype'] in svtype_counts:
                svtype_counts[r['svtype']] += 1
    if svtype_counts != {'DEL':48, 'INS':48}:
        raise RuntimeError(f'Expected 48 DEL and 48 INS, found {svtype_counts}')
    cohort = list_cohort()
    ref = reference_contigs()
    for infile in PREVIOUS:
        if not infile.is_file():
            raise RuntimeError('Existing 803m VCF missing; no 803m rerun authorized: ' + str(infile))

    mapped, issues = {}, []
    for bam in sorted(BAM_DIR.glob('*.dedup.bam')):
        basename = bam.name.removesuffix('.dedup.bam')
        if basename not in cohort:
            continue
        if basename in mapped:
            raise RuntimeError('Duplicate sample BAM: ' + basename)
        ids, sq = sample_ids_from_bam(bam)
        if ids != {basename}:
            issues.append(f'{basename}: RG SM={sorted(ids)}; expected exactly {basename}')
        autosomes = [c for c in ref if re.fullmatch(r'NC_0308(?:0[8-9]|[1-2][0-9]|3[0-6])\.1', c)]
        absent = [c for c in autosomes if c not in sq]
        if absent:
            issues.append(f'{basename}: {len(absent)} ARS1 autosomes missing from BAM header')
        mismatch = [c for c in autosomes if c in sq and ref[c] != sq[c]]
        if mismatch:
            issues.append(f'{basename}: reference lengths mismatch: {mismatch[:3]}')
        idx = find_index(bam)
        if not idx:
            issues.append(f'{basename}: BAM index unavailable')
        else:
            # Read-only index and file integrity preflight; no large BAM traversal.
            try:
                call(['samtools','quickcheck','-v',bam],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            except subprocess.CalledProcessError as exc:
                issues.append(f'{basename}: BAM quickcheck/index preflight failed: {exc}')
                continue
            mapped[basename] = {'sample': basename, 'source_bam':str(bam), 'source_index':str(idx)}

    missing = sorted(set(cohort) - set(mapped))
    if missing:
        issues.append('Cohort samples without indexed BAM: ' + ', '.join(missing))
    ROOT.mkdir(parents=True, exist_ok=True)
    (ROOT/'logs').mkdir(exist_ok=True)
    (ROOT/'runs').mkdir(exist_ok=True)
    (ROOT/'inputs').mkdir(exist_ok=True)
    (ROOT/'manifests').mkdir(exist_ok=True)
    (ROOT/'bam_links').mkdir(exist_ok=True)
    if issues:
        (ROOT/'preflight_issues.txt').write_text('\n'.join(issues)+'\n')
        print('PREFLIGHT_FAIL:')
        print('\n'.join(issues))
        print('Diagnostic:', ROOT/'preflight_issues.txt')
        raise SystemExit(2)

    # Create links only; never copy genome-size BAMs or indexes.
    for sample in cohort:
        r = mapped[sample]
        link = ROOT/'bam_links'/f'{sample}.bam'
        index_suffix = '.csi' if r['source_index'].endswith('.csi') else '.bai'
        ilink = Path(str(link) + index_suffix)
        for target, source in ((link, Path(r['source_bam'])), (ilink, Path(r['source_index']))):
            if target.is_symlink():
                if target.resolve() != source.resolve():
                    raise RuntimeError(f'Unexpected symlink target: {target}')
            elif target.exists():
                raise RuntimeError(f'Refusing to overwrite non-symlink {target}')
            else:
                target.symlink_to(source)
        # Verify the standard symlinked BAM/index pair, not the original index naming.
        try:
            call(['samtools','idxstats',link],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(f'Linked BAM index failed for {sample}: {exc}') from exc
        r['linked_bam'] = str(link)
        r['linked_index'] = str(ilink)

    atomic_tsv(ROOT/'cohort_inventory.tsv',
               ['sample','source_bam','source_index','linked_bam','linked_index'],
               [mapped[s] for s in cohort])
    snap = {'cohort': cohort, 'source_vcf_sha256':sha256(SOURCE),
            'ref':str(REF), 'joint_header_source':str(JOINT), 'n_sites':len(expected)}
    snapfile = ROOT/'cohort_snapshot.json'
    if snapfile.exists():
        old = json.loads(snapfile.read_text())
        if old != snap:
            raise RuntimeError('Cohort/source configuration drift; stop before changing snapshot')
    else:
        snapfile.write_text(json.dumps(snap, indent=2)+'\n')
    print('PREFLIGHT_PASS')
    print('Samples:',len(cohort),'including existing 803m; new genotyping targets:',len(cohort)-1)
    print('Fixed SVs:',len(expected), 'DEL/INS expected 48/48')
    print('Existing 803m result files:',len(PREVIOUS))
    print('Inventory:',ROOT/'cohort_inventory.tsv')
    print('Output root:',ROOT)
    print('Note: no Paragraph run has started.')


def load_inventory():
    p = ROOT/'cohort_inventory.tsv'
    if not p.exists():
        raise RuntimeError('Run preflight first')
    with p.open() as f:
        rows = {r['sample']:r for r in csv.DictReader(f, delimiter='\t')}
    snapshot = json.loads((ROOT/'cohort_snapshot.json').read_text())
    if snapshot['source_vcf_sha256'] != sha256(SOURCE):
        raise RuntimeError('Source 96-site VCF changed since preflight')
    if set(snapshot['cohort']) != set(rows):
        raise RuntimeError('Cohort inventory/sample mismatch')
    for sample in rows:
        r = rows[sample]
        if not Path(r['linked_bam']).exists() or not Path(r['linked_index']).exists():
            raise RuntimeError('Broken BAM/index link for '+sample)
    return rows, snapshot['cohort']


def create_input(sample):
    path = ROOT/'inputs'/f'{sample}.96SV.vcf'
    if path.exists():
        # Existing output must contain exactly the expected 96 fixed alleles.
        records = read_vcf(path)
        if dict((x[2], tuple(x[:5])) for x in records) != parse_source():
            raise RuntimeError('Existing sample input changed: '+sample)
        return path
    temp = path.with_suffix('.vcf.partial')
    if temp.exists():
        raise RuntimeError('Stale partial input found: '+str(temp))
    with SOURCE.open() as src, temp.open('x') as dst:
        have_gt = False
        for line in src:
            if line.startswith('##FORMAT=<ID=GT,'):
                have_gt = True
                dst.write(line)
            elif line.startswith('#CHROM'):
                if not have_gt:
                    dst.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Unassigned pre-genotyping allele">\n')
                dst.write(line.rstrip('\n') + '\tFORMAT\t'+sample+'\n')
            elif line.startswith('#'):
                dst.write(line)
            else:
                dst.write(line.rstrip('\n') + '\tGT\t./.\n')
    os.replace(temp, path)
    return path


def read_vcf(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    rows = []
    with opener(path, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                rows.append(line.rstrip('\n').split('\t'))
    return rows


def depth_readlen(bam):
    contigs = reference_contigs()
    autosomes = [(name,length) for name,length in contigs.items()
                 if re.fullmatch(r'NC_0308(?:0[8-9]|[1-2][0-9]|3[0-6])\.1',name) and length>=1000000]
    if len(autosomes) != 29:
        raise RuntimeError(f'Expected 29 goat ARS1 autosomes; found {len(autosomes)}')
    selected = [autosomes[i] for i in (0,2,4,6,8,10,12,14,16,18,20,22,24,26,27,28)]
    depths = []
    lengths = []
    for j, (chrom, length) in enumerate(selected):
        mid = int(length * (0.25 if j%2==0 else 0.65))
        start,end = mid, mid+2999
        region = f'{chrom}:{start}-{end}'
        p = call(['samtools','depth','-a','-s','-q','20','-Q','20','-G','SUPPLEMENTARY',
                  '-r',region,bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        vals = [int(line.rsplit('\t',1)[1]) for line in p.stdout.splitlines()]
        if len(vals) != 3000:
            raise RuntimeError(f'Incomplete depth region {region}: {len(vals)}')
        depths.append(sum(vals)/3000)
        if len(lengths)<1000:
            v = call(['samtools','view','-F','2308','-q','20',bam,region],
                     stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            for line in v.stdout.splitlines():
                a = line.split('\t')
                if len(a)>10 and a[9] != '*':
                    lengths.append(len(a[9]))
                    if len(lengths)>=1000: break
    d = statistics.median(depths)
    readlen = statistics.mean(lengths) if lengths else 0
    if not (8 <= d <= 120) or not (70 <= readlen <= 300):
        raise RuntimeError(f'Implausible provisional depth/readlen ({d}, {readlen}); manual review required')
    return {'depth':round(d,2),'readlen':round(readlen,2),
            'depth_min':round(min(depths),2),'depth_max':round(max(depths),2),
            'num_depth_windows':len(depths), 'reads_sampled':len(lengths)}


def is_gatk_using_bam(sample):
    # Fail-safe skip sample actively being processed by GATK, with no process modification.
    try:
        ps = call(['ps','-eo','args'],stdout=subprocess.PIPE,stderr=subprocess.PIPE).stdout
    except Exception:
        return True
    needle = f'{sample}.dedup.bam'
    return any('HaplotypeCaller' in line and needle in line for line in ps.splitlines())


def run(limit, allow_gatk_concurrency):
    inventory, cohort = load_inventory()
    required = [RUNNER,PYTHON,REF]
    if not all(x.exists() for x in required):
        raise RuntimeError('Paragraph environment/reference missing')
    if limit < 1:
        raise RuntimeError('--limit must be at least 1')
    lockfile = ROOT/'RUNNING.lock'
    fd = None
    try:
        fd = os.open(lockfile, os.O_CREAT|os.O_EXCL|os.O_WRONLY, 0o600)
        os.write(fd,f'pid={os.getpid()} time={time.ctime()}\n'.encode())
        os.close(fd)
        fd = None
    except FileExistsError:
        raise RuntimeError(f'RUNNING.lock exists; do not start concurrent cohorts: {lockfile}')
    started = 0
    try:
        for sample in cohort:
            if sample == EXISTING_SAMPLE:
                print('REUSE_803m: three existing 72/6/18-SV VCFs; no rerun')
                continue
            out = ROOT/'runs'/sample
            log = ROOT/'logs'/f'{sample}.log'
            if (out/'genotypes.vcf.gz').exists():
                print('ALREADY_COMPLETE_OR_REVIEW:', sample, '(check summarize validation)')
                continue
            if out.exists() or log.exists():
                print('STOP_FOR_MANUAL_REVIEW:',sample,'existing incomplete output or log')
                continue
            if is_gatk_using_bam(sample):
                print('SKIP_ACTIVE_GATK_BAM:',sample)
                continue
            free = shutil.disk_usage('/data').free
            if free < MIN_FREE_BYTES:
                print('STOP_LOW_DATA_DISK_SPACE:',round(free/1024**3,2),'GiB')
                break
            homefree = shutil.disk_usage(ROOT).free
            if homefree < 5 * 1024**3:
                print('STOP_LOW_HOME_SPACE:',round(homefree/1024**3,2),'GiB')
                break
            if not allow_gatk_concurrency:
                ps = call(['ps','-eo','args'],stdout=subprocess.PIPE).stdout
                if any('HaplotypeCaller' in x for x in ps.splitlines()):
                    print('STOP: GATK recovery active. Use --allow-gatk-concurrency only after resource review; default protects GATK.')
                    break
            bam = Path(inventory[sample]['linked_bam'])
            print('DEPTH_CALIBRATION:',sample,flush=True)
            estimates = depth_readlen(bam)
            manifest = ROOT/'manifests'/f'{sample}.tsv'
            manifest_text = ('id\tpath\tdepth\tread length\n'
                             + f"{sample}\t{bam}\t{estimates['depth']:.2f}\t{estimates['readlen']:.2f}\n")
            if manifest.exists():
                if manifest.read_text() != manifest_text:
                    raise RuntimeError('Existing manifest does not match recalibration: '+sample)
            else:
                manifest.write_text(manifest_text)
            (ROOT/'manifests'/f'{sample}.calibration.json').write_text(json.dumps(estimates,indent=2)+'\n')
            vcf = create_input(sample)
            # Preflight quick ref validation; fail if REF mismatch.
            call(['bcftools','norm','-c','e','-f',REF,'-Ov','-o','/dev/null',vcf],
                 stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
            maxreads = min(2000, max(200, math.ceil(20*estimates['depth'])))
            command = ['timeout','--signal=TERM','--kill-after=15s','45m',
                       'nice','-n','15','ionice','-c3',PYTHON,RUNNER,
                       '-i',vcf,'-m',manifest,'-r',REF,'-o',out,'--threads','1','-M',str(maxreads)]
            print('START_PARAGRAPH:',sample,'depth=',estimates['depth'],
                  'readlen=',estimates['readlen'],'maxreads=',maxreads, flush=True)
            (ROOT/'manifests'/f'{sample}.run_command.json').write_text(
                json.dumps({'command':[str(x) for x in command],
                            'source_vcf_sha256':sha256(SOURCE),
                            'calibration':estimates, 'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())},indent=2)+'\n')
            with log.open('x') as f:
                p = subprocess.run([str(x) for x in command],stdout=f,stderr=subprocess.STDOUT)
            print('PARAGRAPH_EXIT_STATUS:',sample,p.returncode,flush=True)
            started += 1
            if p.returncode != 0:
                print('STOP_AFTER_FAILURE:',sample,'log=',log)
                break
            try:
                got = read_vcf(out/'genotypes.vcf.gz')
                expected = parse_source()
                check = {a[2]:tuple(a[:5]) for a in got}
                if len(got) != 96 or check != expected:
                    raise RuntimeError(f'Output allele/count mismatch: {sample}')
            except Exception as exc:
                print('STOP_OUTPUT_VALIDATION:',sample,exc)
                break
            print('SAMPLE_PASS:',sample,'96 records; output=',out)
            if started >= limit:
                break
        print('RUN_COMMAND_FINISHED: newly attempted',started)
    finally:
        if fd is not None: os.close(fd)
        lockfile.unlink(missing_ok=True)


def load_output_file(path, sample, expected):
    if not path.is_file():
        raise RuntimeError('Missing result VCF: '+str(path))
    rows={}
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:
        names = None
        for line in f:
            if line.startswith('#CHROM'):
                names = line.rstrip('\n').split('\t')[9:]
            elif not line.startswith('#'):
                a = line.rstrip('\n').split('\t')
                if names != [sample] or len(a) != 10:
                    raise RuntimeError(f'Sample column mismatch: {path}: {names}')
                pid = a[2]
                if pid not in expected or tuple(a[:5]) != expected[pid] or pid in rows:
                    raise RuntimeError('Output ID/allele mismatch: '+pid+' in '+str(path))
                fmt = dict(zip(a[8].split(':'),a[9].split(':')))
                rows[pid]={'gt':fmt.get('GT','.'),'ft':fmt.get('FT','.'),
                           'dp':fmt.get('DP','.'),'gq':fmt.get('GQ','.'),
                           'pl':fmt.get('PL','.'),'ad':fmt.get('AD','.')}
    return rows


def normalize_gt(gt):
    a=gt.replace('|','/').split('/')
    if len(a)!=2 or any(x not in ('0','1') for x in a): return None
    return '/'.join(sorted(a))


def summarize():
    inventory, cohort = load_inventory()
    expected = parse_source()
    pid_type={}
    with (PILOT/'pilot_candidates.tsv').open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['pilot_id'] in expected:
                pid_type[r['pilot_id']]={'type':r['svtype'],'length':r['svlen_bp']}
    if set(pid_type)!=set(expected):
        raise RuntimeError('Missing metadata for some SVs')
    try:
        pilot = {}
        with (PILOT/'BBG_96SV_final_integrated_audit.tsv').open() as f:
            for r in csv.DictReader(f, delimiter='\t'):
                pilot[r['pilot_id']]=r
    except FileNotFoundError:
        raise RuntimeError('Original 96-SV audit missing')
    if set(pilot)!=set(expected):
        raise RuntimeError('Original pilot audit IDs do not match panel')
    by_sample={}
    for sample in cohort:
        if sample == EXISTING_SAMPLE:
            merged = {}
            for prev in PREVIOUS:
                for pid,callrow in load_output_file(prev,sample,expected).items():
                    if pid in merged:
                        raise RuntimeError('Duplicate 803m SV call: '+pid)
                    merged[pid]=callrow
            if set(merged)!=set(expected):
                raise RuntimeError('Existing 803m outputs do not form complete 96-SV panel')
            by_sample[sample]=merged
        else:
            path = ROOT/'runs'/sample/'genotypes.vcf.gz'
            if path.exists():
                calls = load_output_file(path,sample,expected)
                if set(calls) != set(expected):
                    raise RuntimeError(f'Incomplete variant set for {sample}')
                by_sample[sample]=calls
    n = len(by_sample)
    print('SAMPLES_WITH_COMPLETE_96_CALLS:',n,'/',len(cohort))
    if not n: raise RuntimeError('No sample results')
    long=[]
    for sample in cohort:
        if sample not in by_sample: continue
        for pid in sorted(expected):
            c=by_sample[sample][pid]
            gt=normalize_gt(c['gt'])
            try: dp=int(c['dp'])
            except ValueError: dp=-1
            ok=(c['ft']=='PASS' and gt is not None)
            watch=(pilot[pid]['comparison']=='DISCORDANT'
                   or pid in EXCLUSION_FLAGS or dp<10)
            long.append({'sample':sample,'pilot_id':pid,'svtype':pid_type[pid]['type'],
                         'svlen_bp':pid_type[pid]['length'],'gt':c['gt'],'ft':c['ft'],
                         'dp':c['dp'],'gq':c['gq'],'ad':c['ad'],'pl':c['pl'],
                         'callable_ft_pass':int(ok),'review_flag':int(watch)})
    matrix=[]
    for pid in sorted(expected):
        r={'pilot_id':pid,'svtype':pid_type[pid]['type'],'svlen_bp':pid_type[pid]['length']}
        for sample in cohort:
            if sample in by_sample:
                c=by_sample[sample][pid]
                r[sample]=f"{c['gt']}|{c['ft']}|DP={c['dp']}"
            else:
                r[sample]='NOT_RUN'
        matrix.append(r)
    atomic_tsv(ROOT/'cohort_96SV_genotype_matrix.tsv',
               ['pilot_id','svtype','svlen_bp']+cohort,matrix)
    atomic_tsv(ROOT/'cohort_genotypes_long.tsv',
               ['sample','pilot_id','svtype','svlen_bp','gt','ft','dp','gq','ad','pl',
                'callable_ft_pass','review_flag'],long)
    sample_rows=[]
    for s in cohort:
        if s not in by_sample: continue
        z=[r for r in long if r['sample']==s]
        ok=sum(r['callable_ft_pass'] for r in z)
        dp_low_pass=sum(r['callable_ft_pass'] and r['dp'].isdigit() and int(r['dp'])<10 for r in z)
        sample_rows.append({'sample':s,'n_sites':96,'ft_pass_callable':ok,
                            'ft_pass_rate':f'{ok/96:.5f}',
                            'pass_dp_below10':dp_low_pass})
    atomic_tsv(ROOT/'cohort_sample_summary.tsv',list(sample_rows[0]),sample_rows)
    locus=[]
    for pid in sorted(expected):
        z=[r for r in long if r['pilot_id']==pid]
        good=[r for r in z if r['callable_ft_pass']]
        gts=[normalize_gt(r['gt']) for r in good]
        alt=sum(gt.count('1') for gt in gts)
        denom=2*len(gts)
        pilot_flags=pilot[pid]['comparison']=='DISCORDANT' or pid in EXCLUSION_FLAGS
        low=sum(1 for r in good if r['dp'].isdigit() and int(r['dp'])<10)
        primary_ready=(n==22 and len(good)>=20 and not pilot_flags and low==0)
        locus.append({'pilot_id':pid,'svtype':pid_type[pid]['type'],
                      'svlen_bp':pid_type[pid]['length'],
                      'samples_processed':n,'ft_pass':len(good),
                      'ft_pass_rate_processed':f'{len(good)/n:.5f}',
                      'n_00':gts.count('0/0'),'n_01':gts.count('0/1'),
                      'n_11':gts.count('1/1'),'alt_allele_count':alt,
                      'called_alleles':denom,
                      'descriptive_af':f'{alt/denom:.6f}' if denom else '',
                      'pilot_discordance_or_special_flag':int(pilot_flags),
                      'low_dp_lt10_pass_count':low,
                      'exploratory_primary_ready':int(primary_ready)})
    atomic_tsv(ROOT/'cohort_locus_summary.tsv',list(locus[0]),locus)
    print('GENOTYPE_MATRIX:', ROOT/'cohort_96SV_genotype_matrix.tsv')
    print('LONG_FORMAT:', ROOT/'cohort_genotypes_long.tsv')
    print('SAMPLE_SUMMARY:', ROOT/'cohort_sample_summary.tsv')
    print('LOCUS_SUMMARY:', ROOT/'cohort_locus_summary.tsv')
    print('COMPLETE_COHORT:', 'YES' if n==22 else 'NO: frequencies interim only')
    print('PRIMARY_READY_SITES:',sum(r['exploratory_primary_ready'] for r in locus))
    print('SUMMARY_PASS')


def status():
    print('===== BBG 96-SV COHORT STATUS =====')
    try:
        _, cohort=load_inventory()
    except Exception as exc:
        print('PREFLIGHT_NOT_READY:',exc)
        return
    for s in cohort:
        if s==EXISTING_SAMPLE:
            print(s,'REUSED_EXISTING')
        else:
            out=ROOT/'runs'/s/'genotypes.vcf.gz'
            log=ROOT/'logs'/f'{s}.log'
            print(s,'OUTPUT_PRESENT' if out.is_file() else
                  'INCOMPLETE_REVIEW' if log.exists() or (ROOT/'runs'/s).exists() else 'PENDING')
    print('DATA_FREE_GIB:',round(shutil.disk_usage('/data').free/1024**3,1))
    print('GATK_ACTIVE:','yes' if any('HaplotypeCaller' in x for x in call(['ps','-eo','args'],stdout=subprocess.PIPE).stdout.splitlines()) else 'no')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['preflight','run','summarize','status'])
    parser.add_argument('--limit',type=int,default=1,
                        help='Max new sample Paragraph runs in this invocation (default 1)')
    parser.add_argument('--allow-gatk-concurrency',action='store_true',
                        help='Explicitly permit low-priority Paragraph while GATK on *other* BAMs is active; skips a BAM active in GATK')
    a=parser.parse_args()
    try:
        if a.mode=='preflight':preflight()
        elif a.mode=='run':run(a.limit,a.allow_gatk_concurrency)
        elif a.mode=='summarize':summarize()
        elif a.mode=='status':status()
    except (OSError, RuntimeError, subprocess.CalledProcessError, AssertionError) as e:
        print('ERROR:',str(e),file=sys.stderr)
        sys.exit(2)

if __name__=='__main__':main()

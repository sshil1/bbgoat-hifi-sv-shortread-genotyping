#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="$HOME/BBG_SV_extension_20261117/cohort96_20261118"
CONTROLLER="$HOME/bbg96_cohort.py"
LONG="$ROOT/cohort_genotypes_long.tsv"

mkdir -p "$ROOT"

LOG="$ROOT/remaining_batch_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

trap 'rc=$?; echo "BBG_BATCH_EXIT_STATUS=$rc"' EXIT

exec 9>"$ROOT/.remaining_batch.lock"
if ! flock -n 9; then
    echo "STOP: Another BBG batch is running"
    exit 1
fi

SAMPLES=(
    584f
    7ps__f_
    814m
    836m
    849m
    855m
    9ks__f_
    Bk_6__m_
    Hir_13_m_
    Hir_5f
    Kh_5__m_
)

echo "========================================="
echo "BBG 96-SV REMAINING COHORT BATCH"
echo "========================================="
date
echo "LOG=$LOG"
echo "EXPECTED_NEW_SAMPLES=${#SAMPLES[@]}"

echo "===== INITIAL STATUS ====="
python3 "$CONTROLLER" status

check_bam() {
    local sample="$1"
    local bam="/data/bam/dedup/${sample}.dedup.bam"

    echo "BAM_PREFLIGHT: $sample"

    if [[ ! -s "$bam" ]]; then
        echo "STOP_MISSING_BAM: $bam"
        return 1
    fi

    samtools quickcheck -v "$bam"

    samtools view -H "$bam" |
        awk -F '\t' -v expected="$sample" '
        $1=="@RG" {
            n++
            found=0
            for(i=2;i<=NF;i++) {
                if(substr($i,1,3)=="SM:") {
                    found=1
                    if(substr($i,4)!=expected) bad=1
                }
            }
            if(!found) bad=1
        }
        END {
            if(n==0 || bad) {
                print "SAMPLE_IDENTITY_FAIL"
                exit 1
            }
        }'

    samtools idxstats "$bam" >/dev/null

    echo "BAM_PREFLIGHT_PASS: $sample"
}

resource_check() {
    local vm idle wait mem disk

    vm=$(vmstat 1 4)

    read -r idle wait < <(
        printf '%s\n' "$vm" |
        tail -n 3 |
        awk '{
            idle+=$15
            wait+=$16
            n++
        } END {
            printf "%.1f %.1f\n",idle/n,wait/n
        }'
    )

    mem=$(free -g | awk '/^Mem:/ {print $7}')

    disk=$(
        df -BG --output=avail /data |
        tail -1 | tr -dc '0-9'
    )

    echo "CPU_IDLE=$idle"
    echo "IOWAIT=$wait"
    echo "MEM_AVAILABLE_GIB=$mem"
    echo "DISK_AVAILABLE_GIB=$disk"

    if (( mem < 24 || disk < 400 )); then
        echo "RESOURCE_BLOCK: Memory/storage"
        return 1
    fi

    if ! awk -v i="$idle" -v w="$wait" \
        'BEGIN {exit !(i>=25 && w<=5)}'
    then
        echo "RESOURCE_BLOCK: CPU/I/O"
        return 1
    fi

    if pgrep -af '[b]wa-mem2 mem' |
        awk '{
            for(i=1;i<=NF;i++) {
                if($i=="-t" && $(i+1)>16)
                    bad=1
            }
        } END {exit bad ? 0 : 1}'
    then
        echo "RESOURCE_BLOCK: BWA exceeds 16 threads"
        return 1
    fi

    if pgrep -f \
        '[m]ultigrmpy.py|[g]rmpy|[b]bg96_cohort.py run' \
        >/dev/null
    then
        echo "RESOURCE_BLOCK: Paragraph already active"
        return 1
    fi

    echo "RESOURCE_GATE_PASS"
    return 0
}

qc_sample() {
    local sample="$1"

    python3 - "$LONG" "$sample" <<'PY'
import csv
import sys
from collections import Counter
from pathlib import Path

path = Path(sys.argv[1])
sample_name = sys.argv[2]

with path.open() as f:
    rows = [
        r for r in csv.DictReader(f, delimiter="\t")
        if r["sample"] == sample_name
    ]

assert len(rows) == 96, (
    f"{sample_name}: expected 96, got {len(rows)}"
)
assert len({r["pilot_id"] for r in rows}) == 96

types = Counter(r["svtype"] for r in rows)
assert types["DEL"] == 48 and types["INS"] == 48

passed = [
    r for r in rows
    if r["callable_ft_pass"] == "1"
]

low_dp = []
for r in passed:
    try:
        dp = int(r["dp"])
    except (ValueError, TypeError):
        low_dp.append((r["pilot_id"], "UNKNOWN"))
        continue
    if dp < 10:
        low_dp.append((r["pilot_id"], str(dp)))

print("SAMPLE_QC:", sample_name)
print("RECORDS:", len(rows))
print("FT_PASS:", len(passed))
print("FT_PASS_PERCENT:", round(100*len(passed)/96, 2))
print("DEL_PASS:", sum(r["svtype"]=="DEL" for r in passed))
print("INS_PASS:", sum(r["svtype"]=="INS" for r in passed))
print("LOW_DP_OR_UNKNOWN_PASS:", low_dp)
print("FILTERS:", dict(Counter(r["ft"] for r in rows)))

for suffix in ("0014","0015","0027","0050","0052","0063"):
    pid = "BBGSV_PILOT_" + suffix
    r = next(x for x in rows if x["pilot_id"] == pid)
    print(
        "REVIEW_LOCUS:",
        pid,
        "GT=" + r["gt"],
        "FT=" + r["ft"],
        "DP=" + r["dp"]
    )

print("SAMPLE_QC_PASS:", sample_name)
PY
}

for SAMPLE in "${SAMPLES[@]}"; do

    echo
    echo "========================================="
    echo "SAMPLE_CHECKPOINT: $SAMPLE"
    echo "========================================="
    date

    STATUS=$(python3 "$CONTROLLER" status)

    if grep -qx "$SAMPLE OUTPUT_PRESENT" <<< "$STATUS"; then
        echo "SKIP_ALREADY_COMPLETE: $SAMPLE"
        continue
    fi

    if ! grep -qx "$SAMPLE PENDING" <<< "$STATUS"; then
        echo "STOP_UNEXPECTED_STATUS: $SAMPLE"
        echo "$STATUS"
        exit 2
    fi

    check_bam "$SAMPLE"

    READY=0

    for ATTEMPT in 1 2 3 4 5; do
        echo "RESOURCE_ATTEMPT=$ATTEMPT SAMPLE=$SAMPLE"

        if resource_check; then
            READY=1
            break
        fi

        if (( ATTEMPT < 5 )); then
            echo "RESOURCE_RECHECK_IN_60_SECONDS"
            sleep 60
        fi
    done

    if (( READY != 1 )); then
        echo "BATCH_PAUSED_RESOURCE_GATE: $SAMPLE"
        echo "Restart this script to resume safely"
        exit 3
    fi

    echo "START_SINGLE_SAMPLE: $SAMPLE"

    if pgrep -f '[H]aplotypeCaller' >/dev/null; then
        nice -n 10 python3 "$CONTROLLER" run \
            --limit 1 \
            --allow-gatk-concurrency
    else
        nice -n 10 python3 "$CONTROLLER" run \
            --limit 1
    fi

    STATUS_AFTER=$(python3 "$CONTROLLER" status)

    if ! grep -qx "$SAMPLE OUTPUT_PRESENT" <<< "$STATUS_AFTER"; then
        echo "STOP_OUTPUT_NOT_CONFIRMED: $SAMPLE"
        echo "$STATUS_AFTER"
        exit 4
    fi

    python3 "$CONTROLLER" summarize
    qc_sample "$SAMPLE"

    echo "CHECKPOINT_COMPLETE: $SAMPLE"
    date

done

echo
echo "========================================="
echo "VERIFY COMPLETE 22-ANIMAL COHORT"
echo "========================================="

python3 "$CONTROLLER" summarize
python3 "$CONTROLLER" status

python3 - "$LONG" "$ROOT" <<'PY'
import csv
import sys
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime

path = Path(sys.argv[1])
outdir = Path(sys.argv[2])

with path.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

expected_samples = {
    "116f", "189f", "1buck_m_", "2183f",
    "246f", "2869f", "3237m", "369m",
    "4_Female", "4lb__f_", "584f", "7ps__f_",
    "803m", "814m", "836m", "849m", "855m",
    "9ks__f_", "Bk_6__m_", "Hir_13_m_",
    "Hir_5f", "Kh_5__m_"
}

sample_counts = Counter(r["sample"] for r in rows)

assert set(sample_counts) == expected_samples, (
    "Incomplete or unexpected sample set"
)
assert all(n == 96 for n in sample_counts.values())
assert len(rows) == 2112

pairs = [(r["sample"],r["pilot_id"]) for r in rows]
assert len(pairs) == len(set(pairs))

groups = defaultdict(list)

for r in rows:
    groups[r["pilot_id"]].append(r)

assert len(groups) == 96

def valid_gt(gt):
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2:
        return None
    if any(a not in ("0", "1") for a in alleles):
        return None
    return sum(int(a) for a in alleles)

def depth_at_least_10(r):
    try:
        return float(r["dp"]) >= 10
    except (ValueError, TypeError, KeyError):
        return False

results = []

for pid, records in sorted(groups.items()):
    assert len(records) == 22

    typ = {r["svtype"] for r in records}
    assert len(typ) == 1
    typ = typ.pop()

    ft_pass = [
        r for r in records
        if r["callable_ft_pass"] == "1"
    ]

    for r in ft_pass:
        assert valid_gt(r["gt"]) is not None, (
            f"Invalid FT-PASS genotype: {pid}, {r['sample']}"
        )

    dp_pass = [
        r for r in ft_pass
        if depth_at_least_10(r)
    ]

    n_ft = len(ft_pass)
    n_dp = len(dp_pass)

    primary = n_ft >= 20
    sensitivity = n_dp >= 20

    alt_ft = sum(valid_gt(r["gt"]) for r in ft_pass)
    alt_dp = sum(valid_gt(r["gt"]) for r in dp_pass)

    results.append({
        "pilot_id": pid,
        "svtype": typ,
        "animals_total": 22,
        "ft_pass_n": n_ft,
        "dp10_ft_pass_n": n_dp,
        "primary_ft20_eligible": int(primary),
        "sensitivity_dp10_20_eligible": int(sensitivity),
        "primary_alt_allele_frequency": (
            f"{alt_ft/(2*n_ft):.6f}"
            if primary else "NA"
        ),
        "sensitivity_alt_allele_frequency": (
            f"{alt_dp/(2*n_dp):.6f}"
            if sensitivity else "NA"
        )
    })

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
output = outdir / f"BBG_22animal_final_locus_QC_{timestamp}.tsv"

with output.open("x", newline="") as f:
    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(results[0])
    )
    writer.writeheader()
    writer.writerows(results)

total_ft = sum(
    r["callable_ft_pass"] == "1"
    for r in rows
)
total_dp = sum(
    r["callable_ft_pass"] == "1"
    and depth_at_least_10(r)
    for r in rows
)

print("\n===== FINAL COHORT STATISTICS =====")
print("ANIMALS:", 22)
print("SV_LOCI:", 96)
print("TOTAL_GENOTYPE_ATTEMPTS:", len(rows))
print("FT_PASS_TOTAL:", total_ft)
print("FT_PASS_PERCENT:", round(100*total_ft/len(rows),2))
print("FT_PASS_DP10_TOTAL:", total_dp)

print(
    "PRIMARY_ELIGIBLE_FT20:",
    sum(r["primary_ft20_eligible"] for r in results)
)
print(
    "SENSITIVITY_ELIGIBLE_DP10_20:",
    sum(r["sensitivity_dp10_20_eligible"] for r in results)
)

print("FINAL_LOCUS_QC_REPORT:", output)
print("BBG_FINAL_22_ANIMAL_QC_COMPLETE")
PY

echo "===== FINAL GATK STATUS ====="
find /data/gatk/Kh5_recovery_20261117/logs \
    -maxdepth 1 -name '*.done' | wc -l

echo "===== FINAL SERVER STORAGE ====="
df -h /data

echo "BBG_REMAINING_COHORT_BATCH_COMPLETE"

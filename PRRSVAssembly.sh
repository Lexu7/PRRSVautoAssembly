#!/usr/bin/env bash
set -euo pipefail

# Fix GLIBCXX / conda runtime
export LD_LIBRARY_PATH="/ifs1/Software/Miniconda3/lib:${LD_LIBRARY_PATH:-}"

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date -Is)] $*"; }

on_error() {
  echo
  echo "ERROR: Pipeline failed near line $1. See log for details."
}
trap 'on_error $LINENO' ERR

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_DIR="${SCRIPT_DIR}/scripts"

# Defaults
THREADS=8
GENOME_SIZE="18k"
MIN_LEN=300
MAX_LEN=3000
COV_WARN_DP=20
COV_WARN_LOW_FRAC=0.10
COV_WARN_ZERO_FRAC=0.01
INDEL_MPILEUP_MAX_DEPTH=1000
OUTDIR_OVERRIDE=""
RESUME=false

# Data paths (default: script dir)
DB_FULL="${SCRIPT_DIR}/references/prrsv_blast"
REF_FASTA="${SCRIPT_DIR}/references/prrsv_refs.fasta"

# Tools
CANU_BIN="/ifs1/Software/canu-2.2/bin/canu"
FASTP_BIN="fastp"
BLASTN_BIN="blastn"
MINIMAP2_BIN="minimap2"
SAMTOOLS_BIN="samtools"
BCFTOOLS_BIN="bcftools"
SEQKIT_BIN="seqkit"
PYTHON_BIN="python3"

usage() {
  cat <<EOF
Usage:
  $0 [options] reads.fastq[.gz]

Options:
  -o OUTDIR      output directory (default: <sample>_prrsv_assembly)
  -t THREADS     threads (default: $THREADS)

  --min-len N    fastp min length (default: $MIN_LEN)
  --max-len N    fastp max length (default: $MAX_LEN)
  --genome-size X  canu genomeSize (default: $GENOME_SIZE)

  --cov-warn-dp N          depth < N considered low coverage (default: $COV_WARN_DP)
  --cov-warn-low-frac X    warn if low-coverage fraction >= X (default: $COV_WARN_LOW_FRAC)
  --cov-warn-zero-frac X   warn if zero-coverage fraction >= X (default: $COV_WARN_ZERO_FRAC)
  --indel-mpileup-max-depth N
                          max per-file mpileup depth for indel QC (default: $INDEL_MPILEUP_MAX_DEPTH)
  --resume                allow reuse of an existing output directory

  --db-full PREFIX       BLAST database prefix for full PRRSV genomes
  --ref-fasta PATH

  -h, --help
EOF
}

ARGS=()
require_value() {
  [[ $# -ge 2 && -n "${2:-}" ]] || die "Option $1 requires a value"
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o) require_value "$@"; OUTDIR_OVERRIDE="$2"; shift 2;;
    -t) require_value "$@"; THREADS="$2"; shift 2;;
    --min-len) require_value "$@"; MIN_LEN="$2"; shift 2;;
    --max-len) require_value "$@"; MAX_LEN="$2"; shift 2;;
    --genome-size) require_value "$@"; GENOME_SIZE="$2"; shift 2;;

    --cov-warn-dp) require_value "$@"; COV_WARN_DP="$2"; shift 2;;
    --cov-warn-low-frac) require_value "$@"; COV_WARN_LOW_FRAC="$2"; shift 2;;
    --cov-warn-zero-frac) require_value "$@"; COV_WARN_ZERO_FRAC="$2"; shift 2;;
    --indel-mpileup-max-depth) require_value "$@"; INDEL_MPILEUP_MAX_DEPTH="$2"; shift 2;;

    --db-full) require_value "$@"; DB_FULL="$2"; shift 2;;
    --ref-fasta) require_value "$@"; REF_FASTA="$2"; shift 2;;
    --resume) RESUME=true; shift;;

    -h|--help) usage; exit 0;;
    -*) die "Unknown option: $1";;
    *) ARGS+=("$1"); shift;;
  esac
done

[[ ${#ARGS[@]} -eq 1 ]] || { usage; die "Exactly one input FASTQ is required"; }
INPUT_FASTQ="$(readlink -f "${ARGS[0]}")"
[[ -s "$INPUT_FASTQ" ]] || die "Input fastq missing or empty: $INPUT_FASTQ"

# Sample name
BASENAME="$(basename "$INPUT_FASTQ")"
SAMPLE="$BASENAME"
SAMPLE=${SAMPLE%.fastq.gz}; SAMPLE=${SAMPLE%.fq.gz}; SAMPLE=${SAMPLE%.fastq}; SAMPLE=${SAMPLE%.fq}

OUTDIR="${OUTDIR_OVERRIDE:-${SAMPLE}_prrsv_assembly}"
if [[ -e "$OUTDIR" && ! -d "$OUTDIR" ]]; then
  die "Output path exists but is not a directory: $OUTDIR"
fi
if [[ -d "$OUTDIR" ]] && [[ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]] && [[ "$RESUME" != true ]]; then
  die "Output directory is not empty: $OUTDIR. Use --resume only when intentionally continuing this run."
fi
mkdir -p "$OUTDIR"
cd "$OUTDIR"

LOG_FILE="${SAMPLE}.run.log"
exec > >(tee -a "$LOG_FILE") 2>&1

log "=== PRRSV Pipeline Start ==="
log "ScriptDir: $SCRIPT_DIR"
log "Sample:    $SAMPLE"
log "Input:     $INPUT_FASTQ"
log "OutDir:    $(pwd)"
log "Resume:    $RESUME"

# Preconditions
command -v "$FASTP_BIN" >/dev/null || die "Missing tool: $FASTP_BIN"
command -v "$BLASTN_BIN" >/dev/null || die "Missing tool: $BLASTN_BIN"
command -v "$MINIMAP2_BIN" >/dev/null || die "Missing tool: $MINIMAP2_BIN"
command -v "$SAMTOOLS_BIN" >/dev/null || die "Missing tool: $SAMTOOLS_BIN"
command -v "$BCFTOOLS_BIN" >/dev/null || die "Missing tool: $BCFTOOLS_BIN"
command -v "$SEQKIT_BIN" >/dev/null || die "Missing tool: $SEQKIT_BIN"
command -v "$PYTHON_BIN" >/dev/null || die "Missing tool: $PYTHON_BIN"
[[ -x "$CANU_BIN" ]] || die "Canu not executable: $CANU_BIN"
[[ -s "$REF_FASTA" ]] || die "Missing REF_FASTA: $REF_FASTA"
ls "${DB_FULL}".n* >/dev/null 2>&1 || die "Missing BLAST db: $DB_FULL"
[[ -s "${PY_DIR}/select_best_reference.py" ]] || die "Missing ${PY_DIR}/select_best_reference.py"
[[ -s "${PY_DIR}/coverage_qc.py" ]] || die "Missing ${PY_DIR}/coverage_qc.py"
[[ -s "${PY_DIR}/mask_zero_coverage.py" ]] || die "Missing ${PY_DIR}/mask_zero_coverage.py"
[[ -s "${PY_DIR}/fill_zero_coverage_with_ref.py" ]] || die "Missing ${PY_DIR}/fill_zero_coverage_with_ref.py"
[[ -s "${PY_DIR}/resolve_degenerate_bases.py" ]] || die "Missing ${PY_DIR}/resolve_degenerate_bases.py"
[[ -s "${PY_DIR}/indel_qc.py" ]] || die "Missing ${PY_DIR}/indel_qc.py"
[[ -s "${PY_DIR}/parse_cs_stats.py" ]] || die "Missing ${PY_DIR}/parse_cs_stats.py"
[[ -s "${PY_DIR}/generate_report.py" ]] || die "Missing ${PY_DIR}/generate_report.py"

WORKDIR="work"
QC_DIR="${WORKDIR}/qc"
mkdir -p "$WORKDIR" "$QC_DIR"

# Key output names (root)
FINAL_FASTA="${SAMPLE}.final_consensus.fasta"
REF_FILLED_FASTA="${SAMPLE}.final_consensus.ref_filled.fasta"
FASTP_HTML="${SAMPLE}.fastp_report.html"
BEST_REF_ID_TXT="${SAMPLE}.best_ref_id.txt"
REPORT_CN="${SAMPLE}.report_CN.txt"
REPORT_EN="${SAMPLE}.report_EN.txt"
WARNING_TXT="${SAMPLE}.warning.txt"

########################################
# Step 1: fastp
########################################
log "[1] fastp..."
TRIMMED_FASTQ="${WORKDIR}/${SAMPLE}.clean.fastq"
FASTP_JSON="${WORKDIR}/${SAMPLE}.fastp.json"
$FASTP_BIN \
  -i "$INPUT_FASTQ" -o "$TRIMMED_FASTQ" \
  --trim_front1 20 --trim_tail1 20 \
  --low_complexity_filter \
  --disable_adapter_trimming --disable_quality_filtering \
  --length_required "$MIN_LEN" --length_limit "$MAX_LEN" \
  --thread "$THREADS" \
  -h "$FASTP_HTML" -j "$FASTP_JSON" \
  > "${WORKDIR}/${SAMPLE}.fastp.log" 2>&1
[[ -s "$TRIMMED_FASTQ" ]] || die "fastp output empty: $TRIMMED_FASTQ"
[[ -s "$FASTP_HTML" ]] || log "WARNING: fastp HTML missing: $FASTP_HTML"

########################################
# Step 2: Canu
########################################
log "[2] Canu de novo..."
CANU_OUTDIR="${WORKDIR}/${SAMPLE}.canu30"
CANU_PREFIX="$SAMPLE"
set +e
$CANU_BIN -nanopore "$TRIMMED_FASTQ" minOverlapLength=30 \
  -d "$CANU_OUTDIR" -p "$CANU_PREFIX" genomeSize="$GENOME_SIZE" \
  > "${WORKDIR}/${SAMPLE}.canu.log" 2>&1
CANU_EXIT=$?
set -e
ASSEMBLY_FASTA="${CANU_OUTDIR}/${CANU_PREFIX}.contigs.fasta"
if [[ $CANU_EXIT -ne 0 || ! -s "$ASSEMBLY_FASTA" ]]; then
  log "ERROR: Canu failed or contigs missing."
  log "  canu exit code: $CANU_EXIT"
  log "  expected contigs: $ASSEMBLY_FASTA"
  log "  see: ${WORKDIR}/${SAMPLE}.canu.log"
  exit 1
fi

########################################
# Step 3: BLAST
########################################
log "[3] BLAST..."
BLAST_FULL="${WORKDIR}/${SAMPLE}.blast_full.tsv"

$BLASTN_BIN -task megablast -query "$ASSEMBLY_FASTA" -db "$DB_FULL" -evalue 1e-20 -max_target_seqs 50 \
  -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
  -num_threads "$THREADS" > "$BLAST_FULL"

########################################
# Step 4: select best full-genome reference
########################################
log "[4] Select best full-genome reference..."
BEST_REF_TXT="${WORKDIR}/${SAMPLE}.best_ref.txt"
BEST_REF_SCORES_TSV="${WORKDIR}/${SAMPLE}.best_ref_scores.tsv"
$PYTHON_BIN "${PY_DIR}/select_best_reference.py" \
  --blast "$BLAST_FULL" --out "$BEST_REF_TXT" --scores "$BEST_REF_SCORES_TSV"
BEST_REF_ID="$(head -n 1 "$BEST_REF_TXT" | tr -d '\r\n')"
[[ -n "$BEST_REF_ID" && "$BEST_REF_ID" != "No_Match_Found" ]] || die "No suitable reference found."
#echo "$BEST_REF_ID" > "$BEST_REF_ID_TXT"
log "Best ref: $BEST_REF_ID"

########################################
# Step 5: extract ref
########################################
log "[5] Extract best reference fasta..."
BEST_REF_FASTA="${WORKDIR}/${SAMPLE}.best_ref.fasta"
$SEQKIT_BIN grep -p "$BEST_REF_ID" "$REF_FASTA" > "$BEST_REF_FASTA"
[[ -s "$BEST_REF_FASTA" ]] || die "Failed to extract best reference: $BEST_REF_ID"
BEST_REF_COUNT="$($SEQKIT_BIN seq -n "$BEST_REF_FASTA" | wc -l)"
[[ "$BEST_REF_COUNT" -eq 1 ]] || die "Expected exactly one best reference for $BEST_REF_ID, got $BEST_REF_COUNT"

########################################
# Step 6: consensus round1
########################################
log "[6] Consensus round1..."
ROUND1_BAM="${WORKDIR}/${SAMPLE}.round1.bam"
ROUND1_VCF="${WORKDIR}/${SAMPLE}.round1.vcf.gz"
CONS1_FASTA="${WORKDIR}/${SAMPLE}.consensus.round1.fa"

$MINIMAP2_BIN -ax map-ont "$BEST_REF_FASTA" "$TRIMMED_FASTQ" \
  | $SAMTOOLS_BIN sort -@ "$THREADS" -o "$ROUND1_BAM" -
$SAMTOOLS_BIN index "$ROUND1_BAM"
$BCFTOOLS_BIN mpileup -f "$BEST_REF_FASTA" "$ROUND1_BAM" | $BCFTOOLS_BIN call --ploidy 1 -mv -Oz -o "$ROUND1_VCF"
$BCFTOOLS_BIN index "$ROUND1_VCF"
$BCFTOOLS_BIN consensus -f "$BEST_REF_FASTA" "$ROUND1_VCF" > "$CONS1_FASTA"
[[ -s "$CONS1_FASTA" ]] || die "Round1 consensus empty."

########################################
# Step 7: consensus round2 (final)
########################################
log "[7] Consensus round2..."
ROUND2_BAM="${WORKDIR}/${SAMPLE}.round2.bam"
ROUND2_VCF="${WORKDIR}/${SAMPLE}.round2.vcf.gz"
CONS2_FASTA="${WORKDIR}/${SAMPLE}.consensus.round2.fa"

$MINIMAP2_BIN -ax map-ont "$CONS1_FASTA" "$TRIMMED_FASTQ" \
  | $SAMTOOLS_BIN sort -@ "$THREADS" -o "$ROUND2_BAM" -
$SAMTOOLS_BIN index "$ROUND2_BAM"
$BCFTOOLS_BIN mpileup -f "$CONS1_FASTA" "$ROUND2_BAM" | $BCFTOOLS_BIN call --ploidy 1 -mv -Oz -o "$ROUND2_VCF"
$BCFTOOLS_BIN index "$ROUND2_VCF"
$BCFTOOLS_BIN consensus -f "$CONS1_FASTA" "$ROUND2_VCF" > "$CONS2_FASTA"
[[ -s "$CONS2_FASTA" ]] || die "Round2 consensus empty."
cp -f "$CONS2_FASTA" "$FINAL_FASTA"

# Rename FASTA header to sample name
$SEQKIT_BIN replace -p '.*' -r "$SAMPLE" "$FINAL_FASTA" -o "${FINAL_FASTA}.tmp"
mv "${FINAL_FASTA}.tmp" "$FINAL_FASTA"

########################################
# Step 8: mask zero-coverage bases, then QC
########################################
log "[8] Mask zero-coverage bases..."
PRE_MASK_BAM="${QC_DIR}/${SAMPLE}.reads_vs_final.premask.bam"
PRE_MASK_DEPTH_TSV="${QC_DIR}/${SAMPLE}.depth.premask.tsv"
MASK_STATS_KV="${QC_DIR}/${SAMPLE}.zero_coverage_mask.kv"
REF_FILL_STATS_KV="${QC_DIR}/${SAMPLE}.ref_fill_zero_coverage.kv"
REF_TO_FINAL_PAF="${QC_DIR}/${SAMPLE}.best_ref_to_round2.paf"
$MINIMAP2_BIN -ax map-ont "$FINAL_FASTA" "$TRIMMED_FASTQ" \
  | $SAMTOOLS_BIN sort -@ "$THREADS" -o "$PRE_MASK_BAM" -
$SAMTOOLS_BIN index "$PRE_MASK_BAM"
$SAMTOOLS_BIN depth -aa "$PRE_MASK_BAM" > "$PRE_MASK_DEPTH_TSV"
[[ -s "$PRE_MASK_DEPTH_TSV" ]] || die "Pre-mask depth file empty: $PRE_MASK_DEPTH_TSV"

PRE_MASK_ZERO_COUNT="$(awk -F '\t' '$3==0{n++} END{print n+0}' "$PRE_MASK_DEPTH_TSV")"
if [[ "$PRE_MASK_ZERO_COUNT" -gt 0 ]]; then
  $MINIMAP2_BIN -cx asm5 --secondary=no -t "$THREADS" "$BEST_REF_FASTA" "$FINAL_FASTA" > "$REF_TO_FINAL_PAF"
  [[ -s "$REF_TO_FINAL_PAF" ]] || die "Reference-to-consensus alignment is empty: $REF_TO_FINAL_PAF"
  $PYTHON_BIN "${PY_DIR}/fill_zero_coverage_with_ref.py" \
    --consensus "$FINAL_FASTA" --best-ref "$BEST_REF_FASTA" --depth "$PRE_MASK_DEPTH_TSV" --alignment-paf "$REF_TO_FINAL_PAF" --sample "$SAMPLE" \
    --out "$REF_FILLED_FASTA" --stats "$REF_FILL_STATS_KV"
  [[ -s "$REF_FILLED_FASTA" ]] || die "Ref-filled consensus empty: $REF_FILLED_FASTA"
  [[ -s "$REF_FILL_STATS_KV" ]] || die "Ref-filled stats empty: $REF_FILL_STATS_KV"
  log "Best-ref-filled zero-coverage version stats: $(tr '\n' ';' < "$REF_FILL_STATS_KV")"
else
  rm -f "$REF_FILLED_FASTA" "$REF_FILL_STATS_KV"
  log "Best-ref-filled zero-coverage version skipped: no depth=0 positions."
fi

$PYTHON_BIN "${PY_DIR}/mask_zero_coverage.py" \
  --fasta "$FINAL_FASTA" --depth "$PRE_MASK_DEPTH_TSV" --sample "$SAMPLE" \
  --out "${FINAL_FASTA}.tmp" --stats "$MASK_STATS_KV"
mv "${FINAL_FASTA}.tmp" "$FINAL_FASTA"
FINAL_LENGTH="$(awk -F '\t' '$1=="final_length"{print $2}' "$MASK_STATS_KV")"
[[ -n "$FINAL_LENGTH" && "$FINAL_LENGTH" -gt 0 ]] || die "Final consensus is empty after zero-coverage trimming."
log "Zero-coverage masking stats: $(tr '\n' ';' < "$MASK_STATS_KV")"

log "[9] QC: coverage + corrected bases..."
# coverage depth after zero-coverage masking/trimming
READS_VS_FINAL_BAM="${QC_DIR}/${SAMPLE}.reads_vs_final.bam"
DEPTH_TSV="${QC_DIR}/${SAMPLE}.depth.tsv"
FINAL_MPILEUP="${QC_DIR}/${SAMPLE}.final.mpileup"
DEGEN_STATS_KV="${QC_DIR}/${SAMPLE}.degenerate_resolve.kv"
$MINIMAP2_BIN -ax map-ont "$FINAL_FASTA" "$TRIMMED_FASTQ" \
  | $SAMTOOLS_BIN sort -@ "$THREADS" -o "$READS_VS_FINAL_BAM" -
$SAMTOOLS_BIN index "$READS_VS_FINAL_BAM"
$SAMTOOLS_BIN depth -aa "$READS_VS_FINAL_BAM" > "$DEPTH_TSV"
[[ -s "$DEPTH_TSV" ]] || die "Depth file empty: $DEPTH_TSV"
if awk '!/^>/ && /[RYSWKMBDHV]/ {found=1; exit} END {exit !found}' "$FINAL_FASTA"; then
  $SAMTOOLS_BIN mpileup -aa -Q 0 -f "$FINAL_FASTA" "$READS_VS_FINAL_BAM" > "$FINAL_MPILEUP"
  [[ -s "$FINAL_MPILEUP" ]] || die "Final mpileup file empty: $FINAL_MPILEUP"
  $PYTHON_BIN "${PY_DIR}/resolve_degenerate_bases.py" \
    --fasta "$FINAL_FASTA" --mpileup "$FINAL_MPILEUP" \
    --out "${FINAL_FASTA}.tmp" --stats "$DEGEN_STATS_KV"
else
  log "No IUPAC degenerate bases in final consensus; skip degenerate-base mpileup."
  $PYTHON_BIN "${PY_DIR}/resolve_degenerate_bases.py" \
    --fasta "$FINAL_FASTA" --no-degenerate \
    --out "${FINAL_FASTA}.tmp" --stats "$DEGEN_STATS_KV"
fi
mv "${FINAL_FASTA}.tmp" "$FINAL_FASTA"
log "Degenerate-base resolution stats: $(tr '\n' ';' < "$DEGEN_STATS_KV")"

INDEL_MPILEUP="${QC_DIR}/${SAMPLE}.indel_qc.d${INDEL_MPILEUP_MAX_DEPTH}.mpileup"
INDEL_QC_TSV="${QC_DIR}/${SAMPLE}.indel_qc.tsv"
INDEL_QC_KV="${QC_DIR}/${SAMPLE}.indel_qc.summary.tsv"
$SAMTOOLS_BIN mpileup -aa -A -d "$INDEL_MPILEUP_MAX_DEPTH" -Q 0 -f "$FINAL_FASTA" "$READS_VS_FINAL_BAM" > "$INDEL_MPILEUP"
[[ -s "$INDEL_MPILEUP" ]] || die "Indel QC mpileup file empty: $INDEL_MPILEUP"
$PYTHON_BIN "${PY_DIR}/indel_qc.py" \
  --mpileup "$INDEL_MPILEUP" --sample "$SAMPLE" \
  --out "$INDEL_QC_TSV" --summary "$INDEL_QC_KV"
[[ -s "$INDEL_QC_TSV" ]] || die "Indel QC TSV empty: $INDEL_QC_TSV"
[[ -s "$INDEL_QC_KV" ]] || die "Indel QC summary empty: $INDEL_QC_KV"
log "Indel QC stats: $(tr '\n' ';' < "$INDEL_QC_KV")"

COV_KV="${QC_DIR}/${SAMPLE}.coverage_qc.kv"
$PYTHON_BIN "${PY_DIR}/coverage_qc.py" \
  --depth "$DEPTH_TSV" --low-dp "$COV_WARN_DP" --low-frac "$COV_WARN_LOW_FRAC" --zero-frac "$COV_WARN_ZERO_FRAC" \
  -o "$COV_KV"

# corrected bases via cs-tag
PAF_REF_R1="${QC_DIR}/${SAMPLE}.ref_to_round1.paf"
PAF_R1_R2="${QC_DIR}/${SAMPLE}.round1_to_round2.paf"
CS1_KV="${QC_DIR}/${SAMPLE}.cs_round1.kv"
CS2_KV="${QC_DIR}/${SAMPLE}.cs_round2.kv"
$MINIMAP2_BIN -cx asm5 --cs -t "$THREADS" "$CONS1_FASTA" "$BEST_REF_FASTA" > "$PAF_REF_R1"
$MINIMAP2_BIN -cx asm5 --cs -t "$THREADS" "$FINAL_FASTA" "$CONS1_FASTA" > "$PAF_R1_R2"
$PYTHON_BIN "${PY_DIR}/parse_cs_stats.py" --paf "$PAF_REF_R1" -o "$CS1_KV"
$PYTHON_BIN "${PY_DIR}/parse_cs_stats.py" --paf "$PAF_R1_R2" -o "$CS2_KV"

########################################
# Step 10: Report (CN/EN + warning)
########################################
log "[10] Report..."
REPORT_ARGS=(
  --sample "$SAMPLE"
  --input-fastq "$INPUT_FASTQ"
  --final-fasta "$FINAL_FASTA"
  --best-ref-id "$BEST_REF_ID"
  --fastp-html "$FASTP_HTML"
  --fastp-json "$FASTP_JSON"
  --coverage-kv "$COV_KV"
  --cs-round1-kv "$CS1_KV"
  --cs-round2-kv "$CS2_KV"
  --mask-stats-kv "$MASK_STATS_KV"
  --degenerate-stats-kv "$DEGEN_STATS_KV"
  --indel-qc-kv "$INDEL_QC_KV"
  --indel-qc-tsv "$INDEL_QC_TSV"
  --out-cn "$REPORT_CN"
  --out-en "$REPORT_EN"
  --out-warning "$WARNING_TXT"
)
if [[ -s "$REF_FILL_STATS_KV" ]]; then
  REPORT_ARGS+=(--ref-fill-stats-kv "$REF_FILL_STATS_KV" --ref-filled-fasta "$REF_FILLED_FASTA")
fi
$PYTHON_BIN "${PY_DIR}/generate_report.py" "${REPORT_ARGS[@]}"

log "=== Finished ==="
log "Key outputs:"
log "  - $FINAL_FASTA"
if [[ -s "$REF_FILLED_FASTA" ]]; then
  log "  - $REF_FILLED_FASTA"
fi
log "  - $FASTP_HTML"
#log "  - $BEST_REF_ID_TXT"
log "  - $REPORT_CN"
log "  - $REPORT_EN"
log "  - $WARNING_TXT"
log "  - $INDEL_QC_TSV"
log "  - $LOG_FILE"
log "Intermediates in: $WORKDIR/"

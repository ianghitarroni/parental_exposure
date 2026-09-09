#!/usr/bin/env bash
set -euo pipefail

# Complete preprocessing and quantification workflow used to generate the
# featureCounts matrix consumed by scripts/01_deseq2_analysis.R.

DIR_RAW="${DIR_RAW:-00_RawData}"
DIR_TRIM="${DIR_TRIM:-01_Trimmed}"
DIR_QC="${DIR_QC:-02_QC}"
DIR_ALIGN="${DIR_ALIGN:-03_Aligned_STAR}"
DIR_COUNTS="${DIR_COUNTS:-04_Counts}"
DIR_GENOME="${DIR_GENOME:-Ref_Genome_mm39}"
STAR_INDEX="${STAR_INDEX:-${DIR_GENOME}/star_index}"
GTF_FILE="${GTF_FILE:-${DIR_GENOME}/gencode.vM33.annotation.gtf}"
THREADS="${THREADS:-4}"
BAM_SORT_RAM="${BAM_SORT_RAM:-2000000000}"

SAMPLES=("A2" "A3" "B3" "B5")

mkdir -p \
    "$DIR_TRIM" \
    "$DIR_QC/raw" \
    "$DIR_QC/trimmed" \
    "$DIR_ALIGN" \
    "$DIR_COUNTS"

for command in fastqc fastp STAR samtools featureCounts; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $command" >&2
        exit 1
    fi
done

if [[ ! -d "$STAR_INDEX" ]]; then
    echo "ERROR: STAR index directory not found: $STAR_INDEX" >&2
    exit 1
fi

if [[ ! -f "$GTF_FILE" ]]; then
    echo "ERROR: GTF annotation not found: $GTF_FILE" >&2
    exit 1
fi

capture_version() {
    local tool="$1"
    shift
    local version
    version="$("$@" 2>&1)"
    version="${version%%$'\n'*}"
    version="${version//$'\t'/ }"
    printf '%s\t%s\n' "$tool" "$version"
}

{
    printf 'item\tvalue\n'
    capture_version "FastQC" fastqc --version
    capture_version "fastp" fastp --version
    capture_version "STAR" STAR --version
    capture_version "samtools" samtools --version
    capture_version "featureCounts" featureCounts -v
    printf 'reference_genome\tGRCm39/mm39\n'
    printf 'annotation\tGENCODE vM33\n'
    printf 'star_index\t%s\n' "$STAR_INDEX"
    printf 'gtf_file\t%s\n' "$GTF_FILE"
    printf 'threads\t%s\n' "$THREADS"
    printf 'run_started_utc\t%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
} > "$DIR_QC/software_versions.tsv"

run_fastqc_pair() {
    local output_dir="$1"
    local read_1="$2"
    local read_2="$3"
    local report_1="$4"
    local report_2="$5"

    if [[ -f "$report_1" && -f "$report_2" ]]; then
        echo "FastQC reports already exist: $report_1, $report_2"
        return
    fi

    fastqc \
        --threads "$THREADS" \
        --outdir "$output_dir" \
        "$read_1" "$read_2"
}

for sample in "${SAMPLES[@]}"; do
    raw_r1="$DIR_RAW/${sample}_R1.fastq.gz"
    raw_r2="$DIR_RAW/${sample}_R2.fastq.gz"
    trim_r1="$DIR_TRIM/${sample}_trimmed_R1.fastq.gz"
    trim_r2="$DIR_TRIM/${sample}_trimmed_R2.fastq.gz"
    bam="$DIR_ALIGN/${sample}_Aligned.sortedByCoord.out.bam"
    prefix="$DIR_ALIGN/${sample}_"
    tmp_dir="$DIR_ALIGN/${sample}_STARtmp"

    echo "========================================================"
    echo "Processing sample: $sample"
    echo "========================================================"

    if [[ -f "$raw_r1" && -f "$raw_r2" ]]; then
        echo "[1/4] Running FastQC on raw reads"
        run_fastqc_pair \
            "$DIR_QC/raw" \
            "$raw_r1" \
            "$raw_r2" \
            "$DIR_QC/raw/${sample}_R1_fastqc.html" \
            "$DIR_QC/raw/${sample}_R2_fastqc.html"
    elif [[ ! -f "$bam" ]]; then
        echo "ERROR: paired FASTQ files not found for $sample" >&2
        exit 1
    else
        echo "WARNING: raw FASTQ files are unavailable; existing BAM will be reused" >&2
    fi

    if [[ ! -f "$bam" ]]; then
        if [[ ! -f "$trim_r1" || ! -f "$trim_r2" ]]; then
            echo "[2/4] Running fastp"
            fastp \
                -i "$raw_r1" \
                -I "$raw_r2" \
                -o "$trim_r1" \
                -O "$trim_r2" \
                -h "$DIR_QC/${sample}_fastp.html" \
                -j "$DIR_QC/${sample}_fastp.json" \
                --detect_adapter_for_pe \
                -w "$THREADS"
        else
            echo "[2/4] Trimmed FASTQ files already exist"
        fi

        echo "[3/4] Running FastQC on filtered reads"
        run_fastqc_pair \
            "$DIR_QC/trimmed" \
            "$trim_r1" \
            "$trim_r2" \
            "$DIR_QC/trimmed/${sample}_trimmed_R1_fastqc.html" \
            "$DIR_QC/trimmed/${sample}_trimmed_R2_fastqc.html"

        echo "[4/4] Running STAR"
        rm -rf "$tmp_dir"
        STAR \
            --runThreadN "$THREADS" \
            --genomeDir "$STAR_INDEX" \
            --readFilesIn "$trim_r1" "$trim_r2" \
            --readFilesCommand zcat \
            --outFileNamePrefix "$prefix" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outSAMattributes Standard \
            --limitBAMsortRAM "$BAM_SORT_RAM" \
            --outTmpDir "$tmp_dir"

        samtools index "$bam"

        # The original execution removed trimmed FASTQ files after successful
        # alignment to reduce local disk usage. Set KEEP_TRIMMED=1 to retain them.
        if [[ "${KEEP_TRIMMED:-0}" != "1" ]]; then
            rm -f "$trim_r1" "$trim_r2"
        fi
        rm -rf "$tmp_dir"
    else
        echo "Existing BAM found; skipping fastp and STAR"
        if [[ ! -f "${bam}.bai" ]]; then
            samtools index "$bam"
        fi
    fi
done

echo "Running featureCounts"
bam_files=()
for sample in "${SAMPLES[@]}"; do
    bam="$DIR_ALIGN/${sample}_Aligned.sortedByCoord.out.bam"
    if [[ ! -f "$bam" ]]; then
        echo "ERROR: expected BAM not found: $bam" >&2
        exit 1
    fi
    bam_files+=("$bam")
done

featureCounts \
    -T "$THREADS" \
    -p \
    --countReadPairs \
    -t exon \
    -g gene_id \
    -a "$GTF_FILE" \
    -o "$DIR_COUNTS/counts_matrix.txt" \
    "${bam_files[@]}"

echo "Pipeline completed: $DIR_COUNTS/counts_matrix.txt"
echo "QC reports: $DIR_QC"

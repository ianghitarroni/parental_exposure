#!/usr/bin/env bash
set -euo pipefail

# Complete preprocessing and quantification workflow used to generate the
# featureCounts matrix consumed by scripts/01_deseq2_analysis.R.

DIR_RAW="${DIR_RAW:-00_RawData}"
DIR_TRIM="${DIR_TRIM:-01_Trimmed}"
DIR_ALIGN="${DIR_ALIGN:-03_Aligned_STAR}"
DIR_COUNTS="${DIR_COUNTS:-04_Counts}"
DIR_GENOME="${DIR_GENOME:-Ref_Genome_mm39}"
STAR_INDEX="${STAR_INDEX:-${DIR_GENOME}/star_index}"
GTF_FILE="${GTF_FILE:-${DIR_GENOME}/gencode.vM33.annotation.gtf}"
THREADS="${THREADS:-4}"
BAM_SORT_RAM="${BAM_SORT_RAM:-2000000000}"

SAMPLES=("A2" "A3" "B3" "B5")

mkdir -p "$DIR_TRIM" "$DIR_ALIGN" "$DIR_COUNTS"

for command in fastp STAR samtools featureCounts; do
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

    if [[ ! -f "$bam" ]]; then
        if [[ ! -f "$trim_r1" || ! -f "$trim_r2" ]]; then
            if [[ ! -f "$raw_r1" || ! -f "$raw_r2" ]]; then
                echo "ERROR: paired FASTQ files not found for $sample" >&2
                exit 1
            fi

            echo "[1/3] Running fastp"
            fastp \
                -i "$raw_r1" \
                -I "$raw_r2" \
                -o "$trim_r1" \
                -O "$trim_r2" \
                -h "$DIR_TRIM/${sample}_fastp.html" \
                -j "$DIR_TRIM/${sample}_fastp.json" \
                --detect_adapter_for_pe \
                -w "$THREADS"
        else
            echo "[1/3] Trimmed FASTQ files already exist"
        fi

        echo "[2/3] Running STAR"
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

        # The original execution removed trimmed FASTQ files after a successful
        # alignment to reduce local disk usage. Set KEEP_TRIMMED=1 to retain them.
        if [[ "${KEEP_TRIMMED:-0}" != "1" ]]; then
            rm -f "$trim_r1" "$trim_r2"
        fi
        rm -rf "$tmp_dir"
    else
        echo "[1-2/3] Existing BAM found; skipping trimming and alignment"
        if [[ ! -f "${bam}.bai" ]]; then
            samtools index "$bam"
        fi
    fi
done

echo "[3/3] Running featureCounts"
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

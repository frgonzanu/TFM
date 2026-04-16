#!/usr/bin/env bash

THREADS=8

BASE_DIR="$HOME/fgonzaleznun/TFM"
RAW_DIR="$BASE_DIR/data/FastQ/"
REF_DIR="$BASE_DIR/ref/hg38/star_index_hg38"
OUT_DIR="$BASE_DIR/results/star/raw_alignment_relaxed"

echo "BASE_DIR: $BASE_DIR"
echo "RAW_DIR:  $RAW_DIR"
echo "REF_DIR:  $REF_DIR"
echo "OUT_DIR:  $OUT_DIR"
echo

########################################
# ALINEAMIENTO
########################################

echo "Buscando muestras SRR..."

shopt -s nullglob

for r1 in "$RAW_DIR"/SRR*_1.fastq.gz; do

    base=$(basename "$r1" _1.fastq.gz)
    r2="$RAW_DIR/${base}_2.fastq.gz"

    sample_out="$OUT_DIR/$base"
    mkdir -p "$sample_out"

    echo ""
    echo "========================================"
    echo "Muestra: $base"
    echo "R1: $r1"
    echo "R2: $r2"
    echo "Salida: $sample_out"
    echo "========================================"

    STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$REF_DIR" \
        --readFilesIn "$r1" "$r2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$sample_out/" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outSAMattributes NH HI AS nM XS \
	--outFilterScoreMinOverLread 0.3 \
	--outFilterMatchNminOverLread 0.3 \
	--outFilterMatchNmin 30 \
        --outStd Log > "$sample_out/run.log" 2>&1

    echo "Done: $base"

done

echo ""
echo "Alignment completed"

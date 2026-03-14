#!/bin/bash

source ~/.local/share/mamba/etc/profile.d/mamba.sh
mamba activate rnaseq

RAW=/mnt/c/Users/User/Documents/Damola/rawData
TRIM=/mnt/c/Users/User/Documents/Damola/trimmed
REPORTS=/mnt/c/Users/User/Documents/Damola/fastp_reports

mkdir -p "$TRIM" "$REPORTS"

echo "Running fastp trimming using 16 CPUs (8 jobs × 2 threads)..."

ls $RAW/*_1.fastq.gz | parallel -j 8 '
    R1={};
    R2=${R1/_1.fastq.gz/_2.fastq.gz};
    base=$(basename $R1 _1.fastq.gz);

    fastp \
        -i $R1 -I $R2 \
        -o '"$TRIM"'/${base}_1_trimmed.fastq.gz \
        -O '"$TRIM"'/${base}_2_trimmed.fastq.gz \
        --html '"$REPORTS"'/${base}.html \
        --json '"$REPORTS"'/${base}.json \
        --thread 2
'

echo "fastp trimming complete."

mamba deactivate

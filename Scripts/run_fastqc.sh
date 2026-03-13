#!/bin/bash

# Load mamba inside script
source ~/.local/share/mamba/etc/profile.d/mamba.sh
mamba activate rnaseq

RAW=/mnt/c/Users/User/Documents/Damola/rawData
OUT=/mnt/c/Users/User/Documents/Damola/fastqc_raw

mkdir -p "$OUT"

echo "Running FastQC on RAW reads using 16 CPUs..."

ls $RAW/*.fastq.gz | parallel -j 16 "fastqc {} -o $OUT"

echo "FastQC (raw) complete."

mamba deactivate

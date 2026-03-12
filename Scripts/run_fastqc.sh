#!/bin/bash
# Load conda/mamba into this script
source ~/.local/share/mamba/etc/profile.d/mamba.sh
mamba activate rnaseq

RAW=/mnt/c/Users/User/Documents/Damola/rawData
OUT=/mnt/c/Users/User/Documents/Damola/fastqc_raw

mkdir -p $OUT

for fq in $RAW/*.fastq.gz
do
    echo "Running FastQC on $fq"
    fastqc $fq -o $OUT
done

mamba deactivate

#!/usr/bin/env bash
#author:    :Gregory Wickham
#date:      :20241113
#version    :1.0.0
#desc       :Script to split reads into separate bam files after demultiplexing
#usage		:bash separate_reads.sh
#===========================================================================================================

#populate barcodes.txt file which contains a list of barcodes from demultiplexed .bam file

source "$(sudo find ~ -maxdepth 4 -name conda.sh)" #find path to conda base environment
echo "conda sourced"
conda activate samtools

for bam in *.bam; do
    echo "Processing $bam..."  # Check which BAM file is being processed
    samtools view $bam | \
        awk '{for(i=12;i<=NF;i++) if($i ~ /^BC:Z:/) print substr($i,6)}' | \
        sort | \
        uniq > barcodes.txt
    echo "Barcodes for $bam written to barcodes.txt"

    while read barcode; do
        samtools view -h $bam | \
            awk -v bc=$barcode '$0 ~ /^@/ || $0 ~ "BC:Z:"bc {print}' | \
                samtools view -bo "${barcode}.bam" -
    done < barcodes.txt
done
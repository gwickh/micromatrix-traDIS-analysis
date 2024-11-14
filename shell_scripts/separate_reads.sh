#!/usr/bin/env bash
#author:    :Gregory Wickham
#date:      :20241113
#version    :1.0.0
#desc       :Script to split reads into separate bam files after demultiplexing
#usage		:bash separate_reads.sh
#===========================================================================================================

#populate barcodes.txt file which contains a list of barcodes from demultiplexed .bam file
samtools view reads_cat_minion.bam | \
    awk '{for(i=12;i<=NF;i++) if($i ~ /^BC:Z:/) print substr($i,6)}' | \
        sort | \
            uniq > barcodes.txt

while read barcode; do
    samtools view -b -h reads_cat_minion.bam -e "BC:Z:${barcode}" -o "${barcode}.bam"
done < barcodes.txt
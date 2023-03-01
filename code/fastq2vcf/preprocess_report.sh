#!/bin/bash

echo -e "SampleID\tTotal_paired_raw\tTotal_paired_filtered\tGC_perc\tDuplication_rate\tInsert_size\tMapping_perc"

while read sample; do 

    if [ -f "data/trim/${sample}_report.json" ]; then
        totalReadsRaw=$(grep "total_reads" data/trim/${sample}_report.json | tail -n 3 | head -n 1 | grep -o "[0-9]*")
        totalReadsFiltered=$(grep "total_reads" data/trim/${sample}_report.json | tail -n 1 | grep -o "[0-9]*")
        GC=$(grep "gc" data/trim/${sample}_report.json | tail -n 1 | grep -o "[0-9.]*")
        duplicationRate=$(grep -A 1 "duplication" data/trim/${sample}_report.json | tail -n 1 | grep -o "[0-9.]*")
        insertSize=$(grep "peak" data/trim/${sample}_report.json | grep -o "[0-9.]*")
    else
        totalReadsRaw="NA"
        totalReadsFiltered="NA"
        GC="NA"
        duplicationRate="NA"
        insertSize="NA"
    fi
    
    if [ -f "data/bam/${sample}.fixmate.markdup.sort.bam.flagstat" ]; then
        readBam=$(grep "0 primary$" data/bam/${sample}.fixmate.markdup.sort.bam.flagstat | grep -o "^[0-9]*")
        mappedPerc=$(grep "0 mapped" data/bam/${sample}.fixmate.markdup.sort.bam.flagstat | sed 's/.*(\([0-9.]*\)%.*/\1/')
    else
        readBam="NA"
        mappedPerc="NA"
    fi

    
    echo -e "$sample\t$totalReadsRaw\t$totalReadsFiltered\t$GC\t$duplicationRate\t$insertSize\t$mappedPerc"


done < $1
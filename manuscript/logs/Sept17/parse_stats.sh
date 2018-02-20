#!/bin/sh

INDR=./data/RSEM/hg19/Sept17
OUTDR=./manuscript/logs/Sept17

files=$(ls $INDR/*bam)

for file in $files
do
    samtools flagstat $file > $OUTDR/$(basename ${file/.transcript.bam/.logs})  &
done

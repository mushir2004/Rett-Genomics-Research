#!/bin/bash

echo " Phase 2: Aligning Patient Data to Reference..."

# 1. Create a folder for the results
mkdir -p results

# 2. Run BWA MEM alignment
# This compares the 'mutated' file to our 'index'
# Change the input file to the shredded one
bwa mem data/reference/chrX.fa data/patient/patient_mecp2_shredded.fa > results/alignment.sam

echo " Alignment complete! Generated results/alignment.sam"

# 3. Convert SAM to BAM (Binary format - smaller and faster)
samtools view -S -b results/alignment.sam > results/alignment.bam

# 4. Sort the BAM file (Required for visualization)
samtools sort results/alignment.bam -o results/alignment_sorted.bam

# 5. Index the BAM file
samtools index results/alignment_sorted.bam

echo " BAM file sorted and indexed. Ready for mutation hunting!"
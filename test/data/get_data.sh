#!/bin/bash

set -euo pipefail

WD=$(dirname $0)

which samtools seqtk

# Subset aligned BAMs to reads
# Cannot use samtools view on uBAM as cannot index.
samtools cat \
    /project/logsdon_shared/long_read_archive/unsorted/20250610_DMLseq_GM20355_ULK114/CENPA/20250610_1553_3C_PAY68234_eb3f5140/pod5/basecalling/20250610_DMLseq_GM20355_ULK114.bam \
    /project/logsdon_shared/long_read_archive/unsorted/20250610_DMLseq_GM20355_ULK114/CENPA_D2/20250611_1324_3C_PAY68234_360fbf6a/pod5/basecalling/20250610_DMLseq_GM20355_ULK114.bam |
    samtools bam2fq -T 'MM,ML' - | \
    seqtk subseq - <(samtools view /project/logsdon_shared/projects/HGSVC3/NA20355_Dimelo/results_NA20355/align/NA20355_treatment_1.0.bam | grep NA20355_chr8_haplotype1-0000017 | cut -f 1) | \
    bgzip > $WD/NA20355_chr8_CENPA.fq.gz

samtools cat \
    /project/logsdon_shared/long_read_archive/unsorted/20250610_DMLseq_GM20355_ULK114/IgG/20250610_1553_3F_PAY24505_6b59e78f/pod5/basecalling/20250610_DMLseq_GM20355_ULK114.bam \
    /project/logsdon_shared/long_read_archive/unsorted/20250610_DMLseq_GM20355_ULK114/IgG_D2/20250611_1324_3F_PAY24505_aba96014/pod5/basecalling/20250610_DMLseq_GM20355_ULK114.bam | \
    samtools bam2fq -T 'MM,ML' - | \
    seqtk subseq - <(samtools view /project/logsdon_shared/projects/HGSVC3/NA20355_Dimelo/results_NA20355/align/NA20355_control_1.0.bam | grep NA20355_chr8_haplotype1-0000017 | cut -f 1) | \
    bgzip > $WD/NA20355_chr8_IgG.fq.gz

# Get sequence of whole contig.
samtools faidx /project/logsdon_shared/projects/HGSVC3/new_65_asms_renamed/NA20355-asm-renamed-reort.fa NA20355_chr8_haplotype1-0000017 | bgzip > $WD/NA20355_chr8_haplotype1-0000017.fa.gz

# Get RM output
grep NA20355_chr8_haplotype1-0000017 /project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/RM/all_cens_chr8.annotation.fa.out | \
awk -v OFS="\t" '{match($1, "^(.+):", arr); $1=arr[1]; print}' | bgzip > test/data/NA20355_chr8_sat_annot.bed.gz

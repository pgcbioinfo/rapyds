#!/usr/bin/env bash
# dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
dir=$(pwd)
#bwa alignment
cd "$dir/$4"
pwd
hpc pgcbioinfo/bwa bwa aln $1 reads/$3_$2_read1.fastq > $3_$2_aln_sa1.sai
hpc pgcbioinfo/bwa bwa aln $1 reads/$3_$2_read2.fastq > $3_$2_aln_sa2.sai

hpc pgcbioinfo/bwa bwa sampe $1 $3_$2_aln_sa1.sai $3_$2_aln_sa2.sai reads/$3_$2_read1.fastq reads/$3_$2_read2.fastq > aligned_pairs_$3_$2.sam

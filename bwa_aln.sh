#!/usr/bin/env bash
# dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
dir=$(pwd)
#bwa alignment
cd "$dir/$4"
pwd
bwa aln $1 reads/$3_$2_read1.fastq > $3_$2_aln_sa1.sai
bwa aln $1 reads/$3_$2_read2.fastq > $3_$2_aln_sa2.sai

bwa sampe $1 $3_$2_aln_sa1.sai $3_$2_aln_sa2.sai reads/$3_$2_read1.fastq reads/$3_$2_read2.fastq > aligned_pairs_$3_$2.sam
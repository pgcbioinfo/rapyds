#!/usr/bin/env bash
dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
#bwa alignment
cd $dir/$4
pwd
bwa aln $1 $dir/reads/$3_$2_read1.fastq > $3_$2_aln_sa1.sai
bwa aln $1 $dir/reads/$3_$2_read2.fastq > $3_$2_aln_sa2.sai

bwa sampe $1 $dir/$4/$3_$2_aln_sa1.sai $dir/$4/$3_$2_aln_sa2.sai $dir/reads/$3_$2_read1.fastq $dir/reads/$3_$2_read2.fastq > $dir/$4/aligned_pairs_$3_$2.sam
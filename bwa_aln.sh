#!/usr/bin/env bash
#
# RApyDS
# Restriction site-associated DNA from Python-implemented Digestion Simulations
# https://github.com/pgcbioinfo/rapyds

# bwa_aln.sh
#

BWA_CMD=bwa

dir=$(pwd)
#bwa alignment
cd "$dir/$4"
pwd
${BWA_CMD} aln $1 reads/$3_$2_read1.fastq > $3_$2_aln_sa1.sai
${BWA_CMD} aln $1 reads/$3_$2_read2.fastq > $3_$2_aln_sa2.sai

${BWA_CMD} sampe $1 $3_$2_aln_sa1.sai $3_$2_aln_sa2.sai reads/$3_$2_read1.fastq reads/$3_$2_read2.fastq > aligned_pairs_$3_$2.sam

#!/usr/bin/env bash
#
# RApyDS
# Restriction site-associated DNA from Python-implemented Digestion Simulations
# https://github.com/pgcbioinfo/rapyds

# bwa_index.sh
#

BWA_CMD=bwa

dir=$(pwd)
#bwa alignment
cd "$dir/$3"
${BWA_CMD} index -p $2 $1

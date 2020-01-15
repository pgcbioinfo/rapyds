#!/usr/bin/env bash
#
# RApyDS
# Restriction site-associated DNA from Python-implemented Digestion Simulations
# https://github.com/pgcbioinfo/rapyds

# bwa_index.sh
#


dir=$(pwd)
#bwa alignment
cd "$dir/$3"
bwa index -p $2 $1

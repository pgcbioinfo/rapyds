#!/usr/bin/env bash
#dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
dir=$(pwd)
#bwa alignment
# rm -rf $dir/$3
# mkdir $dir/$3
cd "$dir/$3"
hpc pgcbioinfo/bwa bwa index -p $2 $1

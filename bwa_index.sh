#!/usr/bin/env bash
dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
#bwa alignment
rm -rf $dir/$3
mkdir $dir/$3
cd "$dir/$3"
bwa index -a bwtsw $dir/$1 -p $2

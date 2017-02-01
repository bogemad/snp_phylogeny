#!/bin/bash

base_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
temp=$base_path/.temp
logs=$base_path/logs

rm -rf $temp
mkdir -p $temp
mkdir -p $logs

cpus=$1
mem=$2

sed "s~xxxcpusxxx~$cpus~g" $base_path/scripts/qsub_submit.sh | sed "s~xxxmemxxx~$mem~g" | sed "s~xxxbasepathxxx~$base_path~g" > $temp/qsub.sh

qsub $temp/qsub.sh

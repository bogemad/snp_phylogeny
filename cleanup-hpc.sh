#!/bin/bash

base_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
temp=$base_path/.temp
logs=$base_path/logs

mkdir -p $temp

mkdir -p $logs

sed "s~xxxbasepathxxx~$base_path~g" $base_path/scripts/qsub_cleanup.sh > $temp/qsub.sh

qsub $temp/qsub.sh
sleep 10
rm -rf $temp/*

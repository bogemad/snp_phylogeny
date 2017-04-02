#!/bin/bash

base_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
temp=$base_path/.temp

mkdir -p $temp

sed "s~xxxbasepathxxx~$base_path~g" $base_path/scripts/qsub_make.sh > $temp/qsub.sh

qsub $temp/qsub.sh
sleep 5
rm -rf $temp/*

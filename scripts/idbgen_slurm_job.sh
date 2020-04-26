#!/usr/bin/env bash

path="$1"
num_cpus="$2"
type="$3"

echo "/software/ligandscout/ligandscout4.ilab/idbgen -i $path -o ${path:0:-4}.ldb -t $type -C $num_cpus"
/software/ligandscout/ligandscout4.ilab/idbgen -i $path -o ${path:0:-4}.ldb -t $type -C $num_cpus

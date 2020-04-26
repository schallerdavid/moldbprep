#!/usr/bin/env bash

directory_path="$1"
pattern="$2"
type="$3"
queue="cn-cpu"
num_cpus="5"

get_abs_filename() {
  # $1 : relative filename
  if [ -d "$(dirname "$1")" ]; then
    echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
  fi
}

for path in $directory_path/*; do
  if [[ $path == *$pattern* ]]; then
    sbatch -p$queue -n $num_cpus $(dirname $(get_abs_filename $0))/idbgen_slurm_job.sh $(get_abs_filename $path) $num_cpus $type
  fi
done
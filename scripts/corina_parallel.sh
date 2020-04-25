#!/usr/bin/env bash

queue="cn-cpu"
msi=8 # max stereo isomers
msc=4 # max stereo centers
directory_path="$1"

get_abs_filename() {
  # $1 : relative filename
  if [ -d "$(dirname "$1")" ]; then
    echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
  fi
}

for path in $directory_path/*; do
  if [[ $path == *".sdf" ]]; then
    sbatch -p$queue -n 1 $(dirname $(get_abs_filename $0))/corina_slurm_job.sh $(get_abs_filename $path) $msi $msc
  fi
done
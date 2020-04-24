#!/usr/bin/env bash

path="$1"
msi="$2"
msc="$3"

echo "corina -d stergen,msi=$msi,msc=$msc,preserve,wh,r2d -i t=sdf $path -o t=sdf ${path:0:-4}_3D.sdf"

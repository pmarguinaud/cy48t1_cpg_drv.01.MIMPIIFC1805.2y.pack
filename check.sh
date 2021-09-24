#!/bin/bash


for f in $(ls slurm-*.out | sort -r)
do
echo $f
grep NORMDIFF $f | grep -v 0.000000000e+00 | grep -v REF
#break
done



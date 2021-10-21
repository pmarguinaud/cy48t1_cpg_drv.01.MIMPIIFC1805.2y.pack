#!/bin/bash

let "i=0"

F=""

for f in $(ls slurm-*.out | sort -r)
do
F="$f $F"
let "i=$i+1"
if [ $i -eq 2 ]
then
  break
fi
done

grep NORMDIFF $F


#!/bin/bash

set -x

cd fypp

for f in *.fypp
do
  b=$(basename $f .fypp)
  /opt/softs/anaconda3/bin/fypp -m os -M . -m yaml -m field_config ./$b.fypp ./$b.F90

  if [ ! -f ../src/local/arpifs/module/$b.F90 ]
  then
    cp $b.F90 ../src/local/arpifs/module/$b.F90
  else
    cmp $b.F90 ../src/local/arpifs/module/$b.F90
    if [ $? -ne 0 ]
    then
      cp $b.F90 ../src/local/arpifs/module/$b.F90
    fi
  fi

done



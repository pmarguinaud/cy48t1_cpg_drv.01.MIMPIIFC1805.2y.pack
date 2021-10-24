#!/bin/bash

for f in \
  arpifs/adiab/cpg_gp_nhee.F90 \
  arpifs/adiab/cpg_gp_nhqe.F90 \
  arpifs/adiab/cpg_gp_hyd.F90
do
  ./pass.pl src/local/$f
  \mv src/local/$f.new src/local/$f
  rm src/local/$f.xml
done



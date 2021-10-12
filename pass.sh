#!/bin/bash

for f in \
   arpifs/phys_dmn/mf_phys.F90    \
   arpifs/adiab/cpg.F90           \
   arpifs/adiab/cpg_gp.F90        \
   arpifs/adiab/cpg_dyn.F90       \
   arpifs/adiab/cpg_dia.F90       \
   arpifs/dia/cpdysldia.F90       \
   arpifs/phys_dmn/aplpar.F90     \
   arpifs/phys_dmn/apl_arome.F90  \
   arpifs/phys_dmn/mf_phys_prep.F90 
do
  ./pass.pl src/local/$f
  \mv src/local/$f.new src/local/$f
  rm src/local/$f.xml
done



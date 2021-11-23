#!/bin/bash

for f in \
arpifs/phys_dmn/mf_phys_nhqe_part1.F90 \
arpifs/phys_dmn/mf_phys_nhqe_part2.F90 \
arpifs/phys_dmn/mf_phys_save_phsurf_part1.F90 \
arpifs/phys_dmn/mf_phys_save_phsurf_part2.F90 \
arpifs/phys_dmn/mf_phys_transfer.F90 \
arpifs/phys_dmn/mf_phys.F90
do
  echo "==> $f <=="
  ./scripts/pass.pl src/local/$f
# \mv src/local/$f.new src/local/$f
# rm src/local/$f.xml
done



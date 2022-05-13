#!/bin/bash

set -x
rm -f save_*.F90 *.xml

mkdir -p tmp mod

for f in type_model.F90 yomtoph.F90 yomphy.F90 yomphy2.F90 \
         yomphy0.F90 model_physics_mf_mod.F90 yomphy.F90 yomphy1.F90 yomphy3.F90 yomphyds.F90 \
         yomcvmnh.F90 yomvdoz.F90 yomsimphl.F90 yomarphy.F90  yomparar.F90 yommse.F90 yomlouis.F90 \
         yomnorgwd.F90 eint_mod.F90 yomdprecips.F90 yomdvisi.F90 yomsta.F90 yomcst.F90

do
  for view in $(cat .gmkview)
  do
    if [ -f "src/$view/arpifs/module/$f" ]
    then
    cp src/$view/arpifs/module/$f tmp/$f
    ./scripts/loadsave.pl --dir mod --types tmp/$f
    break
    fi
  done
done


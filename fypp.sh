#!/bin/bash

set -x
set -e

cd fypp

declare -A dir

dir[gprcp_expl]=arpifs/adiab
dir[test_ydvars]=arpifs/adiab
dir[gpmpfc_expl_part2]=arpifs/adiab
dir[gpinislb_part2_expl]=arpifs/adiab
dir[gptf2_expl_3tl_part2]=arpifs/adiab
dir[gpnspng_exp_expl]=arpifs/adiab
dir[lattex_expl_2tl]=arpifs/adiab
dir[lattex_expl_3tl]=arpifs/adiab
dir[lattex_expl_vspltrans]=arpifs/adiab
dir[lavabo_expl_lphy]=arpifs/adiab
dir[gprcp_expl]=arpifs/adiab
dir[test_ydvars]=arpifs/adiab
dir[gpmpfc_expl_part2]=arpifs/adiab
dir[gpinislb_part2_expl]=arpifs/adiab
dir[gptf2_expl_3tl_part2]=arpifs/adiab
dir[gpnspng_exp_expl]=arpifs/adiab
dir[lattex_expl_2tl]=arpifs/adiab
dir[lattex_expl_3tl]=arpifs/adiab
dir[lattex_expl_vspltrans]=arpifs/adiab
dir[lavabo_expl_lphy]=arpifs/adiab
dir[lavabo_expl_laitvspcqm_part1]=arpifs/adiab
dir[lavabo_expl_laitvspcqm_part2]=arpifs/adiab

for f in *.fypp
do

  if [ $f = "cpg_macros.fypp" ]
  then
    continue
  fi

  if [ $f = "field_definition.fypp" ]
  then
    continue
  fi

  b=$(basename $f .fypp)

  if [ "x${dir[$b]}" != "x" ]
  then
    d=${dir[$b]}
  else
    d=arpifs/module
  fi

  /opt/softs/anaconda3/bin/fypp -m os -M . -m yaml -m field_config ./$b.fypp ./$b.F90

  if [ ! -f ../src/local/$d/$b.F90 ]
  then
    cp $b.F90 ../src/local/$d/$b.F90
  else
    set +e
    cmp $b.F90 ../src/local/$d/$b.F90
    c=$?
    set -e
    if [ $c -ne 0 ]
    then
      cp $b.F90 ../src/local/$d/$b.F90
    fi
  fi

done



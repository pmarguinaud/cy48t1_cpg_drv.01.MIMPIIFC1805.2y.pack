#!/bin/bash

set -x
set -e

cd fypp

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

  d=arpifs/module

  if [ "$b" = "gprcp_expl" ]
  then
    d=arpifs/adiab
  fi

  if [ "$b" = "test_ydvars" ]
  then
    d=arpifs/adiab
  fi

  if [ "$b" = "gpmpfc_expl_part2" ]
  then
    d=arpifs/adiab
  fi

  if [ "$b" = "gpinislb_part2_expl" ]
  then
    d=arpifs/adiab
  fi

  if [ "$b" = "gptf2_expl_3tl_part2" ]
  then
    d=arpifs/adiab
  fi

  if [ "$b" = "gpnspng_exp_expl" ]
  then
    d=arpifs/adiab
  fi

  if [ "$b" = "lattex_expl_2tl" ]
  then
    d=arpifs/adiab
  fi

  if [ "$b" = "lattex_expl_3tl" ]
  then
    d=arpifs/adiab
  fi

  if [ "$b" = "lattex_expl_vspltrans" ]
  then
    d=arpifs/adiab
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



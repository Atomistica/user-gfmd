#!/bin/bash

foldername='01_test'
lmp_exe="$@"

mkdir -p $foldername

for dstep in '0.0005' '0.005' '0.05'
do
  echo $dstep
  $lmp_exe -in in.gf -var dstep \-$dstep 
  sleep 1
  $lmp_exe -in in.full -var dstep \-$dstep 
  subfolder=dstep$dstep 
  mkdir -p $foldername/$subfolder 
  cp -p in.gf in.full $foldername/$subfolder 
  mv log.gf log.full dumpinitprobegf dumpinitprobefull dumpprobegf dumpprobefull $foldername/$subfolder
done


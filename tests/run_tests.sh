#! /bin/bash

touch ../../../tools/python/__init__.py
touch ../../../tools/python/pizza/__init__.py
export PYTHONPATH="$PYTHONPATH:$(readlink -f ../../../tools/python)"

if [ -z "$1" ]; then
  echo "Syntax: run_tests.sh <path-to-LAMMPS-executable> [-n<number-of-processes>]"
  exit 999
fi

if [ ! -e "$1" ]; then
  echo "Executable $1 does not exist."
  exit 1
fi

CMD=$(readlink -f $1)

np=""
OPTIND=2
while getopts ":n:" opt; do
  case $opt in
    n)
      np=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument."
      exit 1
      ;;
  esac
done

nok=0
nfailed=0

if [ -n "$np" ]; then
  CMD="mpirun -np $np $CMD"
fi

for i in TEST_*; do

  echo $i

  cd $i

  if [ -e lammps.in ]; then
    $CMD -in lammps.in > OUT
  else
    sh run_test.sh $CMD > OUT
  fi
  if [ $? -eq 0 ]; then
    python eval.py
    if [ $? -eq 0 ]; then
      echo ".ok."
      let nok=$nok+1
    else
      echo "     .failed."
      let nfailed=$nfailed+1
    fi
  else
    echo "     .failed."
    let nfailed=$nfailed+1
  fi

  cd ..

done

let ntot=$nok+$nfailed

echo "Ran $ntot tests; $nfailed failures, $nok successes."
if [ $nfailed -gt 0 ]; then
  echo "WARNING: Some tests failed."
fi

Running the tests
-----------------

Execute

run_test.sh <LAMMPS-exectuable>

to run all tests.  

Required: 
  Python and module numpy
  LAMMPS folder ../../../tools/python 
  Ensure LAMMPS not missing the empty file ../../../tools/python/pizza/__init__.py
  
May have to give execute permissions: 
  chmod run_test.sh u+x

Old: Delete the following:
For pizza python, may have to set, from the /USER-GFMD/Tests directory: 
  bash or variant: 
    append to existing variable:
      PYTHONPATH=${PYTHONPATH}:`pwd`/../../../tools/python ; export PYTHONPATH
    create:
      PYTHONPATH=`pwd`/../../../tools/python ; export PYTHONPATH
  csh or variant:
    append to existing variable:
      setenv PYTHONPATH ${PYTHONPATH}:`pwd`/../../../tools/python
    create:
      setenv PYTHONPATH `pwd`/../../../tools/python
  

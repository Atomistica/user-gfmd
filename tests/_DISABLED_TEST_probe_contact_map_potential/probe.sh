#! /bin/sh

../../../../lmp_ser_debug_cygwin -in probe.in | grep -A 1 Step | grep -v Step | sed '/--/d' > pot.out
python ../der.py

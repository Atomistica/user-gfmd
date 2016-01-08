#! /bin/sh

tar -cvzf USER-GFMD-r`cat USER-GFMD/REV | awk '{ print $5 }'`.tar.gz --exclude .svn --exclude "*~" --exclude "*.o" --exclude "log.lammps*" --exclude "*.out" --exclude restart* --exclude OUT USER-GFMD


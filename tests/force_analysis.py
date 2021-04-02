# ======================================================================
# USER-GFMD - Green's function molecular dynamics for LAMMPS
# https://github.com/Atomistica/user-gfmd
#
# Copyright (2011-2016) Lars Pastewka <lars.pastewka@kit.edu>, Tristan A. Sharp
# and others. See the AUTHORS file in the top-level USER-GFMD directory.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ======================================================================

# Quick script to plot the force difference between GF and full dump files.
# A dump file exists for each "probe displacement magnitude"
# Tristan Feb 14 2014

# GFMD is 
# Order u**2 accurate in energy
# Order u accurate in force

# So if GFMD is working, error in force should be quadratic with displacement u

# To match up probe atoms in the two probe dump files based on location, use this algorithm:
# Travervse atom1 array: for jj=0...N-1      * * * * * * *  and find closest atom2 (at location "here")
# put the "here" values into sortorder  =   [0,1,3,3,2,4,X]
# Look for duplicate pairing and record where no pairing (X) was found in recallbadindices

import numpy
import glob

# A function to find the magnitude of the differce in two 3D vectors
def diffmag(fx1,fy1,fz1,fx2,fy2,fz2):
  return ((fx1-fx2)**(2.)+(fy1-fy2)**(2.)+(fz1-fz2)**(2.))**(0.5)

# These are the relative step sizes of probe atoms in the test
#dstep = (1e-6, 1e-4, 1e-2, 1e-0)
#folders = ('dstep0.000001/','dstep0.0001/','dstep0.01/','dstep1.0/')
folders = glob.glob('dstep*')
dstep = [float(x[5:]) for x in folders]
print(folders, dstep)
err=numpy.zeros((len(dstep),))
absmag=numpy.zeros((len(dstep),))

print("Now we will print the contents of the probe dump files, with each")
print("probe atom identified in the GF run with the corresponding one in the full run")
  
for i in numpy.arange(len(dstep)):

  print(folders[i])
  
  # read in the two files to compare
  pos_and_forces_gf   = numpy.loadtxt(folders[i]+'/dumpprobegf', skiprows=9)
  pos_and_forces_full = numpy.loadtxt(folders[i]+'/dumpprobefull', skiprows=9)
  
  x1 =pos_and_forces_gf[:,0]
  y1 =pos_and_forces_gf[:,1]
  z1 =pos_and_forces_gf[:,2]
  fx1=pos_and_forces_gf[:,3]
  fy1=pos_and_forces_gf[:,4]
  fz1=pos_and_forces_gf[:,5]
  
  x2 =pos_and_forces_full[:,0]
  y2 =pos_and_forces_full[:,1]
  z2 =pos_and_forces_full[:,2]
  fx2=pos_and_forces_full[:,3]
  fy2=pos_and_forces_full[:,4]
  fz2=pos_and_forces_full[:,5]
  
  if len(x1) != len(x2):
    print("File data not even the same length.  (Diff num probe atoms?)")
  
  # Approx distance scale between atoms for error checking
  nn_scale = 1.0 
  
  # The tests can be offset in the z-direction
  z1-=min(z1)
  z2-=min(z2)
  
  sortorder=numpy.zeros((len(x1),), dtype=numpy.int)
  recallbadindices=numpy.zeros((len(x1),), dtype=numpy.int)
  
  # Loop through all GF atoms
  no_err=1
  for jj in numpy.arange(len(x1)):
    
    # Find the distance of all other atoms to this GF atom
    distsq_of_atoms_to_atomjj = (x2-x1[jj])**(2.) + (y2-y1[jj])**(2.) + (z2-z1[jj])**(2.)
    min_distsq_of_atoms_to_atomjj = min(distsq_of_atoms_to_atomjj)
    
    # Print error message if atoms in the two files aren't even close together
    if min_distsq_of_atoms_to_atomjj > (nn_scale/4.)**(2.):
      if no_err:
        print("Warning: atoms not even close to co-located; ID-ing may break down.  Posible due to PBCs.")
      no_err = 1
      recallbadindices[jj] = 1 # Atom1[jj] has nothing near it and should never be used
    
    # Record this atom's spatially colocated atom
    else:
      here = numpy.nonzero(distsq_of_atoms_to_atomjj == min(distsq_of_atoms_to_atomjj))[0]
      if len(here) > 1:
        print("Error: Possibly multiple atoms in file2 equally near atom jj in file 1")
        exit
      sortorder[jj]=here
  
  # Sort the second file's array to the same order as the first's
  x2 = x2[sortorder]
  y2 = y2[sortorder]
  z2 = z2[sortorder]
  fx2=fx2[sortorder]
  fy2=fy2[sortorder]
  fz2=fz2[sortorder]
  
  # The force vector (and error) are 3N dimensional.  Calculate the magnitude.
  for jj in numpy.arange(len(x1)):
    if not (recallbadindices[jj]):
      err[i]    += diffmag(fx1[jj],fy1[jj],fz1[jj],fx2[jj],fy2[jj],fz2[jj])
      absmag[i] += diffmag(fx1[jj],fy1[jj],fz1[jj],0,0,0)
  
  # Print for visual check; useful if need to debug this script
  print("        ", "     x", "   y", "  z","  fx","  fy","  fz")
  for jj in numpy.arange(len(x1)):
    if not (recallbadindices[jj]):
      if (diffmag(x1[jj],y1[jj],z1[jj],x2[jj],y2[jj],z2[jj]) > (nn_scale/4.)):
        print("Fatal error: The comparison is invalid but was not flagged!!")
      print("full", x1[jj], y1[jj], z1[jj], fx1[jj], fy1[jj], fz1[jj])
      print("gf  ", x2[jj], y2[jj], z2[jj], fx2[jj], fy2[jj], fz2[jj])
  
  print(err[i], "/", absmag[i])
    
relerr = err/absmag
print("relative error", relerr)

numpy.savetxt('err.out', numpy.transpose([dstep, err, relerr]))

testpass=1
if (relerr[2] > 0.1):
  print("Greater than 10% error found when displacing probes 1e-2")
  testpass=0
if (relerr[1] > 0.1):
  print("Greater than 10% error found when displacing probes 1e-4")
  testpass=0

if ((err[2] / err[3]) > (3e-4)):
  print("dstep shrunk by 2 O.o.M.  Force error should shrink by 4 O.o.M. ")
  print("...But shrunk only by ", err[2]/err[3])
  testpass=0
if ((err[1] / err[2]) > (3e-4)):
  print("dstep shrunk by 2 O.o.M.  Force error should shrink by 4 O.o.M. ")
  print("...But shrunk only by ", err[1]/err[2])
  testpass=0

if testpass:
  print("Seems to testpass!  Check criterion v00")
  print("return 1")
 
 
  
  #cd, '..'
  #end
  
  #; plot the error in the gfmd data(black)
  #window, 0, xsize=500, ysize=500
  #plot, u, err, psym=4, thick=4, background='ffffff'x, color=0, /xlog, /ylog,$
  #title="Force Response (and Error) vs Displacement",charsize=1.3,$
  #yrange=[min(err)/10d,max(absmag)]
  #
  #print, "The grey line is a linear guideline, red is quadratic."
  #print, "Restoring force should be approximately linear, then noisy at large displacements"
  #print, "If GFMD is working, there should be a large quadratic range of the error (large data points)"
  #print, "Or check agains subversioned .png files of output"
  #; overlay a quadratic guide on error and linear guide on magnitdue
  #modelx =u;findgen(1000)/1000. * max(u)
  #modely =modelx*modelx / (max(modelx * modelx)) * max(err)
  #oplot, color='ff'x, modelx, modely
  ##modely = modelx /max(modelx) * max(absmag)
  #oplot, color='777777'x, modelx, modely
  #
  #; what percentage error is that?
  #oplot, color='555555'x, u, absmag, psym=4, thick=2
  #
  #; print also the convergence criterion; should be force tolerance
  #for i = 0, n_elements(a)-1 do $
  #  print, "Convergence criterion: ",  find_str_at_end_of_file('criterion', a[i]+'/log.full')
  #
  #stop
  #
                                                                                                                                                                                #end

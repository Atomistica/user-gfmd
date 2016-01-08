
# Quick script to plot the force difference between GF and full dump files.
# A dump file exists for each "probe displacement magnitude"

# TAS Feb 14 2014 - IDL to Python and boolean test criterion v00
# TAS May 12 2014 - Update comments and update test criterion to v01

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
import pdb 


# A function to find the magnitude of the difference of two 3D vectors
def diffmag(fx1,fy1,fz1,fx2,fy2,fz2):
  return ((fx1-fx2)**(2.)+(fy1-fy2)**(2.)+(fz1-fz2)**(2.))**(0.5)

# Read energies from LAMMPS log file
def read_energy_difference(fn):
  e = []
  f = open(fn)
  l = f.readline()
  last_l = l
  while l:
    if l.startswith('Loop time'):
      step, energy = map(float, last_l.split())
      e += [energy]
    last_l = l
    l = f.readline()
  if len(e) == 2:
    return e[1]-e[0]
  elif len(e) == 3:
    return e[2]-e[1]
  elif len(e) == 4:
    return e[3]-e[1]
  else:
    raise RuntimeError('Expected two, three or four runs in output.')


# Approx distance scale between atoms for error checking in this test
nn_scale = 2.35212725006408
  
# These are the relative step sizes of probe atoms in the test
#dstep = (1e-6, 1e-4, 1e-2, 1e-0)
#folders = ('dstep0.000001/','dstep0.0001/','dstep0.01/','dstep1.0/')
folders = glob.glob('01_test/dstep*')
dsteplist = map(lambda x: float(x[x.find('dstep')+5:]), folders)
dstep = numpy.array(dsteplist)

interactivemode = 0
if interactivemode:
  print folders, dstep

# ---------------------------------------------- # 
# Nothing below here should need to be changed
# ---------------------------------------------- # 

# The force on probe atoms (squared) and 
# the error in the GF-calculated force (squared) at that dstep
err2=numpy.zeros((len(dstep),))
absmag2=numpy.zeros((len(dstep),))
eerr=numpy.zeros((len(dstep),))

if interactivemode:
  print "Now we will print the contents of the probe dump files, with each"
  print "probe atom identified in the GF run with the corresponding one in the full run"
  
for i in numpy.arange(len(dstep)):

  if interactivemode:
    print folders[i]

  # read the energies
  denergy_gf = read_energy_difference(folders[i]+'/log.gf')
  denergy_full = read_energy_difference(folders[i]+'/log.full')
  eerr[i] = denergy_full-denergy_gf
  
  # read in the init forces
  init_pos_and_forces_gf   = numpy.loadtxt(folders[i]+'/dumpinitprobegf', skiprows=9).reshape(-1,6)
  init_pos_and_forces_full = numpy.loadtxt(folders[i]+'/dumpinitprobefull', skiprows=9).reshape(-1,6)

  # read in the final forces
  pos_and_forces_gf   = numpy.loadtxt(folders[i]+'/dumpprobegf', skiprows=9).reshape(-1,6)
  pos_and_forces_full = numpy.loadtxt(folders[i]+'/dumpprobefull', skiprows=9).reshape(-1,6)

  x1 =init_pos_and_forces_gf[:,0]
  y1 =init_pos_and_forces_gf[:,1]
  z1 =init_pos_and_forces_gf[:,2]
  dx1=pos_and_forces_gf[:,0]-init_pos_and_forces_gf[:,0]
  dy1=pos_and_forces_gf[:,1]-init_pos_and_forces_gf[:,1]
  dz1=pos_and_forces_gf[:,2]-init_pos_and_forces_gf[:,2]
  fx1=pos_and_forces_gf[:,3]-init_pos_and_forces_gf[:,3]
  fy1=pos_and_forces_gf[:,4]-init_pos_and_forces_gf[:,4]
  fz1=pos_and_forces_gf[:,5]-init_pos_and_forces_gf[:,5]
 
  x2 =init_pos_and_forces_full[:,0]
  y2 =init_pos_and_forces_full[:,1]
  z2 =init_pos_and_forces_full[:,2]
  dx2=pos_and_forces_full[:,0]-init_pos_and_forces_full[:,0]
  dy2=pos_and_forces_full[:,1]-init_pos_and_forces_full[:,1]
  dz2=pos_and_forces_full[:,2]-init_pos_and_forces_full[:,2]
  fx2=pos_and_forces_full[:,3]-init_pos_and_forces_full[:,3]
  fy2=pos_and_forces_full[:,4]-init_pos_and_forces_full[:,4]
  fz2=pos_and_forces_full[:,5]-init_pos_and_forces_full[:,5]

  if len(x1) != len(x2):
    print "File data not even the same length!  Invalid test."\
    "(Different number of probe atoms in in.gf and in.full?)"
  
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
        print "Warning: atoms not even close to co-located; ID-ing may break down and "\
        "artifically increase the error.  May be due to PBCs, if only 1 atom has wrapped."
      no_err = 1
      recallbadindices[jj] = 1 # Atom1[jj] has nothing near it and should never be used
    
    # Record this atom's spatially colocated atom
    else:
      here = numpy.nonzero(distsq_of_atoms_to_atomjj == min(distsq_of_atoms_to_atomjj))[0]
      if len(here) > 1:
        print "Error: Possibly multiple atoms in file2 equally near atom jj in file 1"
        exit
      sortorder[jj]=here
  
  # Sort the second file's array to the same order as the first's
  x2 = x2[sortorder]
  y2 = y2[sortorder]
  z2 = z2[sortorder]
  fx2=fx2[sortorder]
  fy2=fy2[sortorder]
  fz2=fz2[sortorder]
  
  # The force vector (and error) are 3N dimensional.  Calculate the magnitude squared.
  for jj in numpy.arange(len(x1)):
    if not (recallbadindices[jj]):
      err2[i]    += (diffmag(fx1[jj],fy1[jj],fz1[jj],fx2[jj],fy2[jj],fz2[jj]))**(2.0)
      absmag2[i] += (diffmag(fx1[jj],fy1[jj],fz1[jj],0,0,0))**(2.0)
  
  # Print for visual check; useful if need to debug this script
  if interactivemode:
    print "%14s%14s%14s%14s%14s%14s%14s%14s" % ( "        ", "    de", "    dx", "   dy", "  dz", "  fx", "  fy", "  fz" )
    for jj in numpy.arange(len(dx1)):
      if not (recallbadindices[jj]):
        if (diffmag(x1[jj],y1[jj],z1[jj],x2[jj],y2[jj],z2[jj]) > (nn_scale/4.)):
          print "Analysis script fatal error: The comparison is invalid but was not flagged!"
        print "%14s%14e%14e%14e%14e%14e%14e%14e" % ( "full", denergy_full, dx2[jj], dy2[jj], dz2[jj], fx2[jj], fy2[jj], fz2[jj] )
        print "%14s%14e%14e%14e%14e%14e%14e%14e" % ( "gf  ", denergy_gf, dx1[jj], dy1[jj], dz1[jj], fx1[jj], fy1[jj], fz1[jj] )
        print "%14s%14e%14e%14e%14e%14e%14e%14e" % ( "diff", denergy_full-denergy_gf, x1[jj]-x2[jj], y1[jj]-y2[jj], z1[jj]-z2[jj], fx1[jj]-fx2[jj], fy1[jj]-fy2[jj], fz1[jj]-fz2[jj] )
    print (err2[i])**(0.5), "/", (absmag2[i])**(0.5)
    
relerr = (err2)**(0.5)/(absmag2)**(0.5)
if interactivemode:
 print "relative error", relerr

numpy.savetxt('err.out', numpy.transpose([dstep, (err2)**(0.5), relerr, eerr]))


# Create an automated test to return a pass/fail boolean
# Check that force error goes like the square of the displacement err ~ u^(2-tolerance)
testpass=1

# Pick two folders (large dstep and small dstep)
largedstepindex = numpy.flatnonzero((dstep >= 1e-2*nn_scale) & (dstep < 5e-2*nn_scale))
if len(largedstepindex) < 1:
  print "No test dstep was between 1e-2 and 5e-2 of the nearest neighbor distance", nn_scale
  testpass=0
else:
  largedstepindex = (largedstepindex)[0]
  if (relerr[largedstepindex] > 0.1):
    raise RuntimeError("Greater than 10% error found when displacing probes "+repr(dstep[largedstepindex]))
    testpass=0

smalldstepindex = numpy.flatnonzero((dstep >= 1e-4*nn_scale) & (dstep < 5e-3*nn_scale))
if len(smalldstepindex) < 1:
  print "No test dstep was between 1e-4 and 5e-3 of the nearest neighbor distance", nn_scale
  testpass=0
else:
  smalldstepindex = (smalldstepindex)[0]

if testpass == 1:
  tolerance = 0.2
  largedsteperr=(err2[largedstepindex])**(0.5)
  smalldsteperr=(err2[smalldstepindex])**(0.5)
  goeslike = numpy.log(smalldsteperr / largedsteperr) / numpy.log(dstep[smalldstepindex] / dstep[largedstepindex]) 
  print "Difference goes like displacement to the "+repr(goeslike)+\
  " between u= "+repr(dstep[smalldstepindex])+" and "+repr(dstep[largedstepindex])
  if (2.0 - goeslike > tolerance):
    raise RuntimeError("Fail: Force error did not shrink (approx) like the square of displacement.")
    testpass=0

if testpass:
  print "Test passed!  Check criterion: v01"
  print "return 1"

f=open('goodness_measure.txt', 'w')
f.write('Difference goes like displacement to the '+repr(goeslike)+'\n')
f.write('  between u= '+repr(dstep[smalldstepindex])+' and '+repr(dstep[largedstepindex])+'\n')
f.write('Fractional error '+repr(relerr[largedstepindex])+' at dstep '+repr(dstep[largedstepindex])+'\n')
f.close()

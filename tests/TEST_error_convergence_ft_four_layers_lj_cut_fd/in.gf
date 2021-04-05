# Check for quadratic convergence of difference between GF and all-atom
# This is the GF run script

log log.gf
shell "rm log.lammps"

units           lj
boundary        p p f
processors      * * 1
atom_style      gfmd

# Create structure
variable        height equal 20
variable        bottom equal 0
include         in.structure

# Setup interaction
variable        rcut equal 1.75 # between 2nd and 3rd
pair_style      lj/cut/gf ${rcut}
pair_coeff      * * 1.0 1.0 ${rcut}
pair_modify	shift yes

# Optional neighbor command
neighbor        1.0 bin
neigh_modify    every 1 delay 0 check yes

# Screen output
thermo          50
thermo_style    custom step pe
thermo_modify   norm no format float "%50.45e"

# Debug dump output
#dump            debug all nc 10 traj_gf.nc type x y z fx fy fz

# Bottom row of atoms is GFMD layer
fix             gfmd bottom gfmd gfmd ft fcc100 ${nndist} 4 lj/cut/fd height ${height} reset_xeq reset_map
fix_modify      gfmd energy yes

# Initialize (before pushing) allowing to equilibriate before holding any atoms rigid
min_style       cg
min_modify      dmax 0.01 line quadratic
minimize        0.0 1.0e-7 1000 1000

# Store the initial probe atom locations
dump            initprobedump probe1 custom 1 dumpinitprobegf x y z fx fy fz
dump_modify     initprobedump format line "%f %f %f %50.45e %50.45e %50.45e"
run             0 
undump          initprobedump

# Make rigid after surface relaxes
fix             2      probe1 setforce 0.0 0.0 0.0
fix_modify      2      energy no
fix             fixID1 probe1 aveforce 0.0 0.0 0.0
fix_modify      fixID1 energy no

# Image Dump (optional)
#undump         3
#dump           3 all image 10 aimag.gf.*.${n}.jpg fz type adiam 0.8 size 1024 1024 view 90.0 0.0 center s 0.5 0.5 0.5 zoom 2 shiny 1 axes yes 0.8 0.02 
#dump_modify    3 backcolor white amap -1.0 1.0 ca 0 5 min red -.0000001 gray 0.0 black 0.0000001 navy max blue
#dump_modify    3 pad 4

# Displace probe group and minimize
displace_atoms  probe1 move ${dstep} 0 ${dstep} units box
min_style       cg
min_modify      dmax 0.3 line quadratic
minimize        0.0 1.0e-7 100000 100000

# Store the final probe atom locations and forces
unfix           2
unfix           fixID1
dump            probedump probe1 custom 1 dumpprobegf x y z fx fy fz
dump_modify     probedump format line "%f %f %f %50.45e %50.45e %50.45e"
run             0



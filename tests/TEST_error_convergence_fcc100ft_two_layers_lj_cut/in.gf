# Check for quadratic convergence of difference between GF and all-atom
# This is the GF run script

log log.gf
shell "rm log.lammps"

units		lj
boundary        p p f
processors      * * 1
atom_style	gfmd
newton          off

read_data       8gfdatafilelj_12e_hiprec.lammps

variable rcut equal 1.75 # between 2nd and 3rd

pair_style lj/cut ${rcut}
pair_coeff  * * 1.0 1.0 ${rcut}
pair_coeff  * 3 1.0 1.0 ${rcut}
pair_coeff  3 * 1.0 1.0 ${rcut}
pair_coeff  2 2 0.0 1.0 0.1
pair_coeff  3 3 1.0 1.0 ${rcut}
pair_modify	shift yes

# Optional neighbor command
neighbor        1.0 bin
neigh_modify    every 1 delay 0 check yes
group           rbgroup type 3
group           gfmdgroup type 2
group           type4 type 4
region          columnregion1 block 0.5 1.6 0.5 1.612 INF INF units box
region          columnregion2 block 1.6 2.7 4.5 5.612 INF INF units box
group           col1group region columnregion1
group           col2group region columnregion2
group           columns union col1group col2group
group           probe1 intersect columns gfmdgroup
#group           probe1 type 4

# GFMD retains linear response; optionally compare to making GF atoms rigid
# fix             bottomfix gfmdgroup setforce 0.0 0.0 0.0

thermo		50


fix             gfmd gfmdgroup gfmd gfmd fcc100ft x 2 68.527040 -0.55454572 -3.3768452 0.55455513 0.0 0.0 height 20 linf 0.87256787 -0.87256787 mulphase 0.5 # 1.7451357 -1.7451357 # TAS2
#fix             gfmd gfmdgroup gfmd gfmd fcc100 2 68.527040 -3.3768452 height 2 linf 1.7451357 -1.7451357
fix_modify      gfmd energy yes


# initialize (before pushing) allowing to equilibriate before holding any atoms rigid
min_style cg
min_modify dmax 0.01 line quadratic
minimize 0.0 1.0e-7 1000 1000

# run 100

# Store the initial probe atom locations
# dump  initprobedump all custom 1 dumpinitallgf x y z fx fy fz
# run 0 
# undump initprobedump

# Make rigid after surface relaxes
fix             2      probe1 setforce 0.0 0.0 0.0
fix_modify      2      energy no
fix             fixID1 probe1 aveforce 0.0 0.0 0.0
fix_modify      fixID1 energy no



#variable n loop 0 3
#label nloop

# Image Dump (optional)
#undump 3
#dump 3 all image 10 aimag.gf.*.${n}.jpg fz type adiam 0.8 size 1024 1024 view 90.0 0.0 center s 0.5 0.5 0.5 zoom 2 shiny 1 axes yes 0.8 0.02 
#dump_modify 3 backcolor white amap -1.0 1.0 ca 0 5 min red -.0000001 gray 0.0 black 0.0000001 navy max blue
#dump_modify 3 pad 4


# Several loops to encourage minimize to hit force tolerance, before pushing on the last loop
#if "$n < 3" then "variable dstep equal 0.0" else "variable dstep equal -0.0001"
displace_atoms probe1 move 0 ${dstep} ${dstep} units box
min_style cg
min_modify dmax 0.3 line quadratic
minimize 0.0 1.0e-7 100000 100000

# # Dump all atoms for posterity
# dump            dumpid all custom 1 dump.ljgf x y z fx fy fz
# run 0
# undump          dumpid

#next n
#jump SELF nloop

#run 100

unfix 2
unfix fixID1
dump           probedump probe1 custom 1 dumpprobegf x y z fx fy fz
run 0



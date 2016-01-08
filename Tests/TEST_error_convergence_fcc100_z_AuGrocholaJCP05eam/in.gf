# Check for quadratic convergence of difference between GF and all-atom
# This is the GF run script


log log.gf
shell "rm log.lammps"

units		metal
boundary        p p f
processors      * * 1
atom_style	gfmd
newton          off

read_data       8gfdatafile

# variable lj_to_metal equal 2.5702 # = (Au cell)/(LJ nn cell)= (4.08 nm)/(2^1/6 * 2^1/2 sigma)

pair_style eam/alloy
pair_coeff   * * Au-Grochola-JCP05.eam.alloy Au Au Au Au Au
pair_modify	shift yes

# Optional neighbor command
neighbor        4.0 bin
neigh_modify    every 1 delay 0 check yes
group           gfmdgroup type 2 4
group           type4 type 4
group           type5 type 5

# Create a probe of atoms to displace
region          columnregion1 block 0.5 4.6 0.5 4.6 INF INF units box
region          columnregion2 block 12.6 14.7 12.4 14.7 INF INF units box
group           col1group region columnregion1
group           col2group region columnregion2
group           columns union col1group col2group
#group           probe intersect columns type4
#group           probe type 2 4 # full 4 planes of gfmd layer
group		probe type 4 # top 2 planes of the gfmd layer
#group		probe type 3 # above gfmd layer
group		probe1 intersect col1group type4
group		probe2 intersect col2group type4
#group		probe union probe1 probe2

# GFMD retains linear response; optionally compare to making GF atoms rigid
# fix             bottomfix gfmdgroup setforce 0.0 0.0 0.0

thermo		50

fix             mygfmdfix gfmdgroup gfmd gfmd fcc100fteam           &
                -1.2303882E-005  8.30808382204118   3               &
                -5.4140458E-002      6.680740E-002                  & 
                -3.0213523E-002     -1.989247E-002                  & 
                -9.5703736E-003      3.781924E-002                  & 
                -0.1021523       1.54040       -3.54941782708E-002  &
                 5.3959908E-002 -7.471582E-002  1.325760451E-002    &
                 1.8473237E-002 -6.431390E-002  3.705879027E-003    & 
                height 10
fix_modify      mygfmdfix energy yes


# Optional dump_image command
#dump 3 all image 1 aimag.gf.*.jpg fz type adiam 0.8 size 1024 1024 view 90.0 0.0 center s 0.5 0.5 0.5 zoom 2 shiny 1 axes yes 0.8 0.02
#dump_modify 3 backcolor white amap -0.01 0.01 ca 0 6 min red -1e-6 black -9e-7 gray 9e-7 gray 1e-6 blue max aqua
#dump_modify 3 pad 4


# initialize (before pushing) allowing to equilibriate before holding any atoms rigid
min_style cg
min_modify dmax 0.3 line quadratic
minimize 0.0 1.0e-9 1000 1000

run 1 # can help with minimization
minimize 0.0 1.0e-9 1000 1000

# Store the initial probe atom locations
# dump  initprobedump all custom 1 dumpinitallgf x y z fx fy fz
# run 0 
# undump initprobedump

# Make rigid after surface relaxes
fix             2      probe setforce 0.0 0.0 0.0
fix_modify      2      energy no
fix             fixID1 probe aveforce 0.0 0.0 0.0
fix_modify      fixID1 energy no



displace_atoms probe1 move ${dstep} 0 ${dstep} units box
displace_atoms probe2 move 0 ${dstep} 0        units box

min_modify dmax 0.3 line quadratic
minimize 0.0 1.0e-9 100000 100000
run 1
minimize 0.0 1.0e-9 100000 100000


# # Dump all atoms for posterity
# dump            dumpid all custom 1 dump.ljgf x y z fx fy fz
# run 0
# undump          dumpid

#run 100

unfix 2
unfix fixID1
dump           probedump probe custom 1 dumpprobegf x y z fx fy fz
run 0



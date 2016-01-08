
units		metal
boundary        p p f
processors      1 1 1
atom_style	gfmd
newton          off

read_data       2255negbare_gfdatafile
  orthogonal box = (0 0 -22.448) to (14.425 14.425 33.464)
  1 by 1 by 1 MPI processor grid
  150 atoms
#read_data       ../rb_fcc_h05brms.25restor_l1024x8_lmin11.8_eb_sq_flat_5typed

variable lj_to_metal equal 2.5702 # = (Au cell)/(LJ nn cell)= (4.08 nm)/(2^1/6 * 2^1/2 sigma)
variable rcut_interface equal 1.1224620483*${lj_to_metal}
variable rcut_interface equal 1.1224620483*2.5701999999999998181

variable rcut_expand equal ${rcut_interface}+4.0
variable rcut_expand equal 2.8849519565406600563+4.0

pair_style hybrid eam/alloy lj/expand ${rcut_interface}
pair_style hybrid eam/alloy lj/expand 2.8849519565406600563
pair_coeff   * * eam/alloy Au-Grochola-JCP05.eam.alloy Au Au NULL Au Au

#pair_coeff   * * lj/cut    0.0 ${lj_to_metal} ${rcut_interface}
#pair_coeff   * * eam/alloy Au-Grochola-JCP05.eam.alloy NULL NULL NULL Au NULL

pair_coeff   1 3 lj/expand    0.0 ${lj_to_metal} 0.0 ${rcut_interface}
pair_coeff   1 3 lj/expand    0.0 2.5701999999999998181 0.0 ${rcut_interface}
pair_coeff   1 3 lj/expand    0.0 2.5701999999999998181 0.0 2.8849519565406600563
pair_coeff   2 3 lj/expand    0.0 ${lj_to_metal} 0.0 ${rcut_interface}
pair_coeff   2 3 lj/expand    0.0 2.5701999999999998181 0.0 ${rcut_interface}
pair_coeff   2 3 lj/expand    0.0 2.5701999999999998181 0.0 2.8849519565406600563
pair_coeff   3 3 lj/expand    0.0 ${lj_to_metal} 0.0 ${rcut_interface}
pair_coeff   3 3 lj/expand    0.0 2.5701999999999998181 0.0 ${rcut_interface}
pair_coeff   3 3 lj/expand    0.0 2.5701999999999998181 0.0 2.8849519565406600563
pair_coeff   3 4 lj/expand    0.0 ${lj_to_metal} 4.0 ${rcut_interface} #expand to decrease corrugation # disabled
pair_coeff   3 4 lj/expand    0.0 2.5701999999999998181 4.0 ${rcut_interface} 
pair_coeff   3 4 lj/expand    0.0 2.5701999999999998181 4.0 2.8849519565406600563 
pair_coeff   3 5 lj/expand    0.0 ${lj_to_metal} 0.0 ${rcut_interface}
pair_coeff   3 5 lj/expand    0.0 2.5701999999999998181 0.0 ${rcut_interface}
pair_coeff   3 5 lj/expand    0.0 2.5701999999999998181 0.0 2.8849519565406600563
pair_modify	shift yes


neighbor        4.0 bin
neigh_modify    every 1 delay 0 check yes
group           rbgroup type 3
50 atoms in group rbgroup
group           gfmdgroup type 2 5
100 atoms in group gfmdgroup
group           type2 type 2
50 atoms in group type2
group           type4 type 4
0 atoms in group type4
group           type5 type 5
50 atoms in group type5
#variable        maxgfht equal bound(type5,zmax)
#region          topgfslice block INF INF INF INF ${maxgfht}-0.1 ${maxgfht}+0.1 units box
#group           surfgfatoms region topgfslice
region          columnregion1 block 0.5 4.6 0.5 4.6 INF INF units box
region          columnregion2 block 2.6 6.7 4.6 6.7 INF INF units box
group           col1group region columnregion1
15 atoms in group col1group
group           col2group region columnregion2
6 atoms in group col2group
group           columns union col1group col2group
21 atoms in group columns
group           probe intersect columns type5 # surfgfatoms
7 atoms in group probe


#fix             bottomfix gfmdgroup setforce 0.0 0.0 0.0

thermo		500



#fix myfixdel type2 addforce 0.0 0.0 -0.01 # del this again


# Idea to have at rho of energy minimum.
# Stiffness matrix (par term) is dFdrho * d2fdr2
#fix             mygfmdfix gfmdgroup gfmd gfmd fcc100fteam #                0.0 8.0 4 #                2 3 #                0 0 #                0 0 #                0 0 #                0 0 #                0 0 #                height 4 # dump_greens_function
#fix_modify      mygfmdfix energy yes

#fix             mygfmdfix gfmdgroup gfmd gfmd fcc100fteam #                -0.075 8.3 4 #                0.001 0.001 #                0 0 #                0 0 #                0004 0 #                000.2 0 #                000.16 0 #                height 4 # dump_greens_function
#
#
#fix_modify      mygfmdfix energy yes

# Caution about rho going negative?
#fix             mygfmdfix gfmdgroup gfmd gfmd fcc100fteam #                 0.0 8.0 4   #                 -0.1 .2    #                 -0.01 .02  #                 -.002 .004  #                 1000.4   0.0   #                 0.0   0.1   #                 0.0   0.05  #                 height 4 dump_greens_function
#fix_modify      mygfmdfix energy yes

# fix id group gfmd solver Dcalculator
#              dF/drho  d2F/drho2 num_atoms_per_gf_layer_cell
#              dfdr     d2fdr2 ( nn)
#              dfdr     d2fdr2 (2nn)
#              dfdr     d2fdr2 (3nn)
#              d2phidr2 dphidr ( nn) is this the same as k
#              d2phidr2 dphidr (2nn)
#              d2phidr2 dphidr (3nn)
#              height value

fix             mygfmdfix gfmdgroup gfmd gfmd fcc100fteam                 -7.521826747793241E-002  8.32896611616331 4                   -5.368227746962540E-002  6.422825369266792E-002                   -3.041321389379986E-002 -2.049745662446555E-002                   -9.126569474822686E-003  3.545931446285652E-002                    1.46855     -0.091631                   -.069519      0.053247                   -.038890      0.017848                   height 5  #dump_greens_function
USER-GFMD rev $Id: REV 357:358 $
Using a 0.015000 % relative thickness for the ghost shell.
Using elasticity kernel fcc100fteam with total substrate height of 5 layers and 3 degrees of freedom per layer.
Surface normal is 1.
Phase multiplier is 0.500000.
Spring constants for 1-th neighbor interaction are -0.053682 (dfdr) and 0.064228 (d2fdr2) and pairpot k 1.468550 and dphidr -0.091631.
Internal force for this layer is 0.000000.
Spring constants for 2-th neighbor interaction are -0.030413 (dfdr) and -0.020497 (d2fdr2) and pairpot k -0.069519 and dphidr 0.053247.
Internal force for this layer is 0.000000.
Spring constants for 3-th neighbor interaction are -0.009127 (dfdr) and 0.035459 (d2fdr2) and pairpot k -0.038890 and dphidr 0.017848.
Internal force for this layer is 0.000000.
Elastic grid has dimension 5 x 5 with 4 atoms per unit cell.
Elastic grid has dimension 5 x 5 per processor.
0.000000 x 0.000000 - 14.424978 x 14.424978 section of the domain is on proc 0.
0 x 0 - 4 x 4 section of the grid is on proc 0.
100 atoms found in the elastic manifold total.
Surface vector along U is ( 2.885, 0 ).
Surface vector along V is ( 0, 2.885 ).
Using static FFT solver.
fix_modify      mygfmdfix energy yes


#fix             mygfmdfix gfmdgroup gfmd gfmd fcc100fteam #                0 0 4   #                0 0 #                0 0 #                0 0 #                 1.5     -0.0916   #                -.07      0.0532   #                -.04      0.0178   #                height 5  #dump_greens_function
#fix_modify      mygfmdfix energy yes


#fix          phigfmd gfmdgroup gfmd gfmd fcc100ft 4 84.9014 -3.601 0.88829  height 4 # dump_greens_function  # shouldnt be able to redefine nu
#fix_modify   phigfmd energy yes
#fix          simplegfmd type2 gfmd gfmd fcc100ft 2 40 -2 0.5 height 6 # dump_greens_function # 84.9014 -3.601 0.88829 height 3  dump_greens_function # shouldnt be able to redefine nu
#fix_modify   simplegfmd energy yes


dump 3 all image 10 aimag.gf.*.jpg fz type adiam 0.8 size 1024 1024 view 90.0 0.0 center s 0.5 0.5 0.5 zoom 2 shiny 1 axes yes 0.8 0.02
dump_modify 3 backcolor white amap -1.0 1.0 ca 0.0 4 min red 0.0 black 0.000001 navy max blue
dump_modify 3 pad 4

# initialize for pushing
min_style cg
min_modify dmax 0.3*${lj_to_metal} line quadratic
min_modify dmax 0.3*2.5701999999999998181 line quadratic
minimize 0.0 1.0e-6 100000 100000
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Grid shift should be 0, 0. (from center of mass)
Reference zero for layer 0 is 1.442500, 1.442500.
Reference zero for layer 1 is 0.000000, 0.000000.
Reference zero for layer 2 is 1.442500, 1.442500.
Reference zero for layer 3 is 0.000000, 0.000000.
Grid shift is 0, 0. ( from grid index )
Hence, grid needs to be shifted by 0, 0 to be aligned with domain boundaries.
Memory usage per processor = 3.98023 Mbytes
Step Temp E_pair E_mol TotEng Press 
       0            0   -186.52402            0   -186.52402   -10313.406 
      59            0   -188.74495            0   -188.61625   -3933.5285 
Loop time of 0.283108 on 1 procs for 59 steps with 150 atoms

Minimization stats:
  Stopping criterion = force tolerance
  Energy initial, next-to-last, final = 
        -186.524022579     -188.616253538     -188.616253538
  Force two-norm initial, final = 4.46589 6.68076e-07
  Force max component initial, final = 0.728508 1.37517e-07
  Final line search alpha, max atom move = 1 1.37517e-07
  Iterations, force evaluations = 59 117

Pair  time (%) = 0.0339158 (11.9798)
Neigh time (%) = 0 (0)
Comm  time (%) = 0.000537157 (0.189736)
Outpt time (%) = 0.189165 (66.8172)
Other time (%) = 0.0594902 (21.0133)
Gffft  time (%) = 0.0205936 (7.27413)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11182 ave 11182 max 11182 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11182
Ave neighs/atom = 74.5467
Neighbor list builds = 0
Dangerous builds = 0
min_style hftn
minimize 0.0 1.0e-6 100000 100000
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Memory usage per processor = 4.90949 Mbytes
Step Temp E_pair E_mol TotEng Press 
      59            0   -188.74495            0   -188.61625   -3933.5285 
Loop time of 0.000582933 on 1 procs for 0 steps with 150 atoms

Minimization stats:
  Stopping criterion = force tolerance
  Energy initial, next-to-last, final = 
        -188.616253538     -188.616253538     -188.616253538
  Force two-norm initial, final = 6.68076e-07 6.68076e-07
  Force max component initial, final = 1.37517e-07 1.37517e-07
  Final line search alpha, max atom move = 0 0
  Iterations, force evaluations = 0 0

Pair  time (%) = 0.000359058 (61.5951)
Neigh time (%) = 0 (0)
Comm  time (%) = 5.00679e-06 (0.858896)
Outpt time (%) = 0 (0)
Other time (%) = 0.000218868 (37.546)
Gffft  time (%) = 0.000442982 (75.9918)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    789 ave 789 max 789 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11689 ave 11689 max 11689 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11689
Ave neighs/atom = 77.9267
Neighbor list builds = 0
Dangerous builds = 0
min_style fire
minimize 0.0 1.0e-6 100000 100000
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Memory usage per processor = 3.28329 Mbytes
Step Temp E_pair E_mol TotEng Press 
      59            0   -188.74495            0   -188.61625   -3933.5285 
      60 1.5242508e-24   -188.74495            0   -188.61625   -3933.5285 
Loop time of 0.038645 on 1 procs for 1 steps with 150 atoms

Minimization stats:
  Stopping criterion = force tolerance
  Energy initial, next-to-last, final = 
        -188.616253538     -188.616253538     -188.616253538
  Force two-norm initial, final = 6.68076e-07 6.68076e-07
  Force max component initial, final = 1.37517e-07 1.37517e-07
  Final line search alpha, max atom move = 0 0
  Iterations, force evaluations = 1 1

Pair  time (%) = 0.000684977 (1.77248)
Neigh time (%) = 0 (0)
Comm  time (%) = 8.10623e-06 (0.0209761)
Outpt time (%) = 0 (0)
Other time (%) = 0.0379519 (98.2065)
Gffft  time (%) = 0.000653982 (1.69228)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    789 ave 789 max 789 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11689 ave 11689 max 11689 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11689
Ave neighs/atom = 77.9267
Neighbor list builds = 0
Dangerous builds = 0

#variable rbminht equal bound(rbgroup,zmin)
#variable ebmaxht equal bound(type4,zmax)
#variable ebrboverlap equal ${ebmaxht}-${rbminht}+${rcut_expand}-0.3075 # 0.3075 for expand of 4.0
#displace_atoms  rbgroup move 0.01 0.0 ${ebrboverlap} units box # break symmetry


# apply these after surface relaxes
fix             2      probe setforce 0.0 0.0 0.0
fix_modify      2      energy no
fix             fixID1 probe aveforce 0.0 0.0 0.0
fix_modify      fixID1 energy no


variable n loop 0 3
label nloop

undump 3
dump 3 all image 1 aimag.gf.*.${n}.jpg fz type adiam 0.8 size 1024 1024 view 90.0 0.0 center s 0.5 0.5 0.5 zoom 2 shiny 1 axes yes 0.8 0.02
dump 3 all image 1 aimag.gf.*.0.jpg fz type adiam 0.8 size 1024 1024 view 90.0 0.0 center s 0.5 0.5 0.5 zoom 2 shiny 1 axes yes 0.8 0.02
dump_modify 3 backcolor white amap -1.0 1.0 ca 0.0 4 min red 0.0 black 0.000001 navy max blue
dump_modify 3 pad 4

if "$n < 3"  then "variable dstep equal 0.0" else "variable dstep equal -1.6"
variable dstep equal 0.0


displace_atoms probe move 0 0 ${dstep} units box
displace_atoms probe move 0 0 0 units box
min_style cg
min_modify dmax 0.3*${lj_to_metal} line quadratic
min_modify dmax 0.3*2.5701999999999998181 line quadratic
minimize 0.0 1.0e-6 100000 100000
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Memory usage per processor = 3.98023 Mbytes
Step Temp E_pair E_mol TotEng Press 
      60 1.5242508e-24   -188.74495            0   -188.61625   -3933.5285 
      61 1.5242508e-24   -188.74495            0   -188.61625   -3933.5283 
Loop time of 0.039443 on 1 procs for 1 steps with 150 atoms

Minimization stats:
  Stopping criterion = force tolerance
  Energy initial, next-to-last, final = 
        -188.616253538     -188.616253538     -188.616253538
  Force two-norm initial, final = 6.44609e-07 4.25777e-07
  Force max component initial, final = 1.37517e-07 7.60257e-08
  Final line search alpha, max atom move = 1 7.60257e-08
  Iterations, force evaluations = 1 2

Pair  time (%) = 0.000976324 (2.47528)
Neigh time (%) = 0 (0)
Comm  time (%) = 1.81198e-05 (0.0459392)
Outpt time (%) = 0 (0)
Other time (%) = 0.0384486 (97.4788)
Gffft  time (%) = 0.000869036 (2.20327)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    789 ave 789 max 789 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11689 ave 11689 max 11689 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11689
Ave neighs/atom = 77.9267
Neighbor list builds = 0
Dangerous builds = 0
min_style hftn
minimize 0.0 1.0e-6 100000 100000
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Memory usage per processor = 4.90949 Mbytes
Step Temp E_pair E_mol TotEng Press 
      61 1.5242508e-24   -188.74495            0   -188.61625   -3933.5283 
Loop time of 0.000566006 on 1 procs for 0 steps with 150 atoms

Minimization stats:
  Stopping criterion = force tolerance
  Energy initial, next-to-last, final = 
        -188.616253538     -188.616253538     -188.616253538
  Force two-norm initial, final = 4.25777e-07 4.25777e-07
  Force max component initial, final = 7.60257e-08 7.60257e-08
  Final line search alpha, max atom move = 0 0
  Iterations, force evaluations = 0 0

Pair  time (%) = 0.000355959 (62.8896)
Neigh time (%) = 0 (0)
Comm  time (%) = 5.00679e-06 (0.884583)
Outpt time (%) = 0 (0)
Other time (%) = 0.00020504 (36.2258)
Gffft  time (%) = 0.000417948 (73.8416)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    789 ave 789 max 789 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11689 ave 11689 max 11689 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11689
Ave neighs/atom = 77.9267
Neighbor list builds = 0
Dangerous builds = 0
min_style fire
minimize 0.0 1.0e-6 100000 100000
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Memory usage per processor = 3.28329 Mbytes
Step Temp E_pair E_mol TotEng Press 
      61            0   -188.74495            0   -188.61625   -3933.5283 
      62 6.1911061e-25   -188.74495            0   -188.61625   -3933.5283 
Loop time of 0.039124 on 1 procs for 1 steps with 150 atoms

Minimization stats:
  Stopping criterion = force tolerance
  Energy initial, next-to-last, final = 
        -188.616253538     -188.616253538     -188.616253538
  Force two-norm initial, final = 4.25777e-07 4.25777e-07
  Force max component initial, final = 7.60257e-08 7.60257e-08
  Final line search alpha, max atom move = 0 0
  Iterations, force evaluations = 1 1

Pair  time (%) = 0.000684023 (1.74835)
Neigh time (%) = 0 (0)
Comm  time (%) = 7.86781e-06 (0.0201099)
Outpt time (%) = 0 (0)
Other time (%) = 0.0384321 (98.2315)
Gffft  time (%) = 0.000610113 (1.55943)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    789 ave 789 max 789 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11689 ave 11689 max 11689 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11689
Ave neighs/atom = 77.9267
Neighbor list builds = 0
Dangerous builds = 0


#dump            dumpid type4 custom 1 dump.melt.*.gf x y z fx fy fz
#run 0
#undump          dumpid

next n
jump SELF nloop

undump 3
dump 3 all image 1 aimag.gf.*.${n}.jpg fz type adiam 0.8 size 1024 1024 view 90.0 0.0 center s 0.5 0.5 0.5 zoom 2 shiny 1 axes yes 0.8 0.02
dump 3 all image 1 aimag.gf.*.1.jpg fz type adiam 0.8 size 1024 1024 view 90.0 0.0 center s 0.5 0.5 0.5 zoom 2 shiny 1 axes yes 0.8 0.02
dump_modify 3 backcolor white amap -1.0 1.0 ca 0.0 4 min red 0.0 black 0.000001 navy max blue
dump_modify 3 pad 4

if "$n < 3"  then "variable dstep equal 0.0" else "variable dstep equal -1.6"
variable dstep equal 0.0


displace_atoms probe move 0 0 ${dstep} units box
displace_atoms probe move 0 0 0 units box
min_style cg
min_modify dmax 0.3*${lj_to_metal} line quadratic
min_modify dmax 0.3*2.5701999999999998181 line quadratic
minimize 0.0 1.0e-6 100000 100000
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Memory usage per processor = 3.98023 Mbytes
Step Temp E_pair E_mol TotEng Press 
      62 6.1911061e-25   -188.74495            0   -188.61625   -3933.5283 
      63 6.1911061e-25   -188.74495            0   -188.61625   -3933.5286 
Loop time of 0.039083 on 1 procs for 1 steps with 150 atoms

Minimization stats:
  Stopping criterion = force tolerance
  Energy initial, next-to-last, final = 
        -188.616253538     -188.616253538     -188.616253538
  Force two-norm initial, final = 4.25777e-07 4.69203e-07
  Force max component initial, final = 7.60257e-08 9.98731e-08
  Final line search alpha, max atom move = 1 9.98731e-08
  Iterations, force evaluations = 1 2

Pair  time (%) = 0.00098896 (2.53041)
Neigh time (%) = 0 (0)
Comm  time (%) = 2.86102e-05 (0.0732038)
Outpt time (%) = 0 (0)
Other time (%) = 0.0380654 (97.3964)
Gffft  time (%) = 0.000912905 (2.33581)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    789 ave 789 max 789 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11689 ave 11689 max 11689 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11689
Ave neighs/atom = 77.9267
Neighbor list builds = 0
Dangerous builds = 0
min_style hftn
minimize 0.0 1.0e-6 100000 100000
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Memory usage per processor = 4.90949 Mbytes
Step Temp E_pair E_mol TotEng Press 
      63 6.1911061e-25   -188.74495            0   -188.61625   -3933.5286 
Loop time of 0.000566959 on 1 procs for 0 steps with 150 atoms

Minimization stats:
  Stopping criterion = force tolerance
  Energy initial, next-to-last, final = 
        -188.616253538     -188.616253538     -188.616253538
  Force two-norm initial, final = 4.69203e-07 4.69203e-07
  Force max component initial, final = 9.98731e-08 9.98731e-08
  Final line search alpha, max atom move = 0 0
  Iterations, force evaluations = 0 0

Pair  time (%) = 0.000334024 (58.9151)
Neigh time (%) = 0 (0)
Comm  time (%) = 5.00679e-06 (0.883095)
Outpt time (%) = 0 (0)
Other time (%) = 0.000227928 (40.2019)
Gffft  time (%) = 0.000442982 (78.1329)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    789 ave 789 max 789 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11689 ave 11689 max 11689 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11689
Ave neighs/atom = 77.9267
Neighbor list builds = 0
Dangerous builds = 0
min_style fire
minimize 0.0 1.0e-6 100000 100000
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Memory usage per processor = 3.28329 Mbytes
Step Temp E_pair E_mol TotEng Press 
      63            0   -188.74495            0   -188.61625   -3933.5286 
      64 7.5184192e-25   -188.74495            0   -188.61625   -3933.5286 
Loop time of 0.038739 on 1 procs for 1 steps with 150 atoms

Minimization stats:
  Stopping criterion = force tolerance
  Energy initial, next-to-last, final = 
        -188.616253538     -188.616253538     -188.616253538
  Force two-norm initial, final = 4.69203e-07 4.69203e-07
  Force max component initial, final = 9.98731e-08 9.98731e-08
  Final line search alpha, max atom move = 0 0
  Iterations, force evaluations = 1 1

Pair  time (%) = 0.00067091 (1.73187)
Neigh time (%) = 0 (0)
Comm  time (%) = 9.05991e-06 (0.0233871)
Outpt time (%) = 0 (0)
Other time (%) = 0.038059 (98.2447)
Gffft  time (%) = 0.000603914 (1.55893)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    789 ave 789 max 789 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11689 ave 11689 max 11689 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11689
Ave neighs/atom = 77.9267
Neighbor list builds = 0
Dangerous builds = 0


#dump            dumpid type4 custom 1 dump.melt.*.gf x y z fx fy fz
#run 0
#undump          dumpid

next n
jump SELF nloop

undump 3
dump 3 all image 1 aimag.gf.*.${n}.jpg fz type adiam 0.8 size 1024 1024 view 90.0 0.0 center s 0.5 0.5 0.5 zoom 2 shiny 1 axes yes 0.8 0.02
dump 3 all image 1 aimag.gf.*.2.jpg fz type adiam 0.8 size 1024 1024 view 90.0 0.0 center s 0.5 0.5 0.5 zoom 2 shiny 1 axes yes 0.8 0.02
dump_modify 3 backcolor white amap -1.0 1.0 ca 0.0 4 min red 0.0 black 0.000001 navy max blue
dump_modify 3 pad 4

if "$n < 3"  then "variable dstep equal 0.0" else "variable dstep equal -1.6"
variable dstep equal 0.0


displace_atoms probe move 0 0 ${dstep} units box
displace_atoms probe move 0 0 0 units box
min_style cg
min_modify dmax 0.3*${lj_to_metal} line quadratic
min_modify dmax 0.3*2.5701999999999998181 line quadratic
minimize 0.0 1.0e-6 100000 100000
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Memory usage per processor = 3.98023 Mbytes
Step Temp E_pair E_mol TotEng Press 
      64 7.5184192e-25   -188.74495            0   -188.61625   -3933.5286 
      65 7.5184192e-25   -188.74495            0   -188.61625   -3933.5285 
Loop time of 0.0389712 on 1 procs for 1 steps with 150 atoms

Minimization stats:
  Stopping criterion = force tolerance
  Energy initial, next-to-last, final = 
        -188.616253538     -188.616253538     -188.616253538
  Force two-norm initial, final = 4.69203e-07 3.706e-07
  Force max component initial, final = 9.98731e-08 8.48995e-08
  Final line search alpha, max atom move = 1 8.48995e-08
  Iterations, force evaluations = 1 2

Pair  time (%) = 0.000973225 (2.49729)
Neigh time (%) = 0 (0)
Comm  time (%) = 1.90735e-05 (0.0489425)
Outpt time (%) = 0 (0)
Other time (%) = 0.0379789 (97.4538)
Gffft  time (%) = 0.000868082 (2.2275)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    789 ave 789 max 789 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11689 ave 11689 max 11689 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11689
Ave neighs/atom = 77.9267
Neighbor list builds = 0
Dangerous builds = 0
min_style hftn
minimize 0.0 1.0e-6 100000 100000
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Memory usage per processor = 4.90949 Mbytes
Step Temp E_pair E_mol TotEng Press 
      65 7.5184192e-25   -188.74495            0   -188.61625   -3933.5285 
Loop time of 0.00054884 on 1 procs for 0 steps with 150 atoms

Minimization stats:
  Stopping criterion = force tolerance
  Energy initial, next-to-last, final = 
        -188.616253538     -188.616253538     -188.616253538
  Force two-norm initial, final = 3.706e-07 3.706e-07
  Force max component initial, final = 8.48995e-08 8.48995e-08
  Final line search alpha, max atom move = 0 0
  Iterations, force evaluations = 0 0

Pair  time (%) = 0.000333786 (60.8167)
Neigh time (%) = 0 (0)
Comm  time (%) = 5.00679e-06 (0.91225)
Outpt time (%) = 0 (0)
Other time (%) = 0.000210047 (38.2711)
Gffft  time (%) = 0.000433683 (79.0182)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    789 ave 789 max 789 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11689 ave 11689 max 11689 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11689
Ave neighs/atom = 77.9267
Neighbor list builds = 0
Dangerous builds = 0
min_style fire
minimize 0.0 1.0e-6 100000 100000
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Memory usage per processor = 3.28329 Mbytes
Step Temp E_pair E_mol TotEng Press 
      65            0   -188.74495            0   -188.61625   -3933.5285 
      66 4.6904454e-25   -188.74495            0   -188.61625   -3933.5285 
Loop time of 0.038888 on 1 procs for 1 steps with 150 atoms

Minimization stats:
  Stopping criterion = force tolerance
  Energy initial, next-to-last, final = 
        -188.616253538     -188.616253538     -188.616253538
  Force two-norm initial, final = 3.706e-07 3.706e-07
  Force max component initial, final = 8.48995e-08 8.48995e-08
  Final line search alpha, max atom move = 0 0
  Iterations, force evaluations = 1 1

Pair  time (%) = 0.000684023 (1.75896)
Neigh time (%) = 0 (0)
Comm  time (%) = 8.10623e-06 (0.0208451)
Outpt time (%) = 0 (0)
Other time (%) = 0.0381958 (98.2202)
Gffft  time (%) = 0.000598907 (1.54008)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    789 ave 789 max 789 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11689 ave 11689 max 11689 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11689
Ave neighs/atom = 77.9267
Neighbor list builds = 0
Dangerous builds = 0


#dump            dumpid type4 custom 1 dump.melt.*.gf x y z fx fy fz
#run 0
#undump          dumpid

next n
jump SELF nloop

undump 3
dump 3 all image 1 aimag.gf.*.${n}.jpg fz type adiam 0.8 size 1024 1024 view 90.0 0.0 center s 0.5 0.5 0.5 zoom 2 shiny 1 axes yes 0.8 0.02
dump 3 all image 1 aimag.gf.*.3.jpg fz type adiam 0.8 size 1024 1024 view 90.0 0.0 center s 0.5 0.5 0.5 zoom 2 shiny 1 axes yes 0.8 0.02
dump_modify 3 backcolor white amap -1.0 1.0 ca 0.0 4 min red 0.0 black 0.000001 navy max blue
dump_modify 3 pad 4

if "$n < 3"  then "variable dstep equal 0.0" else "variable dstep equal -1.6"
variable dstep equal -1.6


displace_atoms probe move 0 0 ${dstep} units box
displace_atoms probe move 0 0 -1.6000000000000000888 units box
min_style cg
min_modify dmax 0.3*${lj_to_metal} line quadratic
min_modify dmax 0.3*2.5701999999999998181 line quadratic
minimize 0.0 1.0e-6 100000 100000
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Memory usage per processor = 3.98023 Mbytes
Step Temp E_pair E_mol TotEng Press 
      66 4.6904454e-25   -107.24965            0   -105.93653    46732.289 
     166 4.6904454e-25   -188.02541            0   -183.23646   -2403.2077 
Loop time of 3.96374 on 1 procs for 100 steps with 150 atoms

Minimization stats:
  Stopping criterion = force tolerance
  Energy initial, next-to-last, final = 
        -105.936527387     -183.236457841     -183.236457841
  Force two-norm initial, final = 96.8155 8.29539e-07
  Force max component initial, final = 37.5675 1.74185e-07
  Final line search alpha, max atom move = 1 1.74185e-07
  Iterations, force evaluations = 100 196

Pair  time (%) = 0.0923729 (2.33045)
Neigh time (%) = 0 (0)
Comm  time (%) = 0.00153708 (0.0387786)
Outpt time (%) = 3.78115 (95.3934)
Other time (%) = 0.0886858 (2.23742)
Gffft  time (%) = 0.0458152 (1.15586)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    789 ave 789 max 789 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11784 ave 11784 max 11784 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11784
Ave neighs/atom = 78.56
Neighbor list builds = 0
Dangerous builds = 0
min_style hftn
minimize 0.0 1.0e-6 100000 100000
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Memory usage per processor = 4.90949 Mbytes
Step Temp E_pair E_mol TotEng Press 
     166 4.6904454e-25   -188.02541            0   -183.23646   -2403.2077 
Loop time of 0.000586987 on 1 procs for 0 steps with 150 atoms

Minimization stats:
  Stopping criterion = force tolerance
  Energy initial, next-to-last, final = 
        -183.236457841     -183.236457841     -183.236457841
  Force two-norm initial, final = 8.29539e-07 8.29539e-07
  Force max component initial, final = 1.74185e-07 1.74185e-07
  Final line search alpha, max atom move = 0 0
  Iterations, force evaluations = 0 0

Pair  time (%) = 0.000361919 (61.6572)
Neigh time (%) = 0 (0)
Comm  time (%) = 5.00679e-06 (0.852965)
Outpt time (%) = 0 (0)
Other time (%) = 0.00022006 (37.4898)
Gffft  time (%) = 0.000442266 (75.3452)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    789 ave 789 max 789 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11245 ave 11245 max 11245 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11245
Ave neighs/atom = 74.9667
Neighbor list builds = 0
Dangerous builds = 0
min_style fire
minimize 0.0 1.0e-6 100000 100000
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Memory usage per processor = 3.28329 Mbytes
Step Temp E_pair E_mol TotEng Press 
     166            0   -188.02541            0   -183.23646   -2403.2077 
     167 2.3500519e-24   -188.02541            0   -183.23646   -2403.2077 
Loop time of 0.0388811 on 1 procs for 1 steps with 150 atoms

Minimization stats:
  Stopping criterion = force tolerance
  Energy initial, next-to-last, final = 
        -183.236457841     -183.236457841     -183.236457841
  Force two-norm initial, final = 8.29539e-07 8.29539e-07
  Force max component initial, final = 1.74185e-07 1.74185e-07
  Final line search alpha, max atom move = 0 0
  Iterations, force evaluations = 1 1

Pair  time (%) = 0.000731945 (1.88252)
Neigh time (%) = 0 (0)
Comm  time (%) = 7.86781e-06 (0.0202356)
Outpt time (%) = 0 (0)
Other time (%) = 0.0381413 (98.0972)
Gffft  time (%) = 0.000586271 (1.50786)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    789 ave 789 max 789 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11245 ave 11245 max 11245 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11245
Ave neighs/atom = 74.9667
Neighbor list builds = 0
Dangerous builds = 0


#dump            dumpid type4 custom 1 dump.melt.*.gf x y z fx fy fz
#run 0
#undump          dumpid

next n
jump SELF nloop

#pair_coeff  * 3 lj/expand 0.0 ${lj_to_metal} 0.0 0.1
#pair_coeff  3 * lj/expand 0.0 ${lj_to_metal} 0.0 0.1
unfix 2
unfix fixID1
dump           surfdump  type5 custom 1 dumpsurfgf x y z fx fy fz
dump           probedump probe custom 1 dumpprobegf x y z fx fy fz
run 0
Center-of-mass of the elastic layer is 6.491250, 6.491250, 3.060000.
Absolute thickness of the ghost shell is 0.624240.
Memory usage per processor = 3.59565 Mbytes
Step Temp E_pair E_mol TotEng Press 
     167 2.3500519e-24   -188.02541            0   -183.23646   -2403.2077 
Loop time of 1.90735e-06 on 1 procs for 0 steps with 150 atoms

Pair  time (%) = 0 (0)
Neigh time (%) = 0 (0)
Comm  time (%) = 0 (0)
Outpt time (%) = 0 (0)
Other time (%) = 1.90735e-06 (100)
Gffft  time (%) = 0.000230074 (12062.5)

Nlocal:    150 ave 150 max 150 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    789 ave 789 max 789 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    11245 ave 11245 max 11245 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 11245
Ave neighs/atom = 74.9667
Neighbor list builds = 0
Dangerous builds = 0

shell "rm log.lammps"
GFMD pre_force time = 0.006982
GFMD post_force time = 0.067555
GFMD total time = 0.074537

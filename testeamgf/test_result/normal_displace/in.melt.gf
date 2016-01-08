# Perform a small simulation using EAM GF

log log.gf

units		metal
boundary        p p f
processors      1 1 1
atom_style	gfmd
newton          off

read_data       2255negbare_gfdatafile
#read_data       ../rb_fcc_h05brms.25restor_l1024x8_lmin11.8_eb_sq_flat_5typed

variable lj_to_metal equal 2.5702 # = (Au cell)/(LJ nn cell)= (4.08 nm)/(2^1/6 * 2^1/2 sigma)
variable rcut_interface equal 1.1224620483*${lj_to_metal}

variable rcut_expand equal ${rcut_interface}+4.0

pair_style hybrid eam/alloy lj/expand ${rcut_interface}
pair_coeff   * * eam/alloy Au-Grochola-JCP05.eam.alloy Au Au NULL Au Au

#pair_coeff   * * lj/cut    0.0 ${lj_to_metal} ${rcut_interface}
#pair_coeff   * * eam/alloy Au-Grochola-JCP05.eam.alloy NULL NULL NULL Au NULL

pair_coeff   1 3 lj/expand    0.0 ${lj_to_metal} 0.0 ${rcut_interface}
pair_coeff   2 3 lj/expand    0.0 ${lj_to_metal} 0.0 ${rcut_interface}
pair_coeff   3 3 lj/expand    0.0 ${lj_to_metal} 0.0 ${rcut_interface}
pair_coeff   3 4 lj/expand    0.0 ${lj_to_metal} 4.0 ${rcut_interface} #expand to decrease corrugation # disabled
pair_coeff   3 5 lj/expand    0.0 ${lj_to_metal} 0.0 ${rcut_interface}
pair_modify	shift yes


neighbor        4.0 bin
neigh_modify    every 1 delay 0 check yes
group           rbgroup type 3
group           gfmdgroup type 2 5
group           type2 type 2
group           type4 type 4
group           type5 type 5
#variable        maxgfht equal bound(type5,zmax)
#region          topgfslice block INF INF INF INF ${maxgfht}-0.1 ${maxgfht}+0.1 units box
#group           surfgfatoms region topgfslice
region          columnregion1 block 0.5 4.6 0.5 4.6 INF INF units box
region          columnregion2 block 2.6 6.7 4.6 6.7 INF INF units box
group           col1group region columnregion1
group           col2group region columnregion2
group           columns union col1group col2group
group           probe intersect columns type5 # surfgfatoms


#fix             bottomfix gfmdgroup setforce 0.0 0.0 0.0

thermo		500



#fix myfixdel type2 addforce 0.0 0.0 -0.01 # del this again


# Idea to have at rho of energy minimum.  
# Stiffness matrix (par term) is dFdrho * d2fdr2
#fix             mygfmdfix gfmdgroup gfmd gfmd fcc100fteam &
#                0.0 8.0 4 &
#                2 3 &
#                0 0 &
#                0 0 &
#                0 0 & 
#                0 0 & 
#                0 0 & 
#                height 4 # dump_greens_function
#fix_modify      mygfmdfix energy yes

#fix             mygfmdfix gfmdgroup gfmd gfmd fcc100fteam &
#                -0.075 8.3 4 &
#                0.001 0.001 &
#                0 0 &
#                0 0 &
#                0004 0 & 
#                000.2 0 & 
#                000.16 0 & 
#                height 4 # dump_greens_function
#
#
#fix_modify      mygfmdfix energy yes

# Caution about rho going negative?
#fix             mygfmdfix gfmdgroup gfmd gfmd fcc100fteam &
#                 0.0 8.0 4   &
#                 -0.1 .2    &
#                 -0.01 .02  &
#                 -.002 .004  &
#                 1000.4   0.0   &
#                 0.0   0.1   &
#                 0.0   0.05  &
#                 height 4 dump_greens_function
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

fix             mygfmdfix gfmdgroup gfmd gfmd fcc100fteam &
                -7.521826747793241E-002  8.32896611616331 4   &
                -5.368227746962540E-002  6.422825369266792E-002   &
                -3.041321389379986E-002 -2.049745662446555E-002   &
                -9.126569474822686E-003  3.545931446285652E-002   &
                 1.46855     -0.091631   &
                -.069519      0.053247   &
                -.038890      0.017848   &
                height 5  #dump_greens_function
fix_modify      mygfmdfix energy yes


#fix             mygfmdfix gfmdgroup gfmd gfmd fcc100fteam &
#                0 0 4   &
#                0 0 &
#                0 0 &
#                0 0 &
#                 1.5     -0.0916   &
#                -.07      0.0532   &
#                -.04      0.0178   &
#                height 5  #dump_greens_function
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
minimize 0.0 1.0e-6 100000 100000
min_style hftn
minimize 0.0 1.0e-6 100000 100000
min_style fire
minimize 0.0 1.0e-6 100000 100000

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
dump_modify 3 backcolor white amap -1.0 1.0 ca 0.0 4 min red 0.0 black 0.000001 navy max blue
dump_modify 3 pad 4

if "$n < 3"  then "variable dstep equal 0.0" else "variable dstep equal -1.6"


displace_atoms probe move 0 0 ${dstep} units box

min_style cg
min_modify dmax 0.3*${lj_to_metal} line quadratic
minimize 0.0 1.0e-6 100000 100000
min_style hftn
minimize 0.0 1.0e-6 100000 100000
min_style fire
minimize 0.0 1.0e-6 100000 100000


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

shell "rm log.lammps"

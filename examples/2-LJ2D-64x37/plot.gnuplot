#!/bin/sh
# To visualize the elastic stiffness coefficient obtained by FixGFC
# for (1+1)D LJ system; 20 atoms per "layer"
# 
# Retrive data for primitive cell
tail -64 LJ2D_prim.log > ESC.LJ2D1.dat
head -1 ESC.LJ2D1.dat >> ESC.LJ2D1.dat
#
# Generate gnuplot script
echo 'set term post enha colo 20' >  gpt.scr
echo 'set out "Phi_q_LJ2D_1.eps"' >> gpt.scr
echo 'set xr [0:2*pi]'            >> gpt.scr
echo 'set yr [-50:300]'           >> gpt.scr
echo 'set xlabel "q"'             >> gpt.scr
echo 'set ylabel "{/Symbol F} ({/Symbol e/s}^2)' >> gpt.scr
echo 'set title "Elastic Stiffness Coefficients for 2D LJ system, primitive cell"' >> gpt.scr
echo 'set xtics ("0" 0, "{/Symbol p}/2" pi/2,"{/Symbol p}" pi,"3{/Symbol p}/2" 1.5*pi,"2{/Symbol p}" 2*pi)' >> gpt.scr
#
echo -e 'plot "ESC.LJ2D1.dat" u ($0/32*pi):5 w lp pt 6 t "{/Symbol F}_{xx}",\c' >> gpt.scr
echo -e '"" u ($0/32*pi):8 w lp pt 5 t "{/Symbol F}_{xy}/i",\c' >> gpt.scr
echo '"" u ($0/32*pi):11 w lp pt 7 t "{/Symbol F}_{yy}"'    >> gpt.scr
#
# plot the figure and remove intermediate files
gnuplot gpt.scr
rm -rf ESC.LJ2D1.dat gpt.scr
#
# Retrive data for complex cell (2 atoms per unit cell)
tail -32 LJ2D_comp.log > ESC.LJ2D2.dat
head -1 ESC.LJ2D2.dat >> ESC.LJ2D2.dat
#
# Generate gnuplot script
echo 'set term post enha colo 20' >  gpt.scr
echo 'set out "Phi_q_LJ2D_2.eps"' >> gpt.scr
echo 'set xr [0:2*pi]'            >> gpt.scr
echo 'set yr [0:200]'             >> gpt.scr
echo 'set xlabel "q"'             >> gpt.scr
echo 'set ylabel "{/Symbol F} ({/Symbol e/s}^2)' >> gpt.scr
echo 'set title "Elastic Stiffness Coefficients for 2D LJ system, complex cell"' >> gpt.scr
echo 'set xtics ("0" 0, "{/Symbol p}/2" pi/2,"{/Symbol p}" pi,"3{/Symbol p}/2" 1.5*pi,"2{/Symbol p}" 2*pi)' >> gpt.scr
#
echo -e 'plot "ESC.LJ2D2.dat" u ($0/16*pi):5 w lp pt 5 t "{/Symbol F}_{11}",\c' >> gpt.scr
echo -e '"" u ($0/16*pi):15 w lp pt 6 t "{/Symbol F}_{22}",\c' >> gpt.scr
echo -e '"" u ($0/16*pi):25 w lp pt 7 t "{/Symbol F}_{33}",\c' >> gpt.scr
echo '"" u ($0/16*pi):35 w lp pt 9 t "{/Symbol F}_{44}"'    >> gpt.scr
#
# plot the figure and remove intermediate files
gnuplot gpt.scr
rm -rf ESC.LJ2D2.dat gpt.scr

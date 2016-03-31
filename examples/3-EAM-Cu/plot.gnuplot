set encoding iso
set xlabel 'q'
set ylabel '{/Symbol F} (eV/\305^2)'
set key 2.2,9.5 spacing 1.2
#
set title "Cu(111) Orthogonal Lattice with EAM Potential"
set xr [0:24]
set yr [0:10]
#
set arrow 1 from 4,0  to 4,10  lt -1 nohead
set arrow 2 from 12,0 to 12,10 lt -1 nohead
set arrow 3 from 16,0 to 16,10 lt -1 nohead
set arrow 4 from 24,0 to 24,10 lt -1 nohead
#
set xtics ("{/Symbol G}" 0, "X" 4, "M" 12, "{/Symbol G}" 16, "Y" 24)
#
plot 'disp.eam3d'  u 0:45 w l s u lt 1 t '',\
     ''            u 0:59 w l s u lt 1 t '',\
     ''            u 0:73 w l s u lt 1 t '',\
     ''            u 0:3  w p pt 4  ps 2 lt 2 t '{/Symbol F}_{11}',\
     ''            u 0:17 w p pt 6  ps 2 lt 2 t '{/Symbol F}_{22}',\
     ''            u 0:31 w p pt 12 ps 2 lt 2 t '{/Symbol F}_{33}',\
     ''            u 0:45 w p pt 5  ps 1 lt 3 t '{/Symbol F}_{44}',\
     ''            u 0:59 w p pt 7  ps 1 lt 3 t '{/Symbol F}_{55}',\
     ''            u 0:73 w p pt 13 ps 1 lt 3 t '{/Symbol F}_{66}'
#
set terminal postscript landscape solid color enhanced lw 3 "Helvetica" 20
set out "Phi_q_EAM3D.eps"
replot
unset out;
set terminal X11
replot

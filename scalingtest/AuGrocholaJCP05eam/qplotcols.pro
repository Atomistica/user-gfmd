; IDL script to plot 
; y = the height difference between layers of the lattice 
;  vs 
; x = height
;
; Created by Tristan 9:47-9:48 6/26/13

;a=read_ascii("dumpinitprobefull",data_start=10)
a=read_ascii("dumpinitallgf",data_start=10)
b=a.field1
z = b[2,*]
z = z(sort(z))
plot, z, z - shift(z,1), psym=4, yrange=[2.03-1.5, 2.03+1.5], /ystyle, thick=5
;plot, z, z - shift(z,1), psym=4, yrange=[-0.1, 40.10], /ystyle, thick=5
;plot, z, z - shift(z,1), psym=4, yrange=[min(z - shift(z,1)), max(z-shift(z,1))], /ystyle, thick=5
oplot, [-10,50], [2.04,2.04]
oplot, [-10,50], [4.07011,4.07011]/2.


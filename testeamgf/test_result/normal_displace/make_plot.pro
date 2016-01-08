

   
a=read_ascii("error_out.txt")
b=a.field1
disp = [0,0.1,0.4,1.6]
tgt= b[0,*]
err= b[1,*]
plot, disp, tgt, psym= 4, background='ffffff'x, color=0, thick=2, yrange=[0,5], xtitle="Displacement magnitude (nm)", ytitle="Force magnitude on displaced test atoms"
xyouts, 1, 0.5, "Length of global force error vector", color='ff'x
oplot, disp, err, psym = 5, color='ff'x, thick=2
xyouts, 0.1, 4, "GF successful at small displacements.  Error (red) is small frac. of full-atom sim (black) and prediction.", color=0
xyouts, 0.1,3.6,"   Closer investigation reveals error is linear rather than quadratic; finding a bug should improve this.", color=0
write_png, "test_result.png", tvrd(true=1)
exit


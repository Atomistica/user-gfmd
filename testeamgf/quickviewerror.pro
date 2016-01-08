

;pro quickviewerror


; Read File 1
file1dat = (read_ascii("dumpprobefull",  data_start=9)).field1
; Assign the data to arrays
x1 =file1dat[0,*] & y1 =file1dat[1,*] & z1 =file1dat[2,*]
fx1=file1dat[3,*] & fy1=file1dat[4,*] & fz1=file1dat[5,*]
; Sort arrays to unique order
order = sort(x1 * max(y1) * max(z1) + y1 * max(z1) + z1)
x1 =x1(order)  & y1 =y1(order)  & z1 =z1(order)
fx1=fx1(order) & fy1=fy1(order) & fz1=fz1(order)



; Read File 2
file2dat = (read_ascii("dumpprobegf",data_start=9)).field1
; Assign the data to arrays
x2 =file2dat[0,*] & y2 =file2dat[1,*] & z2 =file2dat[2,*]
fx2=file2dat[3,*] & fy2=file2dat[4,*] & fz2=file2dat[5,*]
; Sort arrays to unique order
order = sort(x2 * max(y2) * max(z2) + y2 *max(z2) + z2)
x2 =x2(order)  & y2 =y2(order)  & z2 =z2(order)
fx2=fx2(order) & fy2=fy2(order) & fz2=fz2(order)


print, '  ', "  x", "   y", "  z","  fx","  fy","  fz"
for i=0,6 do begin & $
  print, "full", x1[i], y1[i], z1[i], fx1[i], fy1[i], fz1[i] & $
  print, "gf  ", x2[i], y2[i], z2[i], fx2[i], fy2[i], fz2[i] & $
end


; Magnitude of vec{f1} - vec{f2}
force_sqerr_arr =(fx2 - fx1)^(2.) +  $
                 (fy2 - fy1)^(2.) +  $
                 (fz2 - fz1)^(2.) 

force_sq_arr =(fx1)^(2.) +  $
              (fy1)^(2.) +  $
              (fz1)^(2.) 

force = sqrt(total(force_sq_arr))
force_err = sqrt(total(force_sqerr_arr))
frac_err = force_err / force

print, "Golbal force vector: ", force
print, "Error in force on probe atoms: ", force_err
print, "Fractional error in force on probe atoms: ", frac_err




;plot,background='ffffff'x,color=0,psym=4,x1,y1
openu ,1,"error_out.txt",/append
printf,1,force, force_err, frac_err
close, 1
;end


exit



; Quick script to plot the force difference between GF and full dump files.
; A dump file exists for each "probe displacement magnitude"
; Tristan 9/18/2013

; GFMD is 
; Order u^2 accurate in energy
; Order u accurate in force

; So if GFMD is working, error in force should be quadratic with displacement u

; To match up probe atoms in the two probe dump files based on location, use this algorithm:
; Travervse atom1 array: for jj=0...N-1      * * * * * * *  and find closest atom2 (at location "here")
; put the "here" values into sortorder  =   [0,1,3,3,2,4,X]
; Look for duplicate pairing and record where no pairing (X) was found in recallbadindices


function diffmag, fx1,fy1,fz1,fx2,fy2,fz2
  return, sqrt((fx1-fx2)^(2d)+(fy1-fy2)^(2d)+(fz1-fz2)^(2d))
end


function find_str_at_end_of_file, str, file
  spawn, "tac "+file+" | grep -m 1 "+str, outputlines
  mystr = outputlines(where(strmatch(outputlines, "*"+str+"*"))) 
  return, strmid(mystr, strpos(mystr,str)+ strlen(str) ) 
end


pro force_analysis

; Find the dstep folders, containing the output logs
a=file_search('*dstep*/')
if n_elements(a) le 1 then print, "No dstep folders found. Please call script from output folder. See README"

; Hard-code the probe displacement magnitude in each folder, ask user to check
u=float(strmid(a,5,15))

; Print info for user to check that hard-coded array above matches log file contents and names
print, '1)These are what this IDL analysis script (force_analysis) has hard-coded as dstep'
print, '2)These are the dstep file folder names:'
print, '3)These are the displacement magnitudes as recorded in the gf and full-atom log files'
dsteparr = strarr(n_elements(u)) & gfdsteparr = dsteparr
for i = 0, n_elements(a)-1 do $
  dsteparr[i] = find_str_at_end_of_file('displace_atoms\ probe1\ move', a[i]+'/log.full')
for i = 0, n_elements(a)-1 do $
  gfdsteparr[i] = find_str_at_end_of_file('displace_atoms\ probe1\ move', a[i]+'/log.gf')
for i = 0, n_elements(a)-1 do $
  print, u[i], ' ', a[i], ' ', dsteparr[i], ' ', gfdsteparr[i]
print, "Please check! It is important that the values in each line above match." & print, "" & wait, 5

; Initialize arrays
err=0d * u
absmag = 0d * u
nn_scale = 1.0 ; approx distance scale between atoms for error checking

print, "Now we will print the contents of the probe dump files, with each"
print, "probe atom identified in the GF run with the corresponding one in the full run"

for i=0, n_elements(a)-1 do begin 
  print, '' & print, a[i] & cd, a[i]
  a1=read_ascii("dumpprobefull",data_start=9)
  a2=read_ascii("dumpprobegf"  ,data_start=9)
  file1dat=a1.field1
  file2dat=a2.field1

  ; Assign the data to arrays
  x1 =file1dat[0,*] & y1 =file1dat[1,*] & z1 =file1dat[2,*]
  fx1=file1dat[3,*] & fy1=file1dat[4,*] & fz1=file1dat[5,*]
  x2 =file2dat[0,*] & y2 =file2dat[1,*] & z2 =file2dat[2,*]
  fx2=file2dat[3,*] & fy2=file2dat[4,*] & fz2=file2dat[5,*]
  z1-=min(z1)
  z2-=min(z2)

  ; Match up each atom in first array with the nearby one in the second file
  sortorder = intarr(n_elements(x1))
  recallbadindices = intarr(n_elements(x1))
  no_err = 1
  for jj=0, n_elements(x1)-1 do begin
    distsq_of_atoms_to_atomjj = (x2-x1[jj])^(2d) + (y2-y1[jj])^(2d) + (z2-z1[jj])^(2d)
    min_distsq_of_atoms_to_atomjj = min(distsq_of_atoms_to_atomjj)

    if min_distsq_of_atoms_to_atomjj gt (nn_scale/4.)^(2.) then begin
      if no_err then print, "Warning: location-based atom ID-ing is breaking down.  Maybe due to PBCs."
      no_err = 0
      recallbadindices[jj] = 1 ; Atom1[jj] has nothing near it and should never be used
    end else begin
      here = (where(distsq_of_atoms_to_atomjj eq min(distsq_of_atoms_to_atomjj)))
      if n_elements(here) ne 1 then stop, "Error: Possibly multiple atoms in file2 equally near atom jj in file 1"
      sortorder[jj]=here
    end
  end

  x2 =x2(sortorder)  & y2 =y2(sortorder)  & z2 =z2(sortorder)
  fx2=fx2(sortorder) & fy2=fy2(sortorder) & fz2=fz2(sortorder)

  ; print for visual check
  print, "        ", "     x", "   y", "  z","  fx","  fy","  fz"
  for jj=0,n_elements(x1)-1 do begin
    if ~recallbadindices[jj] then begin ;print, "ignore the following line:"
      print, "full", x1[jj], y1[jj], z1[jj], fx1[jj], fy1[jj], fz1[jj] 
      print, "gf  ", x2[jj], y2[jj], z2[jj], fx2[jj], fy2[jj], fz2[jj] 
    end
    if ~(recallbadindices[jj]) && $
        (diffmag(x1[jj],y1[jj],z1[jj],x2[jj],y2[jj],z2[jj]) gt (nn_scale/4.)) $
           then stop, "Fatal error: The comparison is invalid but was not flagged!!"
  end

  ; Find the 2-norm of the force vector and error vector
  for jj=0,n_elements(x1)-1 do begin
    if ~(recallbadindices(jj)) then begin
      err[i]    += diffmag(fx1[jj],fy1[jj],fz1[jj],fx2[jj],fy2[jj],fz2[jj])
      absmag[i] += diffmag(fx1[jj],fy1[jj],fz1[jj],0,0,0)
    end
  end

  
cd, '..'
end

; plot the error in the gfmd data(black)
window, 0, xsize=500, ysize=500
plot, u, err, psym=4, thick=4, background='ffffff'x, color=0, /xlog, /ylog,$
title="Force Response (and Error) vs Displacement",charsize=1.3,$
yrange=[min(err)/10d,max(absmag)]

print, "The grey line is a linear guideline, red is quadratic."
print, "Restoring force should be approximately linear, then noisy at large displacements"
print, "If GFMD is working, there should be a large quadratic range of the error (large data points)"
print, "Or check agains subversioned .png files of output"
; overlay a quadratic guide on error and linear guide on magnitdue
modelx =u;findgen(1000)/1000. * max(u)
modely =modelx*modelx / (max(modelx * modelx)) * max(err)
oplot, color='ff'x, modelx, modely
modely = modelx /max(modelx) * max(absmag)
oplot, color='777777'x, modelx, modely

; what percentage error is that?
oplot, color='555555'x, u, absmag, psym=4, thick=2

; print also the convergence criterion; should be force tolerance
for i = 0, n_elements(a)-1 do $
  print, "Convergence criterion: ",  find_str_at_end_of_file('criterion', a[i]+'/log.full')

stop

end

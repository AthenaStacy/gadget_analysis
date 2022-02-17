pro turb_cosmo
close,2

filename0 = '/nobackupp12/astacy/turbbin_zoom10_ref4_'

filenum0 = '0537'
;filenum0 = '0538'
;filenum0 = '0300'

arraysize0 = 3213997
;arraysize0 = 2031372

;;;;;;;;;;;;;;;;;;;;;;;;
;read in redshifts
;;;;;;;;;;;;;;;;;;;;;;;
Tcmb = 2.7e0 * 22.0
h=0.7

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;read in zeroth data file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
arr1=dblarr(20,arraysize0)

close,1
openr, 1, filename0+filenum0
readf,1,arr1
close,1

cs0 = arr1(0,*)
id0 = arr1(1,*)
x0 = arr1(2,*)
y0 = arr1(3,*)
z0 = arr1(4,*)
temp0 = arr1(5,*)
nh0 = arr1(6,*) * 1.22 * .76
rad0 = arr1(11,*) 
vrad0 = arr1(12,*)
vrot0 = arr1(13,*)
mass0 = arr1(15,*);*1.e10/h
rho0 = nh0*1.67e-24
vx0 = arr1(16,*)
vy0 = arr1(17,*)
vz0 = arr1(18,*)
turb0 = arr1(19,*)
mach0 = turb0/cs0


;v_theta = theta_hat dot v
;theta_hat = r_hat cross vrot_hat
;vrot_hat = vrot / |vrot|

indexS = where(nh0 eq max(nh0))
xS = x0[indexS]
yS = y0[indexS]
zS = z0[indexS]

xS = xS[0]
yS = yS[0]
zS = zS[0]

vrotx = (y0 - yS)*vz0 - (z0 - zS)*vy0
vroty = (z0 - zS)*vx0 - (x0 - xS)*vz0
vrotz = (x0 - xS)*vy0 - (y0 - yS)*vx0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

A = FINDGEN(16) * (!PI*2/16.)
USERSYM, 0.8*COS(A),0.8*SIN(A),/FILL
set_plot,'ps'
device,/encapsulated
device,bits_per_pixel=8,/color,/portrait
device,yoffset=8.4
;device,xsize=13.0,scale_factor=1.0
;device,ysize=12.0,scale_factor=1.0
!P.charsize=1.5
!P.charthick=2.5
loadct,39,/silent
distinct_colors=10
col=0.6*bytscl((indgen(distinct_colors) mod 16)+1)+64
;!p.multi=[0,2,2]
!p.multi=0

!P.charsize=1.5
!p.multi = 0
;device,xoffset=-5.0
device,/encapsulated

nmin = 1.e-1
nmax = max(nh0)

rmin = 1.e1
rmax = max(rad0)

imin  = alog10(nmin)
imax = alog10(nmax)
imin2  = alog10(rmin)

print, 'imin = ', imin
print, 'imax = ', imax

len = fix(imax) - fix(imin) + 1
len = len*10

nbin = make_array(len, /double, value=0)
rbin = make_array(len, /double, value=0)
print, nbin 

for i = 0,len-1 do begin
  expo = float(imin) + float(i)/10.0
  ;print, 'expo = ', expo
  nbin[i] = 10.^expo
  expo2 = float(imin2) + float(i)/5.0
  rbin[i] = 10.^expo2
endfor
print, 'rbin = ', rbin

mach_bin = make_array(len, /double, value=1e-10)
turb_bin = make_array(len, /double, value=1e-10)
vrot_bin = make_array(len, /double, value=1e-10)
vrad_bin = make_array(len, /double, value=1e-10)
temp_bin = make_array(len, /double, value=1e-10)

for i = 0,len-2 do begin
   index = where(nh0 gt nbin[i] and nh0 lt nbin[i+1])
   ;index = where(rad0 lt rbin[i+1] and rad0 gt 0)
   mach_bin[i] = mean(mach0[index])
   turb_bin[i] = sqrt(mean(turb0[index]*turb0[index]))
   vrot_bin[i] = mean(vrot0[index])
   vrad_bin[i] = sqrt(mean(vrad0[index]*vrad0[index])) 
   temp_bin[i] = mean(temp0[index])

   ;index_alt = where(rad0 gt rbin[i] and rad0 lt rbin[i+1])
   ;nbin[i] = mean(nh0[index_alt]) 
endfor

turb_bin[len-1] = turb_bin[len-2]
temp_bin[len-1] = temp_bin[len-2]

print, 'mach = ', mach_bin
print, 'vrot = ', vrot_bin
print, 'vrad = ', vrad_bin


print, 'nbin_gadget = ['
for i = 0, len-1 do begin
    if nbin[i] gt 1.e-3 then begin 
       print, nbin[i], Format='(F16.6, $)' 
       print, ', ',  Format = '(A2, $)'
    endif
endfor

print, ''
print, 'vturb_gadget = ['
for i = 0, len-1 do begin
    if nbin[i] gt 1.e-3 then begin
       print, turb_bin[i],  Format='(F16.6, $)'
       print, ', ',  Format = '(A2, $)' 
    endif
endfor

print, ''
print, 'nbin = ', nbin
print, 'turb_bin = ', turb_bin

setup = 2

if setup eq 1 then begin

plot, nbin, abs(vrad_bin), xtitle='!6 n!DH!N [cm!U-3!N]', ytitle='!6velocity [km s!U-1!N]',title='',/xstyle,/ystyle, /xlog, xrange=[nmin,nmax],yrange=[0,5], thick=4.0, xmargin = [7,7]
axis, yaxis=1, yrange=[0,2],ystyle=1, ytitle='M!Dturb!N', ycharsize=1.0, color=50;
;plots, [nmin,nmax], [4.0*(1.0-1), 4.0*(1.0-1)], thick=3.0, color=50
;oplot, nbin, vrot_bin, linestyle=1, thick=5.0
oplot, nbin, turb_bin, linestyle=1, thick=6.0
oplot, nbin, 2.5*mach_bin, linestyle=5, thick=4.0, color=50

fac  = 1e5
fac2 = 0.
x0 = 1.5e0*fac
x1 = 1.5e1*fac
y2 = 4.7 - fac2
y1 = 4.2 - fac2
y0 = 3.7 - fac2
csize = 1.1
plots, [x0,x1], [y2, y2],  linestyle=0, thick=4.0
plots, [x0,x1], [y1, y1],  linestyle=1, thick=6.0
plots, [x0,x1], [y0, y0],  linestyle=5, thick=4.0, color=50

x3=3e1*fac

abs_val_string = string(33B)
xyouts, x3, y2, '!6v!Drad!N!6', charsize=csize
;xyouts, x3, y1, 'v!Drot!N', charsize=csize
xyouts, x3, y1, 'v!Dturb!N', charsize=csize
xyouts, x3, y0, 'M!Dturb!N', charsize=csize

xyouts, 1e-4*x3, y2, 'Gadget', charsize=1.3

endif

if setup eq 2 then begin

plot, nbin, temp_bin, xtitle='!6 n!DH!N [cm!U-3!N]', ytitle='!6Temp [K]',title='',/xstyle,/ystyle, /xlog, xrange=[nmin,nmax],yrange=[0,1200], thick=3.0

endif


USERSYM, 0.3*COS(A),0.3*SIN(A)

!P.charsize=1.5
!P.charsize=1.
!P.charthick=1.0
device,/close_file
set_plot,'x'
end

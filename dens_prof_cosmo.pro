pro dens_prof_cosmo

arraysize=500
ncol=11

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;read file
arraysize2=100
arr1=dblarr(ncol,arraysize)

close,1
openr,1,'bin_zoom10_ref4_newtstep_dens_0437.dat'
readf,1,arr1
close,1
rad1a = arr1(0,*)   ;radius in AU
menc1a = arr1(1,*)
menc_nosink1a = arr1(2,*)
nh1a = arr1(3,*)
h21a = arr1(4,*)
c_s1a = arr1(5,*)
temp1a = arr1(6,*)
mbe1a = arr1(7,*)
;vrot_loc1a = arr1(8,*)
curl1a = arr1(8,*)
vrad1a = arr1(9,*)
vrot1a = arr1(10,*)

;vort1a = vrot1a / (rad1a * 1.5e8) 
;vort1a = vrot_loc1a / (rad1a * 1.5e8)  
;vort1a = vort1a*vort1a  ;vorticity^2 in s^-1

vort1a = curl1a*curl1a

mrat1a = menc1a/mbe1a

entr1a = (temp1a/300.d)*((1.d4/nh1a)^0.1d)
vkep1a = (6.67d-8*menc1a*2d33/(rad1a*1.5d13))^0.5
vkep1a = vkep1a/1d5
ind1a = where(nh1a gt 1d-2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;read file
arraysize2=100
arr1=dblarr(ncol,arraysize)

close,1
openr,1,'bin_zoom10_ref3_newtstep_dens_0260.dat'
readf,1,arr1
close,1
rad1b = arr1(0,*)   ;radius
menc1b = arr1(1,*)
menc_nosink1b = arr1(2,*)
nh1b = arr1(3,*)
h21b = arr1(4,*)
c_s1b = arr1(5,*)
temp1b = arr1(6,*)
mbe1b = arr1(7,*)
;vrot_loc1b = arr1(8,*)
curl1b = arr1(8,*)
vrad1b = arr1(9,*)
vrot1b = arr1(10,*)

;vort1b = vrot1b / (rad1b * 1.5e8)
;vort1b = vrot_loc1b / (rad1b * 1.5e8)  
;vort1b = vort1b*vort1b  ;vorticity^2 in s^-1

vort1b = curl1b*curl1b

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;read file
arraysize2=100
arr1=dblarr(ncol,arraysize)

close,1
openr,1,'bin_zoom1_ref4_newtstep_dens_0439.dat'
readf,1,arr1
close,1
rad2a = arr1(0,*)   ;radius
menc2a = arr1(1,*)
menc_nosink2a = arr1(2,*)
nh2a = arr1(3,*)
h22a = arr1(4,*)
c_s2a = arr1(5,*)
temp2a = arr1(6,*)
mbe2a = arr1(7,*)
;vrot_loc1b = arr1(8,*)
curl2a = arr1(8,*)
vrad2a = arr1(9,*)
vrot2a = arr1(10,*)

vort2a = curl2a*curl2a

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;read file
arraysize2=100
arr1=dblarr(ncol,arraysize)

close,1
openr,1,'bin_zoom1_ref3_newtstep_dens_0230.dat'
readf,1,arr1
close,1
rad2b = arr1(0,*)   ;radius
menc2b = arr1(1,*)
menc_nosink2b = arr1(2,*)
nh2b = arr1(3,*)
h22b = arr1(4,*)
c_s2b = arr1(5,*)
temp2b = arr1(6,*)
mbe2b = arr1(7,*)
;vrot_loc1b = arr1(8,*)
curl2b = arr1(8,*)
vrad2b = arr1(9,*)
vrot2b = arr1(10,*)


vort2b = curl2b*curl2b

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;values from Turk et al

arr1=dblarr(2,39)
close,1
openr,1,'turk_h2.csv'
readf,1,arr1
close,1
rad_h2 = arr1(0,*) / 1.5e13
h2_turk = arr1(1,*) / 2

arr1=dblarr(2,16)
close,1
openr,1,'turk_h2b.csv'
readf,1,arr1
close,1
rad_h2b = arr1(0,*) / 1.5e13
h2_turkb = arr1(1,*) / 2


arr1=dblarr(2,26)
close,1
openr,1,'turk_menc.csv'
readf,1,arr1
close,1
rad_menc = arr1(0,*) / 1.5e13
menc_turk = arr1(1,*)

arr1=dblarr(2,24)
close,1
openr,1,'turk_mencb.csv'
readf,1,arr1
close,1
rad_mencb = arr1(0,*) / 1.5e13
menc_turkb = arr1(1,*)


arr1=dblarr(2,45)
close,1
openr,1,'turk_rho.csv'
readf,1,arr1
close,1
rad_rho = arr1(0,*) / 1.5e13
rho_turk = arr1(1,*)
rad_nh = rad_rho
nh_turk = rho_turk / 1.67e-24

arr1=dblarr(2,28)
close,1
openr,1,'turk_rhob.csv'
readf,1,arr1
close,1
rad_rhob = arr1(0,*) / 1.5e13
rho_turkb = arr1(1,*)
rad_nhb = rad_rhob
nh_turkb = rho_turkb / 1.67e-24


arr1=dblarr(2,61)
close,1
openr,1,'turk_temp.csv'
readf,1,arr1
close,1
rad_temp = arr1(0,*) / 1.5e13
temp_turk = arr1(1,*)

arr1=dblarr(2,41)
close,1
openr,1,'turk_tempb.csv'
readf,1,arr1
close,1
rad_tempb = arr1(0,*) / 1.5e13
temp_turkb = arr1(1,*)


arr1=dblarr(2,46)
close,1
openr,1,'turk_vort.csv'
readf,1,arr1
close,1
rad_vort = arr1(0,*) / 1.5e13
vort_turk = arr1(1,*)

arr1=dblarr(2,29)
close,1
openr,1,'turk_vortb.csv'
readf,1,arr1
close,1
rad_vortb = arr1(0,*) / 1.5e13
vort_turkb = arr1(1,*)


arr1=dblarr(2,102)
close,1
openr,1,'turk_vrad.csv'
readf,1,arr1
close,1
rad_vrad = arr1(0,*) / 1.5e13
vrad_turk = arr1(1,*)

arr1=dblarr(2,56)
close,1
openr,1,'turk_vradb.csv'
readf,1,arr1
close,1
rad_vradb = arr1(0,*) / 1.5e13
vrad_turkb = arr1(1,*)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

set_plot,'ps'
device,/encapsul
device,bits_per_pixel=8,/color,/portrait
device,yoffset=8.4
device,xsize=13.0,scale_factor=1.0
device,ysize=12.0,scale_factor=1.0

;!P.charsize=0.8
!P.charsize=1.3
!P.charthick=3.0
loadct,12,/silent
distinct_colors=10
col=0.6*bytscl((indgen(distinct_colors) mod 16)+1)+64
;!p.multi=[0,2,2]
!p.multi=0

xmin = 1.e1
xmax = 1.e7

mmin=1.e0
mmax=3.e3

smin = 1e1
smax = 1e5

qmin = 1e-1
qmax = 1e1

setup=1
;setup = 2
;setup=3
;setup=5
;setup=5.5
;setup = 7.5
;setup=8


if setup eq 1 then begin

!p.multi=[0,2,3]
!P.charsize=1.3

plot, [0,0], [0,0], xtitle =  '!6r [AU]', ytitle = 'M!Denc!N [M!I!9n!6!N]', /xlog, /ylog, thick=3.0, /xstyle,/ystyle, xrange=[xmin,xmax], yrange=[1.e-3,1.e5]
oplot, rad1a, menc1a, linestyle=0, thick=3.0
oplot, rad1b, menc1b, linestyle=2, thick=3.0
oplot, rad2a, menc2a, linestyle=0, thick=3.0, color=200
oplot, rad2b, menc2b, linestyle=2, thick=3.0, color=200
oplot, rad_menc, menc_turk, linestyle=0, thick=2.0, color=100
oplot, rad_mencb, menc_turkb, linestyle=2, thick=2.0, color=100

plot, [0,0], [0,0], xtitle =  '!6r [AU]', ytitle = 'v!Drad!N [km s!U-1!N]', /xlog, thick=3.0, /xstyle,/ystyle, xrange=[xmin,xmax], yrange=[-6,0]
oplot, rad1a, vrad1a, linestyle=0, thick=3.0
oplot, rad1b, vrad1b, linestyle=2, thick=3.0
oplot, rad2a, vrad2a, linestyle=0, thick=3.0, color=200
oplot, rad2b, vrad2b, linestyle=2, thick=3.0, color=200
oplot, rad_vrad, vrad_turk, linestyle=0, thick=2.0, color=100
oplot, rad_vradb, vrad_turkb, linestyle=2, thick=2.0, color=100

plot, [0,0], [0,0], xtitle =  '!6r [AU]', ytitle = 'n [cm!U-3!N]', /xlog, /ylog, thick=3.0, /xstyle,/ystyle, xrange=[xmin,xmax], yrange=[1.e2,1.e14]
oplot, rad1a, nh1a, linestyle=0, thick=3.0
oplot, rad1b, nh1b, linestyle=2, thick=3.0
oplot, rad2a, nh2a, linestyle=0, thick=3.0, color=200
oplot, rad2b, nh2b, linestyle=2, thick=3.0, color=200
oplot, rad_nh, nh_turk, linestyle=0, thick=2.0, color=100
oplot, rad_nhb, nh_turkb, linestyle=2, thick=2.0, color=100

plot, [0,0], [0,0], xtitle =  '!6r [AU]', ytitle = 'T [k]', /xlog, thick=3.0, /xstyle,/ystyle, xrange=[xmin,xmax], yrange=[0,2500]
oplot, rad1a, temp1a, linestyle=0, thick=3.0
oplot, rad1b, temp1b, linestyle=2, thick=3.0
oplot, rad2a, temp2a, linestyle=0, thick=3.0, color=200
oplot, rad2b, temp2b, linestyle=2, thick=3.0, color=200
oplot, rad_temp, temp_turk, linestyle=0, thick=2.0, color=100
oplot, rad_tempb, temp_turkb, linestyle=2, thick=2.0, color=100

plot, [0,0], [0,0], xtitle =  '!6r [AU]', ytitle = 'f!DH2!N', /xlog, /ylog, thick=3.0, /xstyle,/ystyle, xrange=[xmin,xmax], yrange=[1.e-4,1.e0]
oplot, rad1a, h21a, linestyle=0, thick=3.0
oplot, rad1b, h21b, linestyle=2, thick=3.0
oplot, rad2a, h22a, linestyle=0, thick=3.0, color=200
oplot, rad2b, h22b, linestyle=2, thick=3.0, color=200
oplot, rad_h2, h2_turk, linestyle=0, thick=2.0, color=100
oplot, rad_h2b, h2_turkb, linestyle=2, thick=2.0, color=100

plot, [0,0], [0,0], xtitle =  '!6r [AU]', ytitle = 'vorticity!U2!N [s!U-2!N]', /xlog, /ylog, thick=3.0, /xstyle,/ystyle, xrange=[xmin,xmax], yrange=[1.e-31,1.e-14]
oplot, rad1a, vort1a, linestyle=0, thick=3.0
oplot, rad1b, vort1b, linestyle=2, thick=3.0
oplot, rad2a, vort2a, linestyle=0, thick=3.0, color=200
oplot, rad2b, vort2b, linestyle=2, thick=3.0, color=200
oplot, rad_vort, vort_turk, linestyle=0, thick=2.0, color=100
oplot, rad_vortb, vort_turkb, linestyle=2, thick=2.0, color=100

print, 'h21a =', h21a

endif



if setup eq 2 then begin
!p.multi=[0,2,2]
!P.charsize=0.85

con_fac = 1.67d-24*3.14d7
con_fac = con_fac/2d33

tmin = 1.e-1
tmax = 1.e3

;vrad1a = -1 * abs(vrad1a)
;vrad1b = -1 * abs(vrad1b)
;vrad1c = -1 * abs(vrad1c)
;vrad1d = -1 * abs(vrad1d)

dims=size(vrad1a, /dimensions)
n1a=dims[0]
mshell1a = dblarr(n1a)
nshell1a = dblarr(n1a)
for i= 1L, n1a-1 do begin
  mshell1a[i] = menc1a[i] - menc1a[i-1]
endfor

dims=size(vrad1b, /dimensions)
n1b=dims[0]
mshell1b = dblarr(n1b)
for i= 1L, n1b-1 do begin
  mshell1b[i] = menc1b[i] - menc1b[i-1]
endfor

dims=size(vrad1c, /dimensions)
n1c=dims[0]
mshell1c = dblarr(n1c)
for i= 1L, n1c-1 do begin
  mshell1c[i] = menc1c[i] - menc1c[i-1]
endfor

dims=size(vrad1d, /dimensions)
n1d=dims[0]
mshell1d = dblarr(n1d)
for i= 1L, n1d-1 do begin
  mshell1d[i] = menc1d[i] - menc1d[i-1]
endfor

mdot_sphere1a =  abs(4.*3.14159*nh1a*((rad1a*1.5d13)^2)*vrad1a*1d5*con_fac)
mdot_sphere1b =  abs(4.*3.14159*nh1b*((rad1b*1.5d13)^2)*vrad1b*1d5*con_fac)
mdot_sphere1c =  abs(4.*3.14159*nh1c*((rad1c*1.5d13)^2)*vrad1c*1d5*con_fac)
mdot_sphere1d =  abs(4.*3.14159*nh1d*((rad1d*1.5d13)^2)*vrad1d*1d5*con_fac)

;mdot_sphere1a =  abs(mshell1a/((2*3.14158*rad1a*1.5d13))*vrad1a*1d5*3.14d7)
;mdot_sphere1b =  abs(mshell1b/((2*3.14158*rad1b*1.5d13))*vrad1b*1d5*3.14d7)
;mdot_sphere1c =  abs(mshell1c/((2*3.14158*rad1c*1.5d13))*vrad1c*1d5*3.14d7)
;mdot_sphere1d =  abs(mshell1d/((2*3.14158*rad1d*1.5d13))*vrad1d*1d5*3.14d7)

print, 'mdot =', mdot_sphere1a

ind1a = where(vrad1a lt 0)
ind1b = where(vrad1b lt 0)
ind1c = where(vrad1c lt 0)
ind1d = where(vrad1d lt 0)

tacc1a = menc_nosink1a/abs(mdot_sphere1a) 
tacc1b = menc_nosink1b/abs(mdot_sphere1b)
tacc1c = menc_nosink1c/abs(mdot_sphere1c)
tacc1d = menc_nosink1d/abs(mdot_sphere1d)

tfrag1a = mbe1a/abs(mdot_sphere1a)
tfrag1b = mbe1b/abs(mdot_sphere1b)
tfrag1c = mbe1c/abs(mdot_sphere1c)
tfrag1d = mbe1d/abs(mdot_sphere1d)

plot, [0,0], [0,0], xtitle =  '!6r [AU]', ytitle = 'dM/dt [M!I!9n!6!N yr!U-1!N]', /xlog, /ylog, thick=3.0, /xstyle,/ystyle, xrange=[xmin,xmax], yrange=[1.e-3,1.e1]
oplot, rad1a(ind1a), mdot_sphere1a(ind1a), linestyle=0, thick=4.0
oplot, rad1b(ind1b), mdot_sphere1b(ind1b), linestyle=1, thick=4.0, color=20
oplot, rad1c(ind1c), mdot_sphere1c(ind1c), linestyle=2, thick=2.0, color=100
oplot, rad1d(ind1d), mdot_sphere1d(ind1d), linestyle=3, thick=2.0, color=200

plot, [0,0], [0,0], xtitle =  '!6r [AU]', ytitle = 't!Dgrowth!N [yr]', /xlog, /ylog, thick=3.0, /xstyle,/ystyle, xrange=[xmin,xmax], yrange=[tmin,tmax]
oplot, rad1a(ind1a), tacc1a(ind1a), linestyle=0, thick=4.0
oplot, rad1b(ind1b), tacc1b(ind1b), linestyle=1, thick=4.0, color=20
oplot, rad1c(ind1c), tacc1c(ind1c), linestyle=2, thick=2.0, color=100
oplot, rad1d(ind1d), tacc1d(ind1d), linestyle=3, thick=2.0, color=200

plot, [0,0], [0,0], xtitle =  '!6r [AU]', ytitle = 't!Dfrag!N [yr]', /xlog, /ylog, thick=3.0, /xstyle,/ystyle, xrange=[xmin,xmax], yrange=[tmin,tmax]
oplot, rad1a(ind1a), tfrag1a(ind1a), linestyle=0, thick=4.0
oplot, rad1b(ind1b), tfrag1b(ind1b), linestyle=1, thick=4.0, color=20
oplot, rad1c(ind1c), tfrag1c(ind1c), linestyle=2, thick=2.0, color=100
oplot, rad1d(ind1d), tfrag1d(ind1d), linestyle=3, thick=2.0, color=200

plot, [0,0], [0,0], xtitle =  '!6r [AU]', ytitle = 't!Dfrag!N / t!Dgrowth!N', /xlog, /ylog, thick=3.0, /xstyle,/ystyle, xrange=[xmin,xmax], yrange=[1.e-1,1.e1]
oplot, rad1a(ind1a), tfrag1a(ind1a)/tacc1a(ind1a), linestyle=0, thick=4.0
oplot, rad1b(ind1b), tfrag1b(ind1b)/tacc1b(ind1b), linestyle=1, thick=4.0, color=20
oplot, rad1c(ind1c), tfrag1c(ind1c)/tacc1c(ind1c), linestyle=2, thick=2.0, color=100
oplot, rad1d(ind1d), tfrag1d(ind1d)/tacc1d(ind1d), linestyle=3, thick=2.0, color=200
endif


if setup eq 3 then begin

plot,  xtitle='!6r [AU]', rad1a, Q1a, xrange=[xmin, xmax],yrange=[qmin, qmax],ytitle='Q', /xlog, /ylog,  /xstyle, /ystyle, thick=2.0, linestyle=0;, xmargin = [10,0], ymargin = [4,4], ycharsize = 1.0;, ytickname = ['10', '10!U2!N', '10!U3!N', '10!U4!N'], ycharsize = 1.0
loadct, 12, /silent
oplot, rad1b, Q1b, linestyle=1, thick=4.0, color=20
oplot, rad1c, Q1c, linestyle=2, thick=2.0, color=100
oplot, rad1d, Q1d, linestyle=3, thick=2.0, color=200
plots, [xmin,xmax], [1, 1], linestyle = 0, color=100

endif


if setup eq 5 then begin

;ind1b = where(vrad1b lt 0)

vrad1a = abs(vrad1a)
vrad1b = abs(vrad1b)

plot, [0,0], [0,0], xtitle =  '!6r [AU]', ytitle = 'c!Ds!N [km s!U-1!N]', /xlog, thick=3.0, /xstyle,/ystyle, xrange=[xmin,xmax], yrange=[0,10]
oplot, rad1a, c_s1a, linestyle=0, thick=2.0
oplot, rad1b, c_s1b, linestyle=1, thick=4.0, color=20 
oplot, rad1c, c_s1c, linestyle=2, thick=2.0, color=100
oplot, rad1d, c_s1d, linestyle=3, thick=2.0, color=200

;oplot, rad1a(ind1a), vrad1a(ind1a), linestyle=1, thick=4.0, color=200
;oplot, rad1b(ind1b), vrad1b(ind1b), linestyle=1, thick=4.0

endif


!P.charsize=1.5
!P.charthick=1.0
device,/close_file
set_plot,'x'
end

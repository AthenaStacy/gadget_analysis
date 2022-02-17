pro g2plot_cosmo
close,2

arraysize = 286000  ;ref4
;arraysize = 1142000
arraysize2 = 89000 ;ref2
;arraysize = 71500  ;ref3
;arraysize = 1039   ;ref0

index_array=indgen(1,arraysize)

;filenum1 = '0030'
filenum1 = '0537'
;filenum1 = '0538'
;filenum1 = '0530'
filenum2 = '0172'

dir = '/nobackupp12/astacy/'

filename =  dir +'snapbin_zoom10_ref4_'
;filename =  dir +'snapbin2_zoom10_ref4_'
;filename =  dir +'snap_chem_bin_zoom10_ref4_'

filename2 =  dir + 'snapbin_zoom10_ref2_'

;filename2 = filename
;filenum2 = filenum1

;;;;;;;;;;;;;;;;;;;;;;;;
;read in redshifts
;;;;;;;;;;;;;;;;;;;;;;;
Tcmb = 2.7e0 * 22.0
h=0.7
G = 6.67e-8

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;read in main data file
;;;;;;;;;;;;;;;;;;;;;;;;;;;
arr1=dblarr(20,arraysize)

close,1
openr, 1, filename+filenum1
readf,1,arr1
close,1

zp1=27

error = arr1(0,*)
id =arr1(1,*)
x=arr1(2,*)
y=arr1(3,*)
z=arr1(4,*)
x1=(1.e0/h)*x/zp1
y1=(1.e0/h)*y/zp1
z1=(1.e0/h)*z/zp1
x1=1.000e3*x1
y1=1.000e3*y1
z1=1.000e3*z1
temp=arr1(5,*)
nh=arr1(6,*) * 1.22 * .76
HII=arr1(7,*)
H2I=arr1(8,*)
HDI=arr1(9,*);*4.e-5
gam= arr1(10,*)
rad=arr1(11,*)    ;dis in AU
vrad=arr1(12,*)
vrot=arr1(13,*)

hsm = arr1(14,*) ;in pc
mass=arr1(15,*);*1.e10/h
rho=nh * 1.67e-24*1.22

print, 'mass =', mass[10]

;bfac = 1.e9/1.06
bfac = 1.e9/1.8 
;bfac = 1.e0
bfieldx = arr1(16,*) * bfac
bfieldy = arr1(17,*) * bfac
bfieldz = arr1(18,*) * bfac
divb = arr1(19,*) * bfac
bfield = (bfieldx*bfieldx + bfieldy*bfieldy + bfieldz*bfieldz)^0.5
ufield = (bfield^2)/8./3.14159

m_B = 3.77d22 * (bfield*bfield*bfield) 
m_B = m_B / nh / nh

m_BE = 596 * (1e4/nh)^ 0.5
;m_BE = 358 * (1e4/nh)^ 0.5
m_BE = m_BE * (temp / 200) ^ 1.5

res_length = hsm*3.108e18  ;resolution length converted from pc to cm
rJfac = 15.*1.38e-16 / (4 * 3.14159 * 6.67e-8 * 1.67e-24)
res_turk =  (1./4.) * (rJfac^0.5) * (temp^0.5) / (rho^0.5)
res_turk2 = (1./64.) * (rJfac^0.5) * (temp^0.5) / (rho^0.5)
res_orion = (1./8) * (rJfac^0.5) * (temp^0.5) / (rho^0.5)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;read in second data file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
arr1=dblarr(20,arraysize2)

close,1

openr, 1, filename2+filenum2
readf,1,arr1
close,1

error_vp = arr1(0,*)
id_vp =arr1(1,*)
x_vp=arr1(2,*)
y_vp=arr1(3,*)
z_vp=arr1(4,*)
temp_vp=arr1(5,*)
nh_vp=arr1(6,*)
afieldx_vp=arr1(7,*)
afieldy_vp=arr1(8,*)
afieldz_vp=arr1(9,*);*4.e-5
gam_vp=arr1(10,*)
rad_vp=arr1(11,*)    ;dis in AU
vrad_vp=arr1(12,*)
vrot_vp=arr1(13,*)

hsm_vp=arr1(14,*) ;in pc
mass_vp=arr1(15,*);*1.e10/h
rho_vp=nh_vp*1.67e-24

bfieldx_vp = arr1(16,*)*bfac
bfieldy_vp = arr1(17,*)*bfac
bfieldz_vp = arr1(18,*)*bfac
nh_test_vp = arr1(19,*)
bfield_vp = (bfieldx_vp*bfieldx_vp + bfieldy_vp*bfieldy_vp + bfieldz_vp*bfieldz_vp)^0.5
ufield_vp = (bfield_vp^2)/8./3.14159

print, 'mass_vp =', mass_vp[10]
print, 'max nh low-res = ', max(nh_vp)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;calculate resolution mass
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
mBE= 700.e0*(temp/200.e0)^1.5/sqrt(nh/1.e4)
mRESval= (1.5*32*mass[10])*(1.e10/h)
print, 'Mres=', mRESval
mRES=replicate(mRESval, 100)


;list = where((index_array mod 10 eq 0) or (abs(error) gt 0.01))
list = where(rho gt 0 and index_array mod 1 eq 0)
print, min(temp), max(temp)
id=id(list)
x=x(list)
y=y(list)
z=z(list)
error = error(list)
nh=nh(list)
temp=temp(list)
HII=HII(list)
H2I=H2I(list)
HDI=HDI(list)
mBE=mBE(list)
gam=gam(list)
rad = rad(list)
vrad = vrad(list)
vrot = vrot(list)
hsm = hsm(list)
mass = mass(list)
rho = rho(list)
bfield = bfield(list)
ufield = ufield(list)
res_length = res_length(list)
res_turk = res_turk(list)
res_turk2 = res_turk2(list)

;hsm = hsm*206264.806    ;convert from pc to AU

list = where(rho_vp > 0 and index_array mod 1 eq 0)
id_vp=id_vp(list)
x_vp=x_vp(list)
y_vp=y_vp(list)
z_vp=z_vp(list)
error_vp = error_vp(list)
nh_vp=nh_vp(list)
temp_vp=temp_vp(list)
afieldx_vp=afieldx_vp(list)
afieldy_vp=afieldy_vp(list)
afieldz_vp=afieldz_vp(list)
gam_vp=gam_vp(list)
rad_vp = rad_vp(list)
vrad_vp = vrad_vp(list)
vrot_vp = vrot_vp(list)
hsm_vp = hsm_vp(list)
mass_vp = mass_vp(list)
rho_vp = rho_vp(list)
bfield_vp = bfield_vp(list)
ufield_vp = ufield_vp(list)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

rmax = 1e8

nmin = 1.e-1
nmax = 1.e12

size=500
size_dbl = 500.0
msink = 7.0*1.98892e33
msink_tot = 22.0*1.98892e33
G = 6.673e-8
r = dblarr(size)
m_nfw = dblarr(size)
n_est = dblarr(size)

rho0 = 1
r_s = 3.e18
conc = 20
rvir = conc * r_s

i_dbl = -1.0
for i = 0,size-1 do begin
  i_dbl = i_dbl+1.0
  r[i] = rmax^(i_dbl/size_dbl)
  r_max = 1.5e13*r[i]
  m_nfw[i] = 4. * 3.14159 * rho0 * (r_s^3) * (alog10((r_s + r_max)/r_s) -  (r_max/(r_s + r_max)))
  n_est[i] =  nmin * nmax^(i_dbl/size_dbl)
endfor

rho_est = n_est*1.67e-24
vkep_tot = dblarr(size)
vkep = (G*msink_tot)^.5
vkep = vkep/((r*1.5e13)^.5)
vkep = vkep/1.e5
rho_prof = 9.7d16/((1.496d13*r)^2)/(1.22*1.67d-24)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;estimate Turk et al. and Schleicher et al. B-fields.  Compare with equipartition.

bfac_turk = 1.e-12
Aturk = 1.78/2.
;U_turk = bfac_turk * 2.d-29 * (n_est/1.e-1)^1.78
B_turk = 1.e-20 * (n_est/n_est[0])^(Aturk)
U_turk = B_turk * B_turk / 8. / 3.14159

B_mach = 1.e-12*(n_est^(2./3.))
U_mach = (B_mach^2) / 8. / 3.14159

B_mach_jet = 2.e-10*(n_est^(2./3.))
U_mach_jet = (B_mach_jet^2) / 8. / 3.14159

B_sch = dblarr(size)
U_sch = dblarr(size)
B_sch_est = dblarr(size)
U_sch_est = dblarr(size)
l_scale = dblarr(size)
n_scale = dblarr(size)
v_edd = dblarr(size)

cs = 2.e5 ;sound speed of 2 km/s
G = 6.67e-8

rj = 3.14159 * (cs^2) / G 
rj = rj^0.5
rj = rj / ((1.e3*1.67e-24)^0.5)

l_edd =  0.1*rj
print, 'l_edd =', l_edd

b_bur = 1./2.   ;Burgers turbulence
b_kol = 1./3. ;Kolmogorov turbulence

l_scale = 25. * 1.5e13 / ((n_est / 1.e12) ^ 0.5)
v_edd = 20. * (l_scale / 200. / 3.18e18) ^ b_kol

for i = 0,size-1 do begin
  ;if n_est[i] lt 1.e3 then l_scale[i] = l_edd * ((1.e3/1.e3)^(1./3.))
  ;if n_est[i] gt 1.e3 then l_scale[i] = l_edd * ((1.e3 / n_est[i])^(1./3.))
  ;v_edd[i] = 0.1 * cs * ((l_scale[i]/l_edd)^b_bur)
endfor

t_edd = l_edd / cs  
t_ff = 1. / (G^0.5)
t_ff = t_ff / rho_est^0.5

t_max = max(t_ff)
t_elapsed = t_max - t_ff
;print, 't_elapsed =', t_elapsed

efac = v_edd / (l_scale) 
efac = efac /  (G^0.5)  ;see eqn 10 of Schleicher et al 2010
rhofac = rho_est^0.5
efac_fin = efac / rhofac

efac_fin = t_elapsed / t_edd  
efac_fin = t_elapsed * v_edd / l_scale

Bmax = 4 * 3.14159 * rho_est * cs^2
Bmax = (Bmax ^ 0.5) / (60^0.5)

for i = 0,size-1 do begin
   if n_est[i] gt 1.e2 then B_sch[i] = 1.d-7*(n_est[i]^0.5)
   if n_est[i] lt 1.e2 then B_sch[i] = 1.d-20*(n_est[i]^7.)

   B_sch_est[i] = 1.e-20 * exp(efac_fin[i])
   if Bmax[i] lt B_sch_est[i] then B_sch_est[i] = Bmax[i]

   U_sch[i] = (B_sch[i]^2) / 8. / 3.14159   ;magnetic energy density
   U_sch_est[i] = (B_sch_est[i]^2) / 8. / 3.14159 
endfor

B_frozen = 1.e-20
B_frozen = B_frozen*(n_est/n_est[0])^(.666666)

;print, 'B_sch_est = ', B_sch_est
;print, 'U_sch_est = ', U_sch_est

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;more values from Turk et al

arr1=dblarr(2,65)
close,1
openr,1,'turk_bfield.csv'
readf,1,arr1
close,1
rho_turk = arr1(0,*) 
u_turk_norm = arr1(1,*)

n_turk = rho_turk / 1.67e-24
u_turk = u_turk_norm * (rho_turk ^ (4./3.))
bfield_turk = (u_turk^0.5) * 8 * 3.14159
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;index = where(rho gt 0)
;index = where(rho gt 1.e-30 and bfield_vp eq bfield_vp)
index = where(bfield gt 0 and rho gt 1.e-24)

rhofit = rho(index)
hsmfit = hsm(index)
massfit = mass(index)
ufit = ufield(index)
bfield_fit = bfield(index)
nfit = nh(index)

;rhofit = alog10(rhofit)
;nfit = alog10(nfit)
;bfield_fit = alog10(bfield_fit)

dims=size(rhofit, /dimensions)
n=dims[0]
print, 'n = ', n

pvolume = hsmfit^3
weight = pvolume / (total(pvolume))

rho_weighted = rhofit

dims=size(rho_weighted, /dimensions)
n=dims[0]
print, 'n = ', n

lnx = alog(nfit)
lny = alog(bfield_fit)
print, 'total_lnx =', total(lnx), 'total_lny =', total(lny)

Bfit = (n*total(lnx*lny) - total(lnx)*total(lny))/(n*total(lnx^2) - (total(lnx))^2)

print, 'B coeff. for pow. law fit = ', Bfit
print, 'diff = ', (total(lny) - Bfit*total(lnx))/n

Afit = (total(lny) - Bfit*total(lnx))/n
Afit = 2.71818^Afit

print, 'A coeff. for pow. law fit = ', Afit

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;index = where(rho lt 1.e40)
index = where(bfield_vp gt 0 and rho_vp gt 1.e-24)

rhofit = rho_vp(index)
hsmfit = hsm_vp(index)
massfit = mass_vp(index)
ufit = ufield_vp(index)
bfield_fit = bfield_vp(index)
nfit = nh_vp(index)

;rhofit = alog10(rhofit)
;nfit = alog10(nfit)
;bfield_fit = alog10(bfield_fit)

dims=size(rhofit, /dimensions)
n=dims[0]
print, 'n = ', n

pvolume = hsmfit^3
weight = pvolume / (total(pvolume))

rho_weighted = rhofit

dims=size(rho_weighted, /dimensions)
n=dims[0]
print, 'n = ', n

lnx = alog(nfit)
lny = alog(bfield_fit)
print, 'total_lnx =', total(lnx), 'total_lny =', total(lny)

Bfit_vp = (n*total(lnx*lny) - total(lnx)*total(lny))/(n*total(lnx^2) - (total(lnx))^2)

print, 'B coeff. for pow. law fit = ', Bfit_vp
print, 'diff = ', (total(lny) - Bfit_vp*total(lnx))/n

Afit_vp = (total(lny) - Bfit_vp*total(lnx))/n
Afit_vp = 2.71818^Afit_vp
print, 'A coeff. for pow. law fit = ', Afit_vp

xcen = x(where(nh eq max(nh)))
ycen = y(where(nh eq max(nh)))
zcen = z(where(nh eq max(nh)))

xcen = xcen(0)
ycen = ycen(0)
zcen = zcen(0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print, 'start plotting'

;A = FINDGEN(16) * (!PI*2/16.)
;USERSYM, 0.8*COS(A),0.8*SIN(A),/FILL
set_plot,'ps'
;;device,encaps=0
device,/encapsul
device,bits_per_pixel=8,/color,/portrait
device,yoffset=8.4
!P.charsize=0.81
;!P.charsize=1.5
!P.charthick=3.0
loadct,12,/silent
;distinct_colors=10
;col=0.6*bytscl((indgen(distinct_colors) mod 16)+1)+64

print, 'start plotting 1'

;rmin = 1.e-1
rmin = 1e5
rmax = 1e8

print, 'start plotting 1a'

bmin = 1.e-22 * bfac
bmax = 1.e-8  * bfac

print, 'start plotting 1b'

rho_min = 1.e-24
rho_max = 1.e-10

umin = 1.e-0
umax = 1.e20

print, 'start plotting 1c'

rho_data_min = min(rho)
rho_data_max = max(rho)

print, 'start plotting 1d'

nmin_vp = min(nh_vp)
nmax_vp = max(nh_vp)

print, 'start plotting 1e'

nmax = 2.e8
nmax_vp = 2.e8

arrnum = 50
arrnum_doub = 50.0

print, 'start plotting 2'

indarr=indgen(1,arrnum)
nh_bin = dblarr(arrnum+1)
nh_bin_vp = dblarr(arrnum+1)

nmin_plot = 1.e-2
nmax_plot = 1.e8

nmin = 1.e-2
nmin_vp = 1.e-1

nmin_arr = nmin
nmin_arr_vp = nmin_vp
i_doub = 0.0

for i = 0,arrnum do begin
   nh_bin[i] = nmin_arr * ((nmax/nmin_arr) ^ (i_doub / arrnum_doub))
   nh_bin_vp[i] = nmin_arr_vp * ((nmax_vp/nmin_arr_vp) ^ (i_doub / arrnum_doub))
   i_doub = i_doub + 1.0
endfor

print, 'nh_bin = ', nh_bin

b_plot = dblarr(arrnum)
nh_plot = dblarr(arrnum)
hii_plot = dblarr(arrnum)

stand_dev = dblarr(arrnum)
stand_devHI = dblarr(arrnum)
stand_devLO = dblarr(arrnum)

b_plot_vp = dblarr(arrnum)
nh_plot_vp = dblarr(arrnum)

pvolume = hsmfit^3
weight = pvolume / (total(pvolume))
variable = alog10(bfield)
variable_vp = alog10(bfield_vp)

for i = 0,arrnum-1 do begin
   index = where(nh gt nh_bin[i] and nh lt nh_bin[i+1] and bfield gt 0)
   index_vp = where(nh_vp gt nh_bin_vp[i] and nh_vp lt nh_bin_vp[i+1] and bfield_vp gt 0)

   nh_plot[i] = (alog10(nh_bin[i]) + alog10(nh_bin[i+1]))/2
   nh_plot_vp[i] = (alog10(nh_bin_vp[i]) + alog10(nh_bin_vp[i+1]))/2

   weight = pvolume(index) / total(pvolume(index))

   b_plot[i] = mean(alog10(bfield(index)))
   b_plot_vp[i] = mean(alog10(bfield_vp(index_vp)))
   hii_plot[i] = mean(alog10(HII(index)))

   var_err = variable(index)
   center = mean(bfield(index))
   center = alog10(center)
   indexHI = where(var_err ne center)
   indexLO = where(var_err ne center)

   stand_dev_arr = var_err - center
   stand_dev_arrHI = var_err(indexHI) - center
   stand_dev_arrLO = var_err(indexLO) - center


   stand_dev_arr = stand_dev_arr*stand_dev_arr
   num = size(stand_dev_arr)
   num = double(num[1])
   stand_dev[i] = sqrt(mean(stand_dev_arr))

   stand_dev_arrHI = stand_dev_arrHI * stand_dev_arrHI
   num = size(stand_dev_arrHI)
   num = double(num[1])
   stand_devHI[i] = sqrt(mean(stand_dev_arrHI))

   stand_dev_arrLO = stand_dev_arrLO * stand_dev_arrLO
   num = size(stand_dev_arrLO)
   num = double(num[1])
   stand_devLO[i] = sqrt(mean(stand_dev_arrLO))

endfor

nfit = 10.^nh_plot
bfield_fit = 10.^b_plot
hii_fit = 10.^hii_plot

dims=size(nfit, /dimensions)
n=dims[0]

lnx = alog(nfit)
lny = alog(bfield_fit)
Bfit = (n*total(lnx*lny) - total(lnx)*total(lny))/(n*total(lnx^2) - (total(lnx))^2)
print, 'B coeff. for pow. law fit = ', Bfit
Afit = (total(lny) - Bfit*total(lnx))/n
Afit = 2.71818^Afit
print, 'A coeff. for pow. law fit = ', Afit

lnx = alog(nfit)
lny = alog(hii_fit)
Dfit = (n*total(lnx*lny) - total(lnx)*total(lny))/(n*total(lnx^2) - (total(lnx))^2)
Cfit = (total(lny) - Dfit*total(lnx))/n
Cfit = 2.71818^Cfit

print, 'C coeff. for pow. law fit = ', Cfit
print, 'D coeff. for pow. law fit = ', Dfit

;setup = 1
;setup = 1.5
setup = 2
;setup = 3
;setup = 3.5
;setup = 4
;setup = 5
;setup = 6 
;setup = 7
;setup = 9
;setup = 10

;device,xsize=13.0,scale_factor=1.0
;device,ysize=12.0,scale_factor=1.0
;!p.multi=0

if setup eq 1 then begin
!p.multi=[0,2,2]
!P.charsize=0.8

nmin = min(nh)
nmax = max(nh)

plot,nh,temp,xtitle='n [cm!U-3!N]', /xlog,/ylog,ytitle='T [K]',psym=3,title='',/xstyle,/ystyle,xrange=[nmin,nmax],yrange=[1.e1,1.e5]


plot,rad,nh,xtitle='r [AU]', /xlog,/ylog,ytitle='n [cm!U-3!N]',psym=3,title='',/xstyle,/ystyle, xrange=[rmin,rmax], yrange=[nmin,nmax];,  xmargin=[6,7]

plot,nh,HII,xtitle='n [cm!U-3!N]',/ylog,ytitle='f!DHII!N',psym=3,title='',/xstyle,/ystyle,/xlog,xrange=[nmin,nmax],yrange=[1.e-10,1.e-3];, ycharsize = 1.25
oplot, nfit, hii_fit, color=150, thick=3
oplot, nfit, Cfit * nfit^Dfit, color=50, thick=3

;plot,nh,H2I,xtitle='n [cm!U-3!N]',/ylog,ytitle='f!DH!l2!N',psym=3,title='',/xstyle,/ystyle,/xlog,xrange=[nmin,nmax],yrange=[1.e-10,1.e-4];, ycharsize = 1.25

plot,rad,H2I,xtitle='r [AU]',/ylog,ytitle='f!DH!l2!N',psym=3,title='',/xstyle,/ystyle,/xlog,xrange=[rmin,rmax],yrange=[1.e-6,1.e0];, ycharsize = 1.25


xyouts, 0.42, 0.92, '(a)', charsize=0.8, /normal
xyouts, 0.92, 0.92, '(b)', charsize=0.8, /normal
xyouts, 0.42, 0.12, '(c)', charsize=0.8, /normal
xyouts, 0.92, 0.12, '(d)', charsize=0.8, /normal
endif

if setup eq 1.5 then begin

nmin = 1.e-2

imin  = alog10(nmin)
imax  = alog10(nmax)

lfac = 5.0

len = fix(imax) - fix(imin) + 1
len = len*lfac

nbin = make_array(len, /double, value=0)
print, nbin

for i = 0,len-1 do begin
  expo = float(imin) + float(i)/lfac
  nbin[i] = 10.^expo
endfor
print, 'nbin = ', nbin

m_B_bin = make_array(len, /double, value=1e-10)
m_BE_bin = make_array(len, /double, value=1e-10)

for i = 0,len-2 do begin
   index = where(nh gt nbin[i] and nh lt nbin[i+1])
   m_B_bin[i] = mean(m_B[index])
   m_BE_bin[i] = mean(m_BE[index])
endfor
print, 'm_B_bin = ', m_B_bin
print, 'm_BE_bin = ', m_BE_bin

print, 'm_b', min(m_B), max(m_B)

print, 'size(nbin) = ', size(nbin)

s = size(nbin)
s = s[1]

nbin = nbin[0: s-7]
m_B_bin = m_B_bin[0: s-7]
m_BE_bin = m_BE_bin[0: s-7]

plot,nh, m_B, xtitle='!6n!DH!N [cm!U-3!N]', ytitle='!6Mass [M!I!9n!6!N]', /xlog,/ylog,psym=3,title='',/xstyle,/ystyle,xrange=[nmin,nmax],yrange=[1.e-6,1.e6]
oplot, nh, m_BE, psym=3, color=100

oplot, nbin, m_B_bin, thick=5.0
oplot, nbin, m_BE_bin, thick=5.0

xa = 1.e7
y1 = 1.e4
y2 = 1.e5 

x1 = 2.e5
x2 = 9.e6

xyouts, xa, y1, 'M!DB!N', charsize=1.0
xyouts, xa, y2, 'M!DBE!N',  charsize=1.0

plots,[x1,x2],[y1,y1]          , thick=4.0
plots,[x1,x2],[y2,y2], color=100, thick=4.0

endif


if setup eq 2 then begin
;device,xsize=13.0,scale_factor=1.0
;device,ysize=6.0,scale_factor=1.0
;!p.multi=[0,2,1]
!p.charsize = 1.5
!p.charthick = 2.0

n_th = [1, 2, 4, 8, 32, 128, 256, 512, 2048, 1.00E+04, 1.00E+06, 1.00E+08, 1.00E+09]

ell = [0, 0.235, 0.348, 0.4466, 0.627, 0.796, 0.877, 0.956, 1.105, 1.262, 1.629, 1.884, 1.979]

B_th = 3.6e-11 * (n_th^(2./3.)) * exp(58. * 0.066 * ell)

B_ex = B_th

len = 13
for i = 0,len-1 do begin
  if n_th[i] gt 0.01 and n_th[i] le 1.0 then B_ex[i] = 3.6e-11 * (n_th[i]^(2./3.))
  if n_th[i] ge 1.0 and n_th[i] le 250 then B_ex[i] = 3.6e-11 * (n_th[i]^(1.43))
  if n_th[i] ge 250 and n_th[i] le 1.e7 then B_ex[i] = 7.3e-10 * (n_th[i]^(0.87))
  if n_th[i] ge 1.e7 and n_th[i] le 1.e8 then B_ex[i] = 1.7e-8 * (n_th[i]^(2./3.))
endfor

print, 'n_th = ', n_th
print, 'B_th = ', B_th
print, 'B_ex = ', B_ex

print, 'start plotting 3'

index = where(indarr mod 1 eq 0)
iplot = where((index_array mod 8 eq 0) or (nh gt 10 and nh lt 100 and index_array mod 4 eq 0) or (nh gt 100 and nh lt 1000 and index_array mod 2 eq 0) or (nh gt 1000.e0))
;iplot = where((index_array mod 16 eq 0))

plot, nh[iplot], bfield[iplot], xtitle='!6n!DH!N [cm!U-3!N]', psym=3, /xlog, /ylog, /xstyle, /ystyle, xrange = [nmin_plot, nmax_plot], yrange = [1e-13, 1e1], ytitle = '!6B [Gauss]'
;oplot, nfit, Afit * nfit^Bfit, color=150, thick=3 

;oplot, 10^nh_plot, 10^b_plot,  color=100, linestyle=3, thick = 10.0
;errplot, 10^nh_plot[index], 10^(b_plot[index] - stand_devLO[index]), 10^(b_plot[index] + stand_devHI[index]), linestyle=3, color=110,thick=3.0

oplot, n_th, B_th, linestyle=3, color=100,thick=10.0 

oplot, n_th, B_ex, linestyle=1, color=150,thick=10.0

oplot, nfit, 5.e-21 * bfac * (nfit/nfit[0])^(2./3.), color=20, thick=5
;oplot, nfit, Afit * nfit^Aturk, linestyle=2, thick=3.0, color=20
oplot, n_turk, bfield_turk*2.e-6 * bfac, linestyle=2, color=200, thick=8.0
;oplot, 10^nh_plot_vp, 10^b_plot_vp,  linestyle=1, color=150, thick=10.0

print, 'n = ', (10^nh_plot) * 2e-24
print, 'bfit = ', 10^b_plot

x1  = 5.e-2 
x2 =  1.e0 
xa =  2.e0

y1 = 1.e-9 * bfac
y2 = 1.e-10 * bfac
y3 = 1.e-11 * bfac
y4 = 1.e-12 * bfac

csize = 1.2
xyouts, xa, y1, 'theory prediction', charsize=csize
xyouts, xa, y2, 'fit to simulation',  charsize=csize
xyouts, xa, y3, 'n!U2/3!N', charsize=csize
xyouts, xa, y4, 'Turk ea 2011 (x 10!U3!N)', charsize=csize

plots,[x1,x2],[y1,y1], color=100, thick=10.0, linestyle=3 
plots,[x1,x2],[y2,y2], color=150, thick=10.0, linestyle=1
plots,[x1,x2],[y3,y3], color=20, thick=5.0
plots,[x1,x2],[y4,y4], color=200, thick=5.0, linestyle=2

print, 'Afit =', Afit, ' Bfit =', Bfit, ' Bcen =', Afit * 1.e10^Bfit
print, 'Afit_vp =', Afit_vp, ' Bbit_vp =', Bfit_vp, ' Bcen_vp =', Afit_vp * 1.e10^Bfit_vp 

;;;;;;;;;DivB calculation;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;divb_err = hsm * abs(divb) / bfield
;divb_err = divb_err / (1.e0/h) * zp1 / 1000. ;convert from phys. pc to comov. kpc

;emin = 1.e-5
;emax = 1.e1

;plot, nh[iplot], divb_err[iplot], xtitle='!6n [cm!U-3!N]', ytitle = '!9G .!6!N B error ', title = '', psym=3, /xlog, /ylog, /xstyle, /ystyle, xrange = [nmin_plot, nmax], yrange = [emin, emax]

;arrnum = 10
;arrnum_doub = 10.0

;indarr=indgen(1,arrnum)
;nh_bin = dblarr(arrnum+1)

;nmin = 1.e-2

;nmin_arr = nmin
;i_doub = 0.0

;for i = 0,arrnum do begin
;   nh_bin[i] = nmin_arr * ((nmax/nmin_arr) ^ (i_doub / arrnum_doub))
;   i_doub = i_doub + 1.0
;endfor

;nh_plot = dblarr(arrnum)
;divb_plot = dblarr(arrnum)

;for i = 0,arrnum-1 do begin
;   index = where(nh gt nh_bin[i] and nh lt nh_bin[i+1] and bfield gt 0)
;   nh_plot[i] = (alog10(nh_bin[i]) + alog10(nh_bin[i+1]))/2
;   divb_plot[i] = mean(alog10(divb_err(index)))
;endfor

;oplot, 10^nh_plot, 10^divb_plot,  thick = 5.0, color=100, linestyle=3

endif


if setup eq 3 then begin

ymin = min (bfield / (nh ^ 0.6666667))
ymax = max (bfield / (nh ^ 0.6666667))

plaw = 2.0 / 3.0

plot, nh, bfield / (nh ^ plaw), xtitle='!6n [cm!U-3!N]', ytitle = 'B [Gauss] / n!U2/3!N', psym=3, /xlog, /ylog, /xstyle, /ystyle, xrange = [nmin, nmax], yrange = [ymin, ymax], xmargin = [12,2]
oplot, n_est, Afit * n_est^(Bfit - plaw), color=200, thick=3
oplot, 10^nh_plot, 10^b_plot / ((10^nh_plot) ^ plaw),  color=100, linestyle=3, thick=3.0
oplot, n_est, B_frozen / (n_est ^ plaw), color=200, thick=3
oplot, n_est, B_turk / (n_est ^ plaw), linestyle=2, thick=3.0

print, 'nh_plot = ', nh_plot, 'b_plot = ', b_plot
print, 'Afit =', Afit, ' Bfit =', Bfit

endif

if setup eq 3.5 then begin

plot, nh, bfield, xtitle='!6n [cm!U-3!N]', ytitle = 'B [Gauss]', psym=3, /xlog, /ylog, /xstyle, /ystyle, xrange = [1.e-4, 1.e4], yrange = [1e-14, 1e-8], xmargin = [12,2]

endif

if setup eq 4 then begin

arrnum = 5
arrnum_doub = 5.0
rho_bin = dblarr(arrnum+1)

i_doub = 0.0
for i = 0,arrnum do begin
rho_bin[i] = rho_data_min * ((rho_data_max/rho_data_min) ^ (i_doub / arrnum_doub))
i_doub = i_doub + 1.0
endfor


u_plot = dblarr(arrnum)
u_plot2 = dblarr(arrnum)
rho_plot = dblarr(arrnum)
stand_dev = dblarr(arrnum)
stand_devHI = dblarr(arrnum)
stand_devLO = dblarr(arrnum)

pvolume = hsmfit^3
weight = pvolume / (total(pvolume))
variable = alog10(ufield/(rho^1.333333333))

for i = 0,arrnum-1 do begin
index = where(rho gt rho_bin[i] and rho lt rho_bin[i+1] and ufield gt 0)

rho_plot[i] = (alog10(rho_bin[i]) + alog10(rho_bin[i+1]))/2

weight = pvolume(index) / total(pvolume(index))

u_plot[i] = total(ufield(index)*weight)
u_plot[i] = mean(alog10(ufield(index)))

print, 'rho_plot = ', rho_plot[i], ' u_plot = ', u_plot[i]

var_err = variable(index)
rho_dum = 10^rho_plot[i]
rho_dum = rho_dum^1.3333333
center =  (10^u_plot[i]) / rho_dum
center = alog10(center)
indexHI = where(var_err gt center)
indexLO = where(var_err lt center)

print, 'center =', center 
print, ' var_err = ', var_err

stand_dev_arr = var_err - center
stand_dev_arrHI = var_err(indexHI) - center
stand_dev_arrLO = var_err(indexLO) - center

stand_dev_arr = stand_dev_arr*stand_dev_arr
num = size(stand_dev_arr)
num = double(num[1])
stand_dev[i] = sqrt(mean(stand_dev_arr))
print, 'dev= ', stand_dev[i], ' num = ', num
;stand_dev[i] = stand_dev[i] / sqrt(num)

stand_dev_arrHI = stand_dev_arrHI * stand_dev_arrHI
num = size(stand_dev_arrHI)
num = double(num[1])
stand_devHI[i] = sqrt(mean(stand_dev_arrHI))
print, 'devHI= ', stand_devHI[i], ' num = ', num
;stand_devHI[i] = stand_devHI[i] / sqrt(num)

stand_dev_arrLO = stand_dev_arrLO * stand_dev_arrLO
num = size(stand_dev_arrLO)
num = double(num[1])
stand_devLO[i] = sqrt(mean(stand_dev_arrLO))
print, 'devLO= ', stand_devLO[i], ' num = ', num
;stand_devLO[i] = stand_devLO[i] / sqrt(num)

endfor


plot, rho, ufield/(rho^1.333333333), xtitle='!7q!6 [g cm!U-3!N]', ytitle = 'U!DB!N / !7q!6!U4/3!N', psym=3,  /xlog, /ylog, /xstyle, /ystyle, xrange = [rho_min, rho_max], yrange = [umin, umax]

oplot, 10^rho_plot, 10^(u_plot - rho_plot*1.3333333), linestyle=0, thick=3.0, color=100
errplot, 10^rho_plot, 10^(u_plot - rho_plot*1.3333333 - stand_devLO), 10^(u_plot - rho_plot*1.3333333 + stand_devHI), thick = 2.0, color=100 


oplot, rho_est, U_turk/(rho_est^1.3333333), linestyle=2, thick=3.0

oplot, rho_est, U_sch/(rho_est^1.3333333), linestyle=1, thick=3.0

endif


if setup eq 5 then begin

lmin = min(res_length)
lmax = max(res_length)

nmin = 1.e-1
nmax = 1.e8

iplot = where((index_array mod 8 eq 0) or (nh gt 10 and nh lt 100 and index_array mod 4 eq 0) or (nh gt 100 and nh lt 1000 and index_array mod 2 eq 0) or (nh gt 1000.e0))

plot, nh[iplot], res_length[iplot], xtitle='!6n!DH!N [cm!U-3!N]', ytitle = '!6resolution length [cm]', psym=3,  /xlog, /ylog, /xstyle, /ystyle, xrange = [nmin, nmax], yrange = [lmin, lmax]
oplot, nh[iplot], res_turk[iplot], psym=3, color=100
oplot, nh[iplot], res_turk2[iplot], psym=3, color=200
oplot, nh[iplot], res_orion[iplot], psym=3, color=50

x1  = 5.e3
x2 =  5.e4
xa =  7.e4

y1 = 4.e19 
y2 = 2.e19 
y3 = 1.e19 

xyouts, xa, y1, 'Gadget', charsize=0.8
xyouts, xa, y2, 'Turk ea,  4 cells',  charsize=0.8
xyouts, xa, y3, 'Turk ea, 64 cells', charsize=0.8

plots,[x1,x2],[y1,y1]          , thick=3.0
plots,[x1,x2],[y2,y2], color=100, thick=3.0
plots,[x1,x2],[y3,y3], color=200, thick=3.0

endif


if setup eq 6 then begin

!p.multi=[0,2,2]

amin = 1.e-22
amax =  1.e-6

print, 'min ax =', min(afieldx_vp)
print, 'max ax =', max(afieldx_vp)

afieldx_vp = abs(afieldx_vp)
afieldy_vp = abs(afieldy_vp)
afieldz_vp = abs(afieldz_vp)

plot, nh_vp, bfield_vp, xtitle='n [cm!U-3!N]', ytitle = 'B [Gauss]', title = 'B from Vec. Pot.', psym=3, /xlog, /ylog, /xstyle, /ystyle, xrange = [nmin, nmax], yrange = [bmin, bmax]

plot, nh_vp, afieldx_vp, xtitle='n [cm!U-3!N]', ytitle = 'Ax', title = ' ', psym=3, /xlog, /ylog, /xstyle, /ystyle, xrange = [nmin, nmax], yrange = [amin, amax]

plot, nh_vp, afieldy_vp, xtitle='n [cm!U-3!N]', ytitle = 'Ay', title = ' ', psym=3, /xlog, /ylog, /xstyle, /ystyle, xrange = [nmin, nmax], yrange = [amin, amax]

plot, nh_vp, afieldz_vp, xtitle='n [cm!U-3!N]', ytitle = 'Az', title = ' ', psym=3, /xlog, /ylog, /xstyle, /ystyle, xrange = [nmin, nmax], yrange = [amin, amax]

endif


if setup eq 7 then begin

!P.charsize=1.5

divb_err = hsm * abs(divb) / bfield
divb_err = divb_err / (1.e0/h) * zp1 / 1000. ;convert from phys. pc to comov. kpc

emin = 1.e-5
emax = 1.e1

plot, nh, divb_err, xtitle='!6n [cm!U-3!N]', ytitle = '!9G .!6!N B error', title = '', psym=3, /xlog, /ylog, /xstyle, /ystyle, xrange = [nmin, nmax], yrange = [emin, emax], xmargin = [15,3]

arrnum = 10
arrnum_doub = 10.0

indarr=indgen(1,arrnum)
nh_bin = dblarr(arrnum+1)

nmin = 1.e-2

nmin_arr = nmin
i_doub = 0.0

for i = 0,arrnum do begin
   nh_bin[i] = nmin_arr * ((nmax/nmin_arr) ^ (i_doub / arrnum_doub))
   i_doub = i_doub + 1.0
endfor

nh_plot = dblarr(arrnum)
divb_plot = dblarr(arrnum)

for i = 0,arrnum-1 do begin
   index = where(nh gt nh_bin[i] and nh lt nh_bin[i+1] and bfield gt 0)
   nh_plot[i] = (alog10(nh_bin[i]) + alog10(nh_bin[i+1]))/2
   divb_plot[i] = mean(alog10(divb_err(index)))
endfor

oplot, 10^nh_plot, 10^divb_plot,  thick = 10.0, color=200, linestyle=3

endif


if setup eq 9 then begin

device,xsize=13.0,scale_factor=1.0
device,ysize=6.0,scale_factor=1.0
!p.multi=[0,2,1]

;plot, rad, vrad, xtitle='r [AU]', psym=3,/xstyle,/ystyle, /xlog, ytitle='v!Drad!N [km/s]', xrange=[rmin,rmax], yrange=[-50, 50],  ycharsize = 1.2

;plot, rad, vrot, xtitle='r [AU]', psym=3,/xstyle,/ystyle, /xlog, ytitle='v!Drot!N [km/s]', xrange=[rmin,rmax], yrange=[-10, 50],  ycharsize = 1.2

plot, nh, vrad, xtitle='n [cm!U-3!N]', psym=3,/xstyle,/ystyle, /xlog, ytitle='v!Drad!N [km/s]', xrange=[nmin,nmax], yrange=[-50, 50],  ycharsize = 1.2

plot, nh, vrot, xtitle='n [cm!U-3!N]', psym=3,/xstyle,/ystyle, /xlog, ytitle='v!Drot!N [km/s]', xrange=[nmin,nmax], yrange=[-10, 50],  ycharsize = 1.2

endif

if setup eq 10 then begin

rotate=0
if rotate eq 1 then begin
  ycopy = y
  zcopy = z
  y = zcopy
  z = ycopy
endif

;xmin = min(x)
;xmax = max(x)
;ymin = min(y)
;ymax = max(y)
;zmin = min(z)
;zmax = max(z) 

xmin = -0.5 & xmax = 0.5
ymin = -0.5 & ymax = 0.5
zmin = -0.5 & zmax = 0.5

;xmin = 1.5*min(x(ind_hiB))
;xmax = 1.5*max(x(ind_hiB))
;ymin = 1.5*min(y(ind_hiB))
;ymax = 1.5*max(y(ind_hiB))

print, 'xmin = ', xmin, 'xmax = ', xmax
print, 'ymin = ', ymin, 'ymax = ', ymax
print, 'zmin = ', zmin, 'zmax = ', zmax
print, 'x_length =', size(x)

zdis = zmax - zmin

xcen = x(where(nh eq max(nh)))
ycen = y(where(nh eq max(nh)))
zcen = z(where(nh eq max(nh)))

xcen = xcen(0)
ycen = ycen(0)
zcen = zcen(0) 

zmin = zcen - 2.0
zmax = zcen + 2.0
print, 'zcen = ', zcen, 'zmin = ', zmin, 'zmax =', zmax

plot, x, y, psym=3, xtitle = 'x [pc]', ytitle = 'y [pc]', xrange = [xmin, xmax], yrange = [ymin, ymax], /xstyle, /ystyle

bfit = Afit * nh^Bfit

xplot = x(ind_hiB)
yplot = y(ind_hiB)
zplot = z(ind_hiB)
hplot = hsm(ind_hiB)
nhplot = nh(ind_hiB)
bplot = bfield(ind_hiB)
bfitplot = bfit(ind_hiB)

arraysize = size(xplot)
arraysize = arraysize(1)
print, 'arraysize=',  arraysize
oplot, xplot, yplot, psym=8, color=100

phi = findgen(90)/89.*2*!pi
circle_x = dblarr(90)
circle_y = dblarr(90)
for i = 0, arraysize-1 do begin
 ;print, nhplot(i), hplot(i), bplot(i), bfitplot(i)
 ;arrow, xplot(i), yplot(i), xplot(i)+hplot(i), yplot(i)+hplot(i), /data, hsize=50.0, color=100
 for j =0, 89 do begin
   circle_x[j] = xplot(i) + hplot(i)*sin(phi[j])
   circle_y[j] = yplot(i) + hplot(i)*cos(phi[j])
 endfor
 oplot, circle_x, circle_y, thick=2.0, color=100 
endfor

endif

;USERSYM, 0.3*COS(A),0.3*SIN(A)

!P.charsize=1.5
!P.charsize=1.
!P.charthick=1.0
device,/close_file
set_plot,'x'
end

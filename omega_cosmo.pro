pro omega_cosmo
close,2

arraysize = 761  ;ref4

index_array=indgen(1,arraysize)

filenum1 = '0537'

dir = '/nobackupp12/astacy/'

filename =  dir +'omega_zoom10_ref4_'

;;;;;;;;;;;;;;;;;;;;;;;;
;read in redshifts
;;;;;;;;;;;;;;;;;;;;;;;
Tcmb = 2.7e0 * 22.0
h=0.7
G = 6.67e-8

set_plot,'ps'
;;device,encaps=0
device,/encapsul
device,bits_per_pixel=8,/color,/portrait
device,yoffset=8.4
!P.charsize=0.81
!P.charthick=3.0
loadct,12,/silent

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;read in main data file
;;;;;;;;;;;;;;;;;;;;;;;;;;;
arr1=dblarr(8,arraysize)

close,1
openr, 1, filename+filenum1
readf,1,arr1
close,1

zp1=27

time = arr1(0,*)
id =arr1(1,*)
nh=arr1(2,*)
den=arr1(3,*)
hsm=arr1(4,*)
delh=arr1(5,*)
m_del_k=arr1(6,*)
omega=arr1(7,*)

nmin = 1.e-1
nmax = 1.e8

plot,nh, omega, xtitle='!6n!DH!N [cm!U-3!N]', ytitle='Omega', /xlog,psym=4,title='',/xstyle,/ystyle,xrange=[nmin,nmax],yrange=[0.5,1.2]

!P.charsize=1.5
!P.charsize=1.
!P.charthick=1.0
device,/close_file
set_plot,'x'
end



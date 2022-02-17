pro divbcom_cosmo

orient=1
which_sim=5

if orient eq 1 then endtype = '.divb'
if orient eq 2 then endtype = '.divbxz'

base = '/scratch/00863/minerva/'

if which_sim eq 1 then sim_name='mono_'
if which_sim eq 2 then sim_name='mono_corr_'
if which_sim eq 3 then sim_name='mono_corr2_'
if which_sim eq 4 then sim_name='mono3d_'
if which_sim eq 5 then sim_name='mono3d_corr2_'


i  = 10
i2 = 40
i3 = 70

;ngrid = 200
ngrid=256

array1=dblarr(ngrid,ngrid)
array2=dblarr(ngrid,ngrid)
array3=dblarr(1,ngrid)
redshift1=dblarr(1,1)
time1=dblarr(1,1)
entries1=dblarr(1,1)
close,1
  openr,1, base + sim_name + string(i,format='(i4.4)')+endtype
readu,1,array1
readu,1,array2
readu,1,array3
readu,1,redshift1
readu,1,time1
readu,1,entries1

array1=rotate(array1,1)
array1=reverse(array1,1)

array2=rotate(array2,1)
array2=reverse(array2,1)

if entries1 gt 0 then begin

array4=dblarr(1,entries1)
array5=dblarr(1,entries1)
array6=dblarr(1,entries1)

readu,1,array4
readu,1,array5
readu,1,array6

endif

close,1

array7=dblarr(ngrid,ngrid)
array8=dblarr(ngrid,ngrid)
array9=dblarr(1,ngrid)
redshift2=dblarr(1,1)
time2=dblarr(1,1)
entries2=dblarr(1,1)
openr,1, base + sim_name + string(i2,format='(i4.4)')+endtype
readu,1,array7
readu,1,array8
readu,1,array9
readu,1,redshift2
readu,1,time2
readu,1,entries2

array7=rotate(array7,1)
array7=reverse(array7,1)

array8=rotate(array8,1)
array8=reverse(array8,1)

if entries2 gt 0 then begin

array10=dblarr(1,entries2)
array11=dblarr(1,entries2)
array12=dblarr(1,entries2)

readu,1,array10
readu,1,array11
readu,1,array12

endif

close,1

array13=dblarr(ngrid,ngrid)
array14=dblarr(ngrid,ngrid)
array15=dblarr(1,ngrid)
redshift3=dblarr(1,1)
time3=dblarr(1,1)
entries3=dblarr(1,1)
openr,1, base + sim_name + string(i3,format='(i4.4)')+endtype
readu,1,array13
readu,1,array14
readu,1,array15
readu,1,redshift3
readu,1,time3
readu,1,entries3

array13=rotate(array13,1)
array13=reverse(array13,1)

array14=rotate(array14,1)
array14=reverse(array14,1)

print, 'entries3=', entries3

if entries3 gt 0 then begin

array16=dblarr(1,entries3)
array17=dblarr(1,entries3)
array18=dblarr(1,entries3)

readu,1,array16
readu,1,array17
readu,1,array18

endif

close,1

print, 'entries1 = ', entries1, 'entries2 =', entries2, 'entries3=', entries3
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;get data to overplot the density contours

array19=dblarr(ngrid,ngrid)
array20=dblarr(ngrid,ngrid)
array21=dblarr(1,ngrid)
close,1

openr,1, base + sim_name + string(i,format='(i4.4)')+endtype
readu,1,array19
readu,1,array20
readu,1,array21

array19=rotate(array19,1)
array19=reverse(array19,1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;get data to overplot the density contours

array22=dblarr(ngrid,ngrid)
array23=dblarr(ngrid,ngrid)
array24=dblarr(1,ngrid)
close,1

openr,1, base + sim_name + string(i2,format='(i4.4)')+endtype
readu,1,array22
readu,1,array23
readu,1,array24

array22=rotate(array22,1)
array22=reverse(array22,1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;get data to overplot the density contours

array25=dblarr(ngrid,ngrid)
array26=dblarr(ngrid,ngrid)
array27=dblarr(1,ngrid)
close,1

openr,1, base + sim_name + string(i3,format='(i4.4)')+endtype
readu,1,array25
readu,1,array26
readu,1,array27

array25=rotate(array25,1)
array25=reverse(array25,1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

set_plot,'ps'
device,encapsulated=1,bits_per_pixel=8,/color,/portrait,xsize=20.0,ysize=7.
!P.charsize=0.7
!P.charthick=2.0
!x.tickinterval=700
!y.tickinterval=700
!P.multi = [0,3,1,0]

loadct,4
;loadct,3

contour, array1, levels=[9,10,11], color=250, /nodata, /normal, xstyle=1,ystyle=1, position=[0, 0.1, 0.224*(42.0/32.0), 0.84+0.1]
loadct, 4
tvscl,array1,0.0,0.1,xsize=0.224*(42.0/32.0),ysize=0.84,/normal
loadct, 0
;contour, array19,/overplot, levels=[7.5, 8., 8.5, 9.], c_colors=[150,150,150,150], /normal, xstyle=1,ystyle=1, thick=3.

contour, array1, levels=[9,10,11], color=250, /nodata, /normal, xstyle=1,ystyle=1, position=[0.3, 0.1, 0.3+0.224*(42.0/32.0), 0.84+0.1]
loadct, 4
tvscl,array7,0.3,0.1,xsize=0.224*(42.0/32.0),ysize=0.84,/normal
loadct, 0
;contour, array22,/overplot, levels=[7.5, 8, 8.5,9], c_colors=[150,150,150,150], /normal, xstyle=1,ystyle=1, thick=3., color=500


contour, array1, levels=[9,10,11], color=250, /nodata, /normal, xstyle=1,ystyle=1, position=[0.6, 0.1, 0.6+0.224*(42.0/32.0), 0.84+0.1], thick=5.0
loadct, 4
tvscl,array13,0.6,0.1,xsize=0.224*(42.0/32.0),ysize=0.84,/normal
loadct, 0
;contour, array25,/overplot, levels=[7.5,8.,8.5,9], c_colors=[150,150,150,150], /normal, xstyle=1,ystyle=1, thick=3., color=500

loadct, 4
tvscl,array3,0.93,0.2,xsize=0.02,ysize=0.5,/normal

loadct,0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
loadct, 12
ssize  = 0.9
ssize2 = 0.6
scol = 100

imax = entries1(0,0)-1
sbegin = 2

if entries1 gt 0 then begin

  if entries1 gt 1 then begin 
    plots,0.224*(42.0/32.0)*array4(0,1), 0.1+0.84*array5(0,1),psym=1,/normal,color=scol,symsize=ssize,thick=2.5
    for i=sbegin,entries1(0,0)-1 do begin
      plots,0.224*(42.0/32.0)*array4(0,i), 0.1 + 0.84*array5(0,i),psym=4,/normal,color=0, symsize=ssize2
    endfor
    ;plots,0.224*(42.0/32.0)*array4(0,imax), 0.1 + 0.84*array5(0,imax),psym=4,/normal,color=scol, symsize=ssize
  endif

  plots,0.224*(42.0/32.0)*array4(0,0), 0.1+0.84*array5(0,0),psym=2,/normal,color=scol, symsize=ssize, thick=2.0

endif

if entries2 gt 0 then begin

  if entries2 gt 1 then begin 
    plots,0.3 + 0.224*(42.0/32.0)*array10(0,1), 0.1+0.84*array11(0,1),psym=1,/normal,color=scol, symsize=ssise, thick=2.5
    for i=sbegin,entries2(0,0)-1 do begin
      plots,0.3 + 0.224*(42.0/32.0)*array10(0,i), 0.1 + 0.84*array11(0,i),psym=4,/normal,color=0, symsize=ssize2
    endfor
  endif

  plots,0.3+0.224*(42.0/32.0)*array10(0,0), 0.1+0.84*array11(0,0),psym=2,/normal,color=scol, symsize=ssize, thick=2.0

endif


if entries3 gt 0 then begin
 
  if entries3 gt 1 then begin 
    plots,0.6+0.224*(42.0/32.0)*array16(0,1), 0.1 + 0.84*array17(0,1),psym=1,/normal,color=scol, symsize=ssize, thick=2.5
    for i=sbegin,entries3(0,0)-1 do begin
      plots,0.6+0.224*(42.0/32.0)*array16(0,i), 0.1 + 0.84*array17(0,i),psym=4,/normal,color=0, symsize=ssize2
    endfor
  endif

  plots,0.6+0.224*(42.0/32.0)*array16(0,0), 0.1 + 0.84*array17(0,0),psym=2,/normal,color=scol, symsize=ssize, thick=2.0

endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

plots,[0.0,1.0],[1.0,1.0],/normal
plots,[0.0,1.0],[0.1,0.1],/normal
plots,[0.0,1.0],[0.94,0.94],/normal
plots,[0.999,0.999],[0.1,0.94],/normal

;xyouts, 0.90, 0.86, '1000-5000 yr', /normal

xyouts,0.90,0.78, 'log error DivB', /normal

xyouts,0.955,0.68,'1.0',/normal & xyouts,0.955,0.44,'-0.5',/normal & xyouts,0.955,0.2,'-2.0',/normal

;xyouts,0.880,0.065,'Size: 10!U4!N AU',/normal

xyouts, 0.1, 0.85, 't = 0.1', color=500, /normal, charthick=3.0
xyouts, 0.4, 0.85, 't = 0.4', color=500, /normal, charthick=3.0
xyouts, 0.7, 0.85, 't = 0.7', color=500, /normal, charthick=3.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


device,/close_file
set_plot,'x'

end

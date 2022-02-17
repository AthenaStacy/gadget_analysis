pro vel_cosmo

;dir='/work/00863/minerva/binHR/'
;dir='/work/00863/minerva/bin_zoom_cut/'
;dir='/work/00863/minerva/bin_zoom/'
dir='/scratch/00863/minerva/'

;name='bin_HR10_ideal_'
;name='bin_zoom10_new_cut_'
;name='bin_zoom10_new_cut_ref3_'
name='bin_zoom10_ref4_'

;which_form = '(i3.3)'
which_form = '(i4.4)'

orient = 1
if orient eq 1 then endtype = '.vel'
if orient eq 2 then endtype = '.velxz'

ngrid=200

N_begin = 445
;N_begin = 679
;N_begin = 5979
;N_begin = 7010

for i=N_begin,N_begin do begin

array1=dblarr(ngrid,ngrid)
array2=dblarr(ngrid,ngrid)
array3=dblarr(1,ngrid)
redshift=dblarr(1,1)
time=dblarr(1,1)
entries=dblarr(1,1)
close,1
openr,1, dir+name+string(i,format=which_form)+endtype
readu,1,array1
readu,1,array2
readu,1,array3
readu,1,redshift
readu,1,time
readu,1,entries

print, time
time0 = time*1.e6

endfor

N_begin=N_begin
N_end=N_begin
jpg_num = N_begin

for i=N_begin,N_end,10 do begin

array1=dblarr(ngrid,ngrid)
array2=dblarr(ngrid,ngrid)
array3=dblarr(1,ngrid)
redshift=dblarr(1,1)
time=dblarr(1,1)
entries=dblarr(1,1)
close,1

openr,1, dir+name+string(i,format=which_form)+endtype

readu,1,array1
readu,1,array2
readu,1,array3
readu,1,redshift
readu,1,time
readu,1,entries

print, time
time = time*1.e6 - time0
print, time
print, entries
print, 'N =', i
print, 'jpg_num =', jpg_num

array1=rotate(array1,1)
array1=reverse(array1,1)

array2=rotate(array2,1)
array2=reverse(array2,1)

if entries gt 0 then begin

array4=dblarr(1,entries)
array5=dblarr(1,entries)
array6=dblarr(1,entries)

readu,1,array4
readu,1,array5
readu,1,array6

endif

close,1

a=findgen(17)*(!pi*2.0/16.0)
usersym,0.8*cos(a),0.8*sin(a),/fill

;array3=congrid(array3,15,180)
;array3=congrid(array3,10,128)
array3=congrid(array3,5,100)

loadct,4
;loadct, 3

set_plot,'x'
device,decomposed=0
window,xsize=ngrid,ysize=ngrid,retain=2
;window,xsize=700,ysize=600,retain=2

tvscl,array1
tvscl,array3,0.90,0.7,/normal

!P.charsize=1.4

loadct,0

if orient eq 1 then begin
 if entries gt 0 then begin
 plots,array4(0,0),array5(0,0),psym=2,/normal, color=0
 endif

 if entries gt 1 then begin
  for j=1,entries(0,0)-1 do begin
  plots,array4(0,j),array5(0,j),psym=4,/normal, color=0
 endfor
 endif
endif


if orient eq 2 then begin
 if entries gt 0 then begin
 plots,array4(0,0),array6(0,0),psym=2,/normal, color=0
 endif

; if entries gt 1 then begin
; plots,array4(0,1),array6(0,1),psym=1,/normal, color=0
; endif

 if entries gt 1 then begin
  for j=1,entries(0,0)-1 do begin
  plots,array4(0,j),array6(0,j),psym=4,/normal, color=0
  endfor
 endif
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

charsize=0.3

;xyouts,0.03,0.95,'z = '+string(redshift,format='(d7.4)'), charsize=csize, /normal
;xyouts,0.03,0.9,'time = '+string(time,format='(d10.4)')+' yr', charsize=csize, /normal

;xyouts,0.03,0.04,'Length: 1 pc (physical)', charsize=csize, /normal

;xyouts,0.70,0.82,'log B-mag [G]',charsize=csize, /normal


;xyouts,0.93,0.94,'-15', charsize=csize, /normal
;xyouts,0.93,0.70,'-12', charsize=csize, /normal

;arrow,0.90,0.7,0.90,0.955,/normalized,hsize=0
;arrow,0.90,0.955,0.92,0.955,/normalized,hsize=0
;arrow,0.92,0.955,0.92,0.7,/normalized,hsize=0
;arrow,0.92,0.7,0.90,0.7,/normalized,hsize=0

!x.margin=0 
!y.margin=0

image=tvrd(true=1)
write_jpeg,name+string(jpg_num,format='(i4.4)')+endtype+'.jpg',image, TRUE=1, QUALITY=75

set_plot,'x'

jpg_num = jpg_num+1

endfor

end

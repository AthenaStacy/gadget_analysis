pro create_ic_nfw, file, outfile, npart0, MC=MC, ICs = ICs, debug = debug, setup = setup, coremc = coremc, rescaleonly = rescaleonly, nocore= nocore, rin = rin



;;@gadget_data

; ==========================================================================
; Read glass file
; ==========================================================================


; ==========================================================================
;  Declare variables to read
; ===================================al=======================================

gadget_npart = LonArr(6)
gadget_mass = DblArr(6)
gadget_time = 0.D0
gadget_redshift = 0.D0
gadget_flag_sfr = 0L            ; flags whether the simulation was
                                ; including star formation
gadget_flag_feedback= 0L        ; flags whether feedback was included
gadget_npart_tot = LonArr(6)    ; total number of particles of each type
gadget_flag_cooling = 0L        ; flags whether cooling was included
gadget_num_files = 0L           ; number of files in multi-file snapshot
gadget_boxsize = 0.D0           ; box-size of simulation in case periodic
                                ; boundaries were used
gadget_omega0 = 0.D0            ; matter density in units of critical density
gadget_omegaLambda = 0.D0       ; cosmological constant parameter
gadget_hubble0 = 0.D0           ; Hubble parameter in units of 100 km/sec/Mpc
gadget_flag_stellarage = 0L     ; flags whether the file contains formation
                                ; times of star particles
gadget_flag_metals = 0L         ; flags whether the file contains metallicity
                                ; values for gas and star

; ==========================================================================
;  Read in the header of the first snapshot and get the number of
;  output sub-files
; ==========================================================================

print, 'Reading glass file...'

openr, unit, file, /f77_unformatted, /get_lun

readu, unit, $
  gadget_npart, $
  gadget_mass, $
  gadget_time, $
  gadget_redshift, $
  gadget_flag_sfr, $
  gadget_flag_feedback, $
  gadget_npart_tot, $
  gadget_flag_cooling, $
  gadget_num_files

close, unit
free_lun, unit

print, 'Number of multiple snapshot files', gadget_num_files
print, 'Total number of particles', gadget_npart_tot
print, 'Number of gas particles', gadget_npart_tot[0]
print, 'Mass array', gadget_mass

ngas = gadget_npart_tot[0]
ngas_original = ngas
ndm = gadget_npart_tot[1]
ntot = total(gadget_npart_tot[0])

; ==========================================================================
;  Read in particles
; ==========================================================================

gadget_ntot  = long(total(gadget_npart_tot))
gadget_pos   = FltArr(3,gadget_ntot)
gadget_vel   = FltArr(3,gadget_ntot)
gadget_id    = LonArr(gadget_ntot)
gadget_pmass = FltArr(gadget_ntot)

if gadget_npart_tot[0] gt 0 then begin
    gadget_therm  = FltArr(gadget_npart_tot[0])
    gadget_rho    = FltArr(gadget_npart_tot[0])
    gadget_Ne     = FltArr(gadget_npart_tot[0])
    gadget_NHI    = FltArr(gadget_npart_tot[0])
    gadget_NHeI   = FltArr(gadget_npart_tot[0])
    gadget_NHeIII = FltArr(gadget_npart_tot[0])
    gadget_h      = FltArr(gadget_npart_tot[0])
    gadget_sfr    = FltArr(gadget_npart_tot[0])
endif

itemp1 = LonArr(6)
itemp2 = LonArr(6)

nstart = LonArr(6)
nend   = LonArr(6)

nstart[0] = 0L

for i = 1, 5 do nstart[i] = nstart[i-1] + gadget_npart_tot[i-1]

for ifile = 0, gadget_num_files - 1 do begin

    if gadget_num_files gt 1 then begin
        prefix = strmid(file, 0, strpos(file,'.'))
        file   = prefix + '.' + strtrim(string(ifile), 1)
    endif

    openr, unit, file, /f77_unformatted, /get_lun

    readu, unit, $
      gadget_npart, $
      gadget_mass, $
      gadget_time, $
      gadget_redshift, $
      gadget_flag_sfr, $
      gadget_flag_feedback, $
      gadget_npart_tot, $
      gadget_flag_cooling, $
      gadget_num_files, $
      gadget_boxsize, $
      gadget_omega0, $
      gadget_omegaLambda, $
      gadget_hubble0, $
      gadget_flag_stellarage, $
      gadget_flag_metals

    gadget_subtot = long(total(gadget_npart))

    nend = nstart + gadget_npart - 1L

    itemp1[0] = 0L

    for i = 1, 5 do itemp1[i] = itemp1[i-1] + gadget_npart[i-1]

    itemp2 = itemp1 + gadget_npart - 1L

    temp_f3x1 = FltArr(3,gadget_subtot)
    temp_f1   = FltArr(gadget_subtot)
    temp_i1   = LonArr(gadget_subtot)

    


; ==========================================================================
    print, 'Output ' + strtrim(string(ifile), 1) + ' - read in particle positions...'

    readu, unit, temp_f3x1

    for i = 0, 5 do begin
        if gadget_npart[i] gt 0L then begin
            gadget_pos[*,nstart[i]:nend[i]] = temp_f3x1[*,itemp1[i]:itemp2[i]]
        endif
    endfor

    print, 'Output ' + strtrim(string(ifile), 1) + ' - read in particle positions... DONE'
; ==========================================================================

; ==========================================================================
    print, 'Output ' + strtrim(string(ifile), 1) + ' - read in particle velocities...'

    readu, unit, temp_f3x1

    for i = 0, 5 do begin
        if gadget_npart[i] gt 0L then begin
            gadget_vel[*,nstart[i]:nend[i]] = temp_f3x1[*,itemp1[i]:itemp2[i]]
        endif
    endfor

    print, 'Output ' + strtrim(string(ifile), 1) + ' - read in particle velocities... DONE'
; ==========================================================================

; ==========================================================================
    print, 'Output ' + strtrim(string(ifile), 1) + ' - read in particle IDs...'

    readu, unit, temp_i1

    for i = 0, 5 do begin
        if gadget_npart[i] gt 0L then begin
            gadget_id[nstart[i]:nend[i]] = temp_i1[itemp1[i]:itemp2[i]]
        endif
    endfor

    print, 'Output ' + strtrim(string(ifile), 1) + ' - read in particle IDs... DONE'
; ==========================================================================

; ==========================================================================
    print, 'Output ' + strtrim(string(ifile), 1) + ' - read in particle masses...'

    imass = where(gadget_npart gt 0L and gadget_mass eq 0.0, n_with_mass)

    if n_with_mass gt 0L then begin
        temp_f1 = FltArr(total(gadget_npart[imass]))
        readu, unit, temp_f1
    endif

    itempm1 = 0L

    for i = 0, 5 do begin
        if gadget_npart[i] gt 0L and gadget_mass[i] eq 0.0 then begin
            itempm2 = itempm1 + gadget_npart[i] - 1

            gadget_pmass[nstart[i]:nend[i]] = temp_f1[itempm1:itempm2]

            itempm1 = itempm1 + gadget_npart[i]
        endif else if gadget_npart[i] gt 0L and gadget_mass[i] gt 0.0 then begin
            gadget_pmass[nstart[i]:nend[i]] = gadget_mass[i]
        endif else if gadget_npart[i] eq 0L and gadget_mass[i] eq 0.0 then begin
            print, 'No particles of type ' + strtrim(string(i), 1)
        endif else begin
            print, 'ERROR: Something went terribly wrong'
            stop
        endelse
    endfor

    print, 'Output ' + strtrim(string(ifile), 1) + ' - read in particle masses... DONE'
; ==========================================================================

    print, 'Output ' + strtrim(string(ifile), 1) + ' - read in gas properties...'

    if gadget_npart[0] gt 0 then begin
        temp_f1 = FltArr(gadget_npart[0])

        readu, unit, temp_f1
        gadget_therm[nstart[0]:nend[0]] = temp_f1
        print, "->Reading energy... DONE"
        if not keyword_set(ICs) then begin
            readu, unit, temp_f1
            gadget_rho[nstart[0]:nend[0]] = temp_f1
            print, "->Reading density... DONE"
            if gadget_flag_cooling gt 0 then begin
                readu, unit, temp_f1
                gadget_Ne[nstart[0]:nend[0]] = temp_f1

                readu, unit, temp_f1
                gadget_NHI[nstart[0]:nend[0]] = temp_f1
                print, "->Reading ne, nHI... DONE"
            endif

;            readu, unit, temp_f1
;            gadget_NHeI[nstart[0]:nend[0]] = temp_f1

;            readu, unit, temp_f1
;            gadget_NHeIII[nstart[0]:nend[0]] = temp_f1

            readu, unit, temp_f1
            gadget_h[nstart[0]:nend[0]] = temp_f1
            print, "->Reading kernel... DONE"
            if gadget_flag_sfr gt 0 then begin
                print, "->Reading star formation... DONE"
                readu, unit, temp_f1
                gadget_sfr[nstart[0]:nend[0]] = temp_f1
            endif

        endif

    endif

    print, 'Output ' + strtrim(string(ifile), 1) + ' - read in gas properties... DONE'
; 
    close, unit
    free_lun, unit

    nstart = nstart + gadget_npart

endfor

print, 'Reading glass file...DONE'

if keyword_set(rescaleonly) then begin

; new box size
    new_gadget_boxsize = 0.02 ;; [kpc]


gadget_pos[0, *] *= new_gadget_boxsize / gadget_boxsize
gadget_pos[1, *] *= new_gadget_boxsize / gadget_boxsize
gadget_pos[2, *] *= new_gadget_boxsize / gadget_boxsize

new_pos = gadget_pos[*,0:ngas - 1]

;; set new size
gadget_boxsize = new_gadget_boxsize    

;; discard particles outside box
good = where(gadget_pos[0,*] gt 0 and gadget_pos[1,*] gt 0 and gadget_pos[2,*] gt 0 and gadget_pos[0,*] lt gadget_boxsize  and gadget_pos[1,*] lt gadget_boxsize  and gadget_pos[2,*] lt gadget_boxsize, count)

new_pos = dblarr(count)
gadget_pos = gadget_pos[*,good]
new_pos = gadget_pos
ngas = count
ntot = count

if keyword_set(debug) then begin
    
    xrange = [-gadget_boxsize, 2 * gadget_boxsize] * 1000 
    yrange = [-gadget_boxsize, 2 * gadget_boxsize] * 1000 
    
    x = new_pos[0, *]
    y = new_pos[1, *]
    z = new_pos[2, *]
    
    good2 = where(z gt gadget_boxsize * 0.45  and z lt  gadget_boxsize * 0.55, count)
    
    xplot = x[good2] * 1000  ;;[pc]
    yplot = y[good2] * 1000  ;;[pc]
        
    title = 'Result 1: particle distribution'
    plot, xplot, yplot, $
      psym=3, /iso, xrange=xrange, yrange=yrange, /xs, /ys, xtitle='pc', ytitle='pc', $
      position = position, color = 200, title = title
    
    stop
    
endif
  
    

print, 'Scaling particle positions... DONE'
    

endif else begin

; ==========================================================================
; ; define problem parameters: cored isothermal (1/r^2) density profile
; ==========================================================================

 case 1 of
    
     (setup eq 'minihalo0'): begin

; ;; minihalo
; ;; ==========================================================================
         Xmassfrac = 1 ;; hydrogen mass fraction
         m_H = 1.67262158d-24   ; proton mass [g]

;; core hydrogen number density
         n0 = 1.5e0  ;; [cm^{-3}]
;;core radius
         r0 = 92 / 1000. ;; [kpc] 
         
; new box size
         new_gadget_boxsize = 1.0 ;; [kpc]
         
;; method parameters           ;;ARS can change powerlaw here
         box = gadget_boxsize
         gadget_boxsize = 1.6
         gadget_pos *= gadget_boxsize / box
         expscale = 2.0 ;; scaling to produce 1/r^2 power law
         expfit = -4.8 ;; empirical exponent to smoothly join the densities at r0, depends on gadget_boxsize 
         ;expfit = -8.0
     end


     (setup eq 'minihalo1'): begin

; ;; minihalo
; ;; ==========================================================================
         Xmassfrac = 1 ;; hydrogen mass fraction
         m_H = 1.67262158d-24   ; proton mass [g]

;; core hydrogen number density
         n0 = 3000  ;; [cm^{-3}]
;;core radius
         r0 = 3 / 1000. ;; [kpc] 
         
; new box size
         new_gadget_boxsize = 0.8 ;; [kpc]
         
;; method parameters
         box = gadget_boxsize
         gadget_boxsize = 1.6
         gadget_pos *= gadget_boxsize / box
         expscale = .5 ;; scaling to produce 1/r^1 power law
         ;expfit = -4.8 ;; empirical exponent to smoothly join the densities at r0, depends on gadget_boxsize 
         expfit = -1
     end

 endcase

; ==========================================================================
; ; assign particle positions
; ==========================================================================
; 

; position of source 
StarPos = dblarr(3)
StarPos[0] = gadget_boxsize / 2.0
StarPos[1] = gadget_boxsize / 2.0
StarPos[2] = gadget_boxsize / 2.0


if not keyword_set(MC) then begin

    print, 'Scaling gas particle positions...'

    new_pos = gadget_pos[*, 0:ngas - 1]
    old_pos = gadget_pos[*, 0:ngas - 1]

    r = sqrt((new_pos[0,*] - StarPos[0])^2 + (new_pos[1,0:*] - StarPos[1])^2 + (new_pos[2,0:*] - StarPos[2])^2)
    good = where(r gt 0)

    idx = sort(r)
    rsort = r[idx]
    r0num = rsort[npart0 - 1]
    
    if not keyword_set(nocore) then begin
        
        ;; temporarily remove particles inside r0num; 
        ;; those are inside the core and hence are not stretched
        ;; ==========================================================================

        good = where(r gt r0num, count)
        new_pos = new_pos[*,good]
        
        if keyword_set(debug) then begin
            
            xrange = [0, gadget_boxsize] * 1000 
            yrange = [0, gadget_boxsize] * 1000 
            
            x = new_pos[0, *]
            y = new_pos[1, *]
            z = new_pos[2, *]
            
            good2 = where(z gt gadget_boxsize * 0.45  and z lt  gadget_boxsize * 0.55, count)
            
            xplot = x[good2] * 1000  ;;[pc]
            yplot = y[good2] * 1000  ;;[pc]
            
            title = 'Step 1: cut out core'
            plot, xplot, yplot, $
              psym=3, /iso, xrange=xrange, yrange=yrange, /xs, /ys, xtitle='pc', ytitle='pc', $
              position = position, color = 200, title = title
            
                                ;image = tvread()
                                ;write_jpeg, 'glass-step1.jpg',;image, quality = 75, true = 1
            
            stop
        endif
        
    endif

    ;;scale to a 1/r^2 density profile
    ; ==========================================================================
    scale = (r[good] / r0num)^expscale

    for i = 0, 2 do begin
        
        new_pos[i,*] = StarPos[i] + $ ; stretching center
          ((new_pos[i,*] - StarPos[i]) * scale ); $ ; stretching
        
    endfor    
    
    if keyword_set(debug) then begin

        xrange = [0, gadget_boxsize] * 1000 
        yrange = [0, gadget_boxsize] * 1000 
        
        x = new_pos[0, *]
        y = new_pos[1, *]
        z = new_pos[2, *]
        
        good2 = where(z gt gadget_boxsize * 0.45  and z lt  gadget_boxsize * 0.55, count)
    
        xplot = x[good2] * 1000  ;;[pc]
        yplot = y[good2] * 1000  ;;[pc]
        
        title = 'Step 2: scale to 1/r^2 profile'
        plot, xplot, yplot, $
          psym=3, /iso, xrange=xrange, yrange=yrange, /xs, /ys, xtitle='pc', ytitle='pc', $
          position = position, color = 300, title = title

        ;image = tvread()
        ;write_jpeg, 'glass-step2.jpg',;image, quality = 75, true = 1

        stop
        
    endif
   
    ;;scale to match core density
    ; ==========================================================================
  ;  if not keyword_set(nocore) then begin   

        for i = 0, 2 do begin
            
            new_pos[i,*] = StarPos[i] + $ ; center
              ((new_pos[i,*] - StarPos[i]) /(gadget_boxsize / 2.)^expfit) ; re-scaling 
        
        endfor    
    
        if keyword_set(debug) then begin
            
            xrange = [0, gadget_boxsize] * 1000 
            yrange = [0, gadget_boxsize] * 1000 
        
            x = new_pos[0, *]
            y = new_pos[1, *]
            z = new_pos[2, *]
            
            good2 = where(z gt gadget_boxsize * 0.45  and z lt  gadget_boxsize * 0.55, count)
    
            xplot = x[good2] * 1000  ;;[pc]
            yplot = y[good2] * 1000  ;;[pc]
        
            title = 'Step 3: rescale by spatially constant factor to match core density'
            plot, xplot, yplot, $
              psym=3, /iso, xrange=xrange, yrange=yrange, /xs, /ys, xtitle='pc', ytitle='pc', $
              position = position, color = 500, title = title
            
                                ;image = tvread()
                                ;write_jpeg, 'glass-step3.jpg',;image, quality = 75, true = 1
            
            stop
            
        endif
        
   ; endif

    if not keyword_set(nocore) then begin   
        
        ;; keep only particles outside r0num; the volume inside r0num
        ;; will be filled below with a uniform glass-like core 
        ;; ==========================================================================
        
        rscaled = sqrt((new_pos[0,*] - StarPos[0])^2 + (new_pos[1,*] - StarPos[1])^2 + (new_pos[2,*] - StarPos[2])^2)
        good = where(rscaled gt r0num, count)
        new_pos = new_pos[*, good]

        
        if keyword_set(debug) then begin
            
            xrange = [0, gadget_boxsize] * 1000 
            yrange = [0, gadget_boxsize] * 1000 
            
            x = new_pos[0, *]
            y = new_pos[1, *]
            z = new_pos[2, *]
            
            good2 = where(z gt gadget_boxsize * 0.45  and z lt  gadget_boxsize * 0.55, count)
            
            xplot = x[good2] * 1000  ;;[pc]
            yplot = y[good2] * 1000  ;;[pc]
            
            title = 'Step 4: clear core'
            plot, xplot, yplot, $
              psym=3, /iso, xrange=xrange, yrange=yrange, /xs, /ys, xtitle='pc', ytitle='pc', $
              position = position, color = 500, title = title
            
                                ;image = tvread()
                                ;write_jpeg, 'glass-step4.jpg',;image, quality = 75, true = 1
            stop
            
        endif

    endif

    ;;add particles inside r0num from original glass file to create the
    ;;uniform core 
    ;;ARS - or a powerlaw core?
    ;; ==========================================================================
        
    if not keyword_set(nocore) then begin   
        
        good = where(r lt r0num, count)

        pos = dblarr(3, n_elements(new_pos[0,*]) + count)
        
        pos[*, 0:n_elements(new_pos[0,*]) - 1] = new_pos
        
        if not keyword_set(coremc) then begin
            
            pos[*, n_elements(new_pos[0,*]):n_elements(new_pos[0,*]) + count - 1] = old_pos[*, good]
            
        endif else begin
            
            ;;Monte Carlo sampling 
            
            for i =  n_elements(new_pos[0,*]), n_elements(new_pos[0,*]) + count - 1 do begin
            
                drandom = randomu(seed, 1, /double) ;
                r = r0num * (drandom)^(1/3.) ;
                
                drandom = randomu(seed, 1, /double) ;
                theta = acos(1 - 2 * drandom) ;
                
                drandom = randomu(seed, 1, /double) ;
                phi = 2 * !PI * drandom ;
                
                ;;transform to cartesian coordiantes
                pos[0,i] = StarPos[0] + r * sin(theta) * cos(phi) ;
                pos[1,i] = StarPos[1] + r * sin(theta) * sin(phi) ;
                pos[2,i] = StarPos[2] + r * cos(theta) ;
            
            endfor    
        
        endelse


        new_pos = pos
        r = sqrt((new_pos[0,*] - StarPos[0])^2 + (new_pos[1,*] - StarPos[1])^2 + (new_pos[2,*] - StarPos[2])^2)
        
        if keyword_set(debug) then begin
            
            xrange = [0, gadget_boxsize] * 1000 
            yrange = [0, gadget_boxsize] * 1000 
            
            x = new_pos[0, *]
            y = new_pos[1, *]
            z = new_pos[2, *]
            
            good2 = where(z gt gadget_boxsize * 0.45  and z lt  gadget_boxsize * 0.55, count)
            
            xplot = x[good2] * 1000  ;;[pc]
            yplot = y[good2] * 1000  ;;[pc]
            
            title = 'Step 5: fill core with glass'
            plot, xplot, yplot, $
              psym=3, /iso, xrange=xrange, yrange=yrange, /xs, /ys, xtitle='pc', ytitle='pc', $
              position = position, color = 500, title = title
            
                                ;image = tvread()
                                ;write_jpeg, 'glass-step5.jpg',;image, quality = 75, true = 1
            
            stop
            
        endif
        

        ;;new particle coordinates
        ;; ==========================================================================
        ngas = n_elements(new_pos[0,*])
        ntot = ngas             ; + ndm
        
        pos = fltarr(3, ntot) ;; note: fltarr
        pos[*, 0:ngas - 1] = new_pos
                                ;pos[*, ngas:ntot - 1] = gadget_pos[*, ngas:ntot - 1]
        
        gadget_pos = pos

    endif else begin
        
        ngas = n_elements(new_pos[0,*])
        ntot = ngas             ; + ndm
        gadget_pos = float(new_pos) 
        
    endelse

        
    ;; rescale problem from core radius r0num to core radius r0
    ;; ==========================================================================
    r = sqrt((gadget_pos[0,0:ngas - 1] - StarPos[0])^2 + (gadget_pos[1,0:ngas-1] - StarPos[1])^2 + (gadget_pos[2,0:ngas-1] - StarPos[2])^2)
        
    
    for i = 0, 2 do begin
        
        gadget_pos[i,*] = StarPos[i] + $ ;  center
          (gadget_pos[i,*] - StarPos[i]) * r0 / r0num ; rescaling
        
    endfor    
    
    new_pos = gadget_pos[*, 0:ngas-1]
    gadget_pos = new_pos 
   
    if keyword_set(debug) then begin
        
        xrange = [0, gadget_boxsize] * 1000 
        yrange = [0, gadget_boxsize] * 1000 
        
        x = new_pos[0, *]
        y = new_pos[1, *]
        z = new_pos[2, *]
        
        good2 = where(z gt gadget_boxsize * 0.45  and z lt  gadget_boxsize * 0.55, count)
        
        xplot = x[good2] * 1000  ;;[pc]
        yplot = y[good2] * 1000  ;;[pc]
        
        title = 'Step 6: rescale to desired core radius'
        plot, xplot, yplot, $
          psym=3, /iso, xrange=xrange, yrange=yrange, /xs, /ys, xtitle='pc', ytitle='pc', $
          position = position, color = 500, title = title
        
                                ;image = tvread()
                                ;write_jpeg, 'glass-step6.jpg',;image, quality = 75, true = 1
        
        stop
        
    endif
        
   
endif else begin

    
    print, 'Monte Carlo sampling particle positions...'
    seed = 101L
    

    ;; particles inside r0 
    new_pos = dblarr(3, npart0) 
    
    for i = 0L, npart0 - 1 do begin
        
        drandom = randomu(seed, 1, /double) ;
        r = r0 * (drandom)^(1/3.) ;
      
        drandom = randomu(seed, 1, /double) ;
        theta = acos(1 - 2 * drandom) ;

        drandom = randomu(seed, 1, /double) ;
        phi = 2 * !PI * drandom ;
    
        ;;transform to cartesian coordiantes
        new_pos[0,i] = StarPos[0] + r * sin(theta) * cos(phi) ;
        new_pos[1,i] = StarPos[1] + r * sin(theta) * sin(phi) ;
        new_pos[2,i] = StarPos[2] + r * cos(theta) ;

    endfor    
    
    if keyword_set(debug) then begin
        
        xrange = [0, gadget_boxsize] * 1000 
        yrange = [0, gadget_boxsize] * 1000 
        
        x = new_pos[0, *]
        y = new_pos[1, *]
        z = new_pos[2, *]
        
        good2 = where(z gt gadget_boxsize * 0.45  and z lt  gadget_boxsize * 0.55, count)
    
        xplot = x[good2] * 1000  ;;[pc]
        yplot = y[good2] * 1000  ;;[pc]
        
        title = 'Step 1: sample core'
        plot, xplot, yplot, $
          psym=3, /iso, xrange=xrange, yrange=yrange, /xs, /ys, xtitle='pc', ytitle='pc', $
          position = position, color = 500, title = title

        ;image = tvread()
        ;write_jpeg, 'mc-step1.jpg',;image, quality = 75, true = 1

        stop
        
    endif


    gadget_pos = new_pos
    ngas = npart0

    ;; particles outside r0 and inside r < sqrt(3) * boxhalf
    ;; (spherical instead of boxy region for
    ;; easier integration of total mass)
    
    Rbox = sqrt(3) * gadget_boxsize / 2.0
    
    
    nleft = (Rbox - r0) * 3. / r0 * double(npart0)
    print, nleft, ngas, nleft+ ngas
  
    if nleft gt 0 then begin
        
        new_pos = dblarr(3, nleft)
        
        ndone = 0L
        for i = 0L, nleft - 1 do begin

            drandom = randomu(seed, 1, /double) ;
            r = Rbox * drandom  ; // particles cover the whole box
            
            if(r lt r0) then begin
                i--
                continue
            endif
            
            drandom = randomu(seed, 1, /double) ;
            theta = acos(1 - 2 * drandom) ;
            
            drandom = randomu(seed, 1, /double) ;
            phi = 2 * !PI * drandom ;
            
            ;;transform to cartesian coordiantes
            new_pos[0,i] = StarPos[0] + r * sin(theta) * cos(phi) ;
            new_pos[1,i] = StarPos[1] + r * sin(theta) * sin(phi) ;
            new_pos[2,i] = StarPos[2] + r * cos(theta) ;
            
            ndone++
        endfor    
        print, ndone
    

        if keyword_set(debug) then begin
            
            xrange = [0, gadget_boxsize] * 1000 
            yrange = [0, gadget_boxsize] * 1000 
            
            x = new_pos[0, *]
            y = new_pos[1, *]
            z = new_pos[2, *]
            
            good2 = where(z gt gadget_boxsize * 0.45  and z lt  gadget_boxsize * 0.55, count)
            
            xplot = x[good2] * 1000  ;;[pc]
            yplot = y[good2] * 1000  ;;[pc]
            
            title = 'Step 2: sample mass outside core'
            plot, xplot, yplot, $
              psym=3, /iso, xrange=xrange, yrange=yrange, /xs, /ys, xtitle='pc', ytitle='pc', $
              position = position, color = 500, title = title
            
            ;image = tvread()
            ;write_jpeg, 'mc-step2.jpg',;image, quality = 75, true = 1
            
            stop
            
        endif

    endif else begin
        nleft = 0
    endelse
    
    ;;new particle coordinates
    ; ==========================================================================
    ngas = n_elements(new_pos[0,*])
    ntot = npart0+nleft

    pos = fltarr(3, ntot) ;; note: fltarr
    pos[*, 0:ngas - 1] = new_pos

    if ntot gt ngas then begin
        pos[*, ngas:ntot - 1] = gadget_pos[*, *]
    endif

    gadget_pos = pos
    ngas = npart0 + nleft

    new_pos = gadget_pos

    if keyword_set(debug) then begin
        
        xrange = [0, gadget_boxsize] * 1000 
        yrange = [0, gadget_boxsize] * 1000 
        
        x = new_pos[0, *]
        y = new_pos[1, *]
        z = new_pos[2, *]
        
        good2 = where(z gt gadget_boxsize * 0.45  and z lt  gadget_boxsize * 0.55, count)
    
        xplot = x[good2] * 1000  ;;[pc]
        yplot = y[good2] * 1000  ;;[pc]
        
        title = 'Step 3: combine'
        plot, xplot, yplot, $
          psym=3, /iso, xrange=xrange, yrange=yrange, /xs, /ys, xtitle='pc', ytitle='pc', $
          position = position, color = 500, title = title

        ;image = tvread()
        ;write_jpeg, 'mc-step3.jpg',;image, quality = 75, true = 1
        
        stop
        
    endif

    print, 'Monte Carlo sampling particle positions...DONE'

endelse


; ==========================================================================
;; define new box
; ==========================================================================

;;shift center

cxold = gadget_boxsize / 2.0
cyold = gadget_boxsize / 2.0
czold = gadget_boxsize / 2.0

cxnew = new_gadget_boxsize / 2.0
cynew = new_gadget_boxsize / 2.0
cznew = new_gadget_boxsize / 2.0

dx = cxnew-cxold
dy = cynew-cyold
dz = cznew-czold

gadget_pos[0, *] += dx
gadget_pos[1, *] += dy
gadget_pos[2, *] += dz

StarPos[0] += dx
StarPos[1] += dy
StarPos[2] += dz

new_pos = gadget_pos[*,0:ngas - 1]

;; set new size
gadget_boxsize = new_gadget_boxsize    

;; discard particles outside box
good = where(gadget_pos[0,*] gt 0 and gadget_pos[1,*] gt 0 and gadget_pos[2,*] gt 0 and gadget_pos[0,*] lt gadget_boxsize  and gadget_pos[1,*] lt gadget_boxsize  and gadget_pos[2,*] lt gadget_boxsize, count)

new_pos = dblarr(count)
gadget_pos = gadget_pos[*,good]
new_pos = gadget_pos
ngas = count
ntot = count


if keyword_set(rin) then begin
;;discard particles inside rin
    r = sqrt((new_pos[0,*] - StarPos[0])^2 + (new_pos[1,0:*] - StarPos[1])^2 + (new_pos[2,0:*] - StarPos[2])^2)
    good = where(r gt rin, count)
    gadget_pos = gadget_pos[*,good]
    ngas -= (ngas - count)
    ntot -= (ngas - count)
   
endif


if keyword_set(debug) then begin
    
    xrange = [0, gadget_boxsize] * 1000 
    yrange = [0, gadget_boxsize] * 1000 
    
    x = new_pos[0, *]
    y = new_pos[1, *]
    z = new_pos[2, *]
    
    good2 = where(z gt gadget_boxsize * 0.45  and z lt  gadget_boxsize * 0.55, count)
    
    xplot = x[good2] * 1000  ;;[pc]
    yplot = y[good2] * 1000  ;;[pc]
        
    title = 'Result 1: particle distribution'
    plot, xplot, yplot, $
      psym=3, /iso, xrange=xrange, yrange=yrange, /xs, /ys, xtitle='pc', ytitle='pc', $
      position = position, color = 500, title = title
    
    ;image = tvread()
    ;write_jpeg, 'result1.jpg',;image, quality = 75, true = 1

    stop
    
endif
  



    
print, 'Old number of gas particles inside box: ', gadget_npart_tot[0]

print, 'New number of gas particles inside box: ', ngas
print, 'New number of  particles inside box: ', ntot

r = sqrt((gadget_pos[0,0:ngas - 1] - StarPos[0])^2 + (gadget_pos[1,0:ngas-1] - StarPos[1])^2 + (gadget_pos[2,0:ngas-1] - StarPos[2])^2)
 
good = where(r lt r0, count)
print, 'New number of gas particles inside r0: ', count

good = where(r gt r0, count)
print, 'New number of gas particles outside r0: ', count

print, 'Scaling particle positions... DONE'

; ==========================================================================
;; cut particles in a ring next to the core 
; ==========================================================================

; r = sqrt((gadget_pos[0,*] - StarPos[0])^2 + (gadget_pos[1,*] - StarPos[1])^2 + (gadget_pos[2,*] - StarPos[2])^2)
; r1 = r0
; r2 = r0 + r0/5.
; good = where(r lt r1 or r gt r2, count)

; gadget_pos = gadget_pos[*,good]
; new_pos = gadget_pos
; ngas = count
; ntot = count

; if keyword_set(debug) then begin
    
;     xrange = [0, gadget_boxsize] * 1000 
;     yrange = [0, gadget_boxsize] * 1000 
    
;     x = new_pos[0, *]
;     y = new_pos[1, *]
;     z = new_pos[2, *]
    
;     good2 = where(z gt gadget_boxsize * 0.45  and z lt  gadget_boxsize * 0.55, count)
    
;     xplot = x[good2] * 1000  ;;[pc]
;     yplot = y[good2] * 1000  ;;[pc]
        
;     title = 'After cutting particle ring'
;     plot, xplot, yplot, $
;       psym=3, /iso, xrange=xrange, yrange=yrange, /xs, /ys, xtitle='pc', ytitle='pc', $
;       position = position, color = fsc_color('black'), title = title
    
;     ;image = tvread()
;     ;write_jpeg, 'cutring.jpg',;image, quality = 75, true = 1

;     stop
    
; endif


; ==========================================================================
; ; assign particle masses
; ==========================================================================

r = sqrt((gadget_pos[0,0:ngas - 1] - StarPos[0])^2 + (gadget_pos[1,0:ngas-1] - StarPos[1])^2 + (gadget_pos[2,0:ngas-1] - StarPos[2])^2)

if not keyword_set(nocore) then begin
    
    good = where(r le r0, count) 
    
    nr0 = count ;; number of particles inside core
    rho_0 = n0 / Xmassfrac * m_H ;; core density
    kpc = 3.0857d21
    mr0 = 4/3. * !Pi * (r0 * kpc)^3 * rho_0 ; core mass
    print, 'total mass: ', mr0 /  1.989d43 
    m_p = mr0 / double(nr0) ;; gas particle mass
    
    if keyword_set(MC) then m_p /= 1.35 ;; to account for Monte Carlo bias
    
endif else begin
    
    r_halo = gadget_boxsize / 10.
    good = where(r le r_halo, count) 
    
    rho_0 = n0 / Xmassfrac * m_H ;; core density
    
    kpc = 3.0857d21
    
    if keyword_set(rin) then begin
        mr0 = 4 * !Pi * rho_0 * (r0 * kpc)^2 * (r_halo - rin) * kpc
    endif else begin
        mr0 = 4 * !Pi * rho_0 * (r0 * kpc)^2 * r_halo * kpc
    endelse

    print, 'total mass: ', mr0 /  1.989d43 
         
    m_p = mr0 / double(count) ;; gas particle mass

endelse

gadget_mass[0] = m_p 
print, 'particle mass [code units]: ', m_p /  1.989d43 
print, 'particle mass [solar mass]: ', m_p /  1.989d33 


; ==========================================================================
; ; compute density profile
; ==========================================================================
; 
;; compute gas density in log bins
lrmax = alog10(max(r))
lrmin = alog10(min(r))
nbin = 100
    
rc = dblarr(nbin-1)
vol = dblarr(nbin-1)
rhoav = dblarr(nbin-1)

dlr = (lrmax - lrmin) / double(nbin)
    

for i = 0l, nbin - 2l do begin

    rmin = 10^(lrmin + dlr * i)
    rmax = 10^(lrmin + dlr * (i + 1))
      
    rc[i] = (rmin + rmax) / 2.
        
    idx = where(r ge rmin and r le rmax, count)
  
    if count gt 1 then begin
        vol[i] = 4 * !pi / 3. * ((rmax * kpc)^3 - (rmin * kpc)^3)
        rhoav[i] = count * gadget_mass[0] / vol[i] ;; equal mass particles
    endif else begin
        ;;   rc[i] = !VALUES.D_NAN
        rhoav[i] = !VALUES.D_NAN
    endelse

endfor


; ==========================================================================
; ; plot profile
; ==========================================================================

;; choose units convenient for plotting
rhoav = rhoav * Xmassfrac / m_H ;; hydrogen number density [cm^{-3}]
rc *= 1000; [pc]
r0 *= 1000; [pc]
boxhalf = gadget_boxsize/2. * 1000 ;[pc]

title = 'Result 2: density profile'
plot, alog10(rc), alog10(rhoav), psym = 1, lines = 0, color = 500, $
  xtitle = 'log!D10!N r [pc]', ytitle = 'log!D10!N n_H [cm^{-3}]', $ 
  xrange = alog10([min(rc),max(rc)]), yrange = alog10([min(rhoav, /nan),max(rhoav, /nan)]), $
  /nodata, title = title

;; input profile
if not keyword_set(nocore) then begin
    
    good = where(rc gt r0, count)
    if count gt 0 then $
      oplot, alog10(rc[good]), alog10(n0 / (rc[good] / r0)^3), psym = 0, lines = 0, color = 500
    
    good = where(rc le r0)
    if count gt 0 then $
      oplot, alog10(rc[good]), alog10(n0 * rc[good]/rc[good]), psym = 0, lines = 0, color = 500

endif else begin
    
    good = where(rc gt 0, count)
    if count gt 0 then $
      oplot, alog10(rc[good]), alog10(n0 / (rc[good] / r0)^3), psym = 0, lines = 0, color = 500
 
endelse

;; box size
oplot, alog10([boxhalf, boxhalf]), [-10,10], lines = 2, color = 500

;;; core radius
oplot, alog10([r0, r0]), [-10,10], lines = 2, color = 500

;; output profile
oplot, alog10(rc), alog10(rhoav), psym = 10, lines = 0, color = 100

if keyword_set(debug) then begin

    ;image = tvread()
    ;write_jpeg, 'result2.jpg',;image, quality = 75, true = 1

endif

stop

; ; ; ; ; plot particle position in a slice
; ; 

if keyword_set(debug) then begin
    stop
    x = gadget_pos[0, 0:ngas - 1]
    y = gadget_pos[1, 0:ngas - 1]
    z = gadget_pos[2, 0:ngas - 1]
    
    xplot = x * 1000
    yplot = y * 1000

    xrange = [-gadget_boxsize, 2*gadget_boxsize] *1000
    yrange = [-gadget_boxsize, 2*gadget_boxsize] * 1000
    
    plot, xplot, yplot, $
      psym=3, /iso, xrange=xrange, yrange=yrange, /xs, /ys, xtitle='pc', ytitle='pc', $
      position = position, color = 500

 stop
endif

endelse

; ==========================================================================
; save modified ICs
; ==========================================================================

print, 'Writing initial conditions to file...'

openw, unit, /get_lun, outfile, /f77_unformatted

massarr = dblarr(6)
flag_sfr = 0L
flag_feedback = 0L
bytesleft = 256-6*4 - 6*8 - 8 - 8 - 2*4
la = intarr(bytesleft/2)
;bytesleft = 256 - 6*4 - 6*8 - 2*8 - 2*4 - 6*4 - 2*4 - 4*8
;la = intarr(bytesleft)

;;gas only
massarr[*] = 0 
massarr[0] = gadget_mass[0] /  1.989d43 ;; code units 

gadget_npart_tot[*] = 0 ;
gadget_npart_tot[0] = ngas ;

writeu, unit, gadget_npart_tot, massarr, gadget_time, gadget_redshift, flag_sfr, flag_feedback, la

writeu, unit, gadget_pos[*,0:ngas-1]
writeu, unit, gadget_pos[*,0:ngas-1] ;; dummy velocity

writeu, unit, lindgen(ngas)+1L ;; ids
writeu, unit, lindgen(ngas)+1L ;; dummy temperature

close,unit
free_lun, unit

print, 'Writing initial conditions to file... DONE'

end


pro makeplot_nfw

;npart0 = 35
;npart0 = 310
;npart0 = 5000
;npart0 = 10000
;npart0= 100000
;create_ic_nfw, "glass128", "minihalo_128", npart0, setup = 'minihalo0', /ICs, /debug;, /nocore;, /rescaleonly

stop

;stop
;npart0 = 800000
;npart0 = 5000
;create_ic_nfw, "glass128", "minihalo_128_nfw2", npart0, setup = 'minihalo1', /ICs, /debug, /nocore

;setup for the dark stars!!
;stop
;npart0 = 800000
;create_ic_nfw, "glass256", "minihalo_256_nfw", npart0, setup = 'minihalo1', /ICs, /debug, /nocore

;setup for the stampede test!
;stop
npart0 = 131000
create_ic_nfw, "glass064", "minihalo_064_stampede", npart0, setup = 'minihalo0', /ICs, /debug;, /nocore

;stop
;npart0 = 14
;create_ic, "glass016", "minihalo-r01e-3n03e3-Ms025-glass", npart0, setup = 'minihalo0', /ICs;, /debug


;stop
;npart0 = 7
;create_ic, "glass016", "minihalo-r01e-3n03e3-Ms050-glass", npart0, setup = 'minihalo0', /ICs;, /debug
;npart0 = 3
;create_ic, "glass016", "minihalo-r01e-3n03e3-Ms100-glass", npart0, setup = 'minihalo0', /ICs;, /debug

; stop





;create_ic, "glass016", "yoshida-N064", 360., setup = 'yoshida', /ICs, /nocore, rin = 1e-4;, /coremc

;stop
;create_ic, "glass016", "yoshida-N032", 45., setup = 'yoshida', /ICs, /nocore, rin = 1e-4, /debug;, /coremc

;stop

;create_ic, "glass016", "yoshida-N048", 151., setup = 'yoshida', /ICs, /nocore, rin = 1e-4 ;, /coremc

;stop



;stop


;stop
;stop


;create_ic, "glass016", "minihalo-r03e-3n03e3-Ms025-glass-nocore", 1000., setup = 'minihalo1', /ICs, /nocore, /debug, rin = 1e-4;, /coremc



;stop
;stop

;create_ic, "glass016", "minihalo-r03e-3n03e3-Ms050-glass-nocore", 500., setup = 'minihalo1', /ICs, /nocore;, /debug;, /coremc
;create_ic, "glass016", "minihalo-r03e-3n03e3-Ms100-glass-nocore", 250., setup = 'minihalo1', /ICs, /nocore;, /debug;, /coremc

;stop

;npart0 = 50
;create_ic, "glass016", "minihalo-r03e-3ne4-Ms025-glass", npart0, setup = 'minihalo4', /ICs ;, /debug

;stop


; stop

; create_ic, "icsL013.2N032-Glass032-Stroemgren", "stroemgren-L1e-2N032-glass", npart0, setup = 'minihalo1', /ICs, /rescaleonly;, /debug
; create_ic, "icsL013.2N032-Glass032-Stroemgren", "stroemgren-L5e-3N032-glass", npart0, setup = 'minihalo1', /ICs, /rescaleonly;, /debug
; create_ic, "icsL013.2N064-Glass032-Stroemgren", "stroemgren-L5e-3N064-glass", npart0, setup = 'minihalo1', /ICs, /rescaleonly;, /debug
; create_ic, "icsL013.2N128-Glass032-Stroemgren", "stroemgren-L5e-3N128-glass", npart0, setup = 'minihalo1', /ICs, /rescaleonly;, /debug

; stop



end

pro modeinfo_uves,mode,newhead,object,orient,xlo,xhi,ylo,yhi,gain,readn,backg,time $
                       ,obslat,obslon,obsalt,ra2000,de2000,jd,NONLINEAR=nonlinear
;===========================================================================================
;The ESO VLT UVES spectrometer (red and blue arms)
    if mode eq 'uves_blue' then begin           ; EEV blue
      id1 = '1'                     ; gain, readn and sizes
      id2 = '5'                     ; drk and sky
      reorient = 6
    endif
    if mode eq 'uves_middle' then begin         ; EEV red
      id1 = '1' ;for the new format (2 CCDs in 2 separate FITS) 
;      id1 = '4' ;for the old format (2 CCDs in one FITS)
      id2 = '6'
      reorient = 1
    endif
    if mode eq 'uves_red' then begin            ; MIT red
      id1 = '1'
      id2 = '6'
      reorient = 1
    endif
    if n_elements(orient) eq 0 then orient = reorient
    if not keyword_set(xr) then begin
      presc = hierarch(newhead, 'HIERARCH ESO DET OUT'+id1+' PRSCX', count=count, errmes=errmes) ; prescan along X axis
      if count eq 0 then presc = 0
      ovrsc = hierarch(newhead, 'HIERARCH ESO DET OUT'+id1+' OVSCX', count=count, errmes=errmes) ; oversan along X axis
      if count eq 0 then ovrsc = 0
      nxa = hierarch(newhead, 'NAXIS1', count=count, errmes=errmes)                              ; full size of the whole frame
      nxd = hierarch(newhead, 'HIERARCH ESO DET OUT'+id1+' NX', count=count, errmes=errmes)      ; valid pixel range for a given CCD
;      if mode eq 'uves_middle' then nxn = nxa - presc - nxd - ovrsc else nxn = 0                 ; offset for 'uves_middle'
;      xlo = nxn + presc
;      xhi = nxn + presc + nxd - 1
      xlo = presc
      xhi = nxa - ovrsc -1
    endif else begin
      xlo = xr(0)
      xhi = xr(1)
    endelse
    if not keyword_set(yr) then begin
      ny = hierarch(newhead, 'NAXIS2', count=count, errmes=errmes)
      ylo = 0
      yhi = ny - 1; - 101
    endif else begin
      ylo = yr(0)
      yhi = yr(1)
    endelse
    if n_elements(gain) eq 0 then begin
      gain = hierarch(newhead, 'HIERARCH ESO DET OUT'+id1+' CONAD', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
    if n_elements(readn) eq 0 then begin
      readn = hierarch(newhead, 'HIERARCH ESO DET OUT'+id1+' RON', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
    if n_elements(backg) eq 0 then begin
      drk = hierarch(newhead, 'HIERARCH ESO INS DET'+id2+' OFFDRK', count=count, errmes=errmes)
      if count eq 0 then message, errmes, /CONTINUE
      sky = hierarch(newhead, 'HIERARCH ESO INS DET'+id2+' OFFSKY', count=count, errmes=errmes)
      if count eq 0 then message, errmes, /CONTINUE
      backg = gain * (drk + sky)        ;convert ADU to electrons
    endif
    if n_elements(time) eq 0 then begin
      time = hierarch(newhead, 'EXPTIME', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
; For object frames, prepare for heliocentric correction
    imtype=hierarch(newhead,'OBJECT',count=count,errmes=errmes)
    obsctg=hierarch(newhead,'HIERARCH ESO DPR CATG',count=count,errmes=errmes)
    if(strmatch(imtype,'object*',/fold_case) or $
       strmatch(obsctg,'science*',/fold_case)) then begin
      if n_elements(ra2000) eq 0 then begin
        ra2000 = hierarch(newhead,'RA',count=count, errmes=errmes)
	ra2000 = ra2000/15.
      endif
      if n_elements(de2000) eq 0 then begin
        de2000 = hierarch(newhead,'DEC',count=count, errmes=errmes)
      endif
      if n_elements(jd) eq 0 then begin
        jd = hierarch(newhead,'MJD-OBS',count=count,errmes=errmes)+time/2./3600./24.+0.5d0
      endif
;observatory coordinates
      if n_elements(obslon) eq 0 then obslon = -hierarch(newhead,'HIERARCH ESO TEL GEOLON',count=count, errmes=errmes)
      if n_elements(obslat) eq 0 then obslat = hierarch(newhead,'HIERARCH ESO TEL GEOLAT',count=count, errmes=errmes)
      if n_elements(obsalt) eq 0 then obsalt = hierarch(newhead,'HIERARCH ESO TEL GEOELEV',count=count, errmes=errmes)
    endif

return
end

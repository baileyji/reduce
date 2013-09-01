pro modeinfo_bhros,mode,newhead,object,orient,xlo,xhi,ylo,yhi,gain,readn,backg,time $
                       ,obslat,obslon,obsalt,ra2000,de2000,jd,NONLINEAR=nonlinear
;===========================================================================================
;Gemini-South bHROS spectrometer (single CCD in the extention)
;Tested for 1x1 and 1x4 binnings
    if n_elements(orient) eq 0 then orient = 1
    datasec=hierarch(newhead, 'DATASEC', count=count, errmes=errmes) ; DATASEC = '[1:2048,1:4608]' format
    datasec=strmid(datasec,1,strlen(datasec)-2)
    datasec=byte(datasec)
    ii=where(datasec eq (byte(':'))[0], nii)
    if(nii gt 0) then datasec(ii)=byte(',')
    datasec=string(datasec)
    xr_ccd=intarr(2)
    yr_ccd=intarr(2)
    reads,datasec,xr_ccd,yr_ccd
    if not keyword_set(xr) then begin
      xlo = xr_ccd[0]-1
      xhi = xr_ccd[1]-1
    endif else begin
      xlo = xr[0]
      xhi = xr[1]
    endelse
    if not keyword_set(yr) then begin
      ylo = yr_ccd[0]-1
      yhi = yr_ccd[1]-1
    endif else begin
      ylo = yr[0]
      yhi = yr[1]
    endelse

; For object frames, read CCD settings and prepare for heliocentric correction
    imtype=hierarch(newhead,'OBSTYPE',count=count,errmes=errmes)
    if strmatch(imtype,'object*',/fold_case) then begin
      if n_elements(gain) eq 0 then begin
        gain = hierarch(newhead, 'GAIN', count=count, errmes=errmes)
        if count eq 0 then message, errmes
      endif else gain=1.0
      if n_elements(readn) eq 0 then begin
        readn = hierarch(newhead, 'RDNOISE', count=count, errmes=errmes)
      endif else readn=0.0
      backg = gain * 0.1        ;convert ADU to electrons
      if n_elements(time) eq 0 then begin
        time = hierarch(newhead, 'EXPTIME', count=count, errmes=errmes)
        if count eq 0 then message, errmes
      endif else time=0.0
      if n_elements(ra2000) eq 0 then begin
        ra = hierarch(newhead,'RA',count=count, errmes=errmes) ; RA in degrees
      endif
      if n_elements(de2000) eq 0 then begin
        dec = hierarch(newhead,'DEC',count=count, errmes=errmes)
      endif
      if n_elements(epoch) eq 0 then begin
        epoch=hierarch(newhead,'EPOCH',count=count,errmes=errmes)
      endif
      precess, ra, dec, epoch, 2000.
      ra2000=ra/15.d0
      de2000=dec
      if n_elements(jd) eq 0 then begin
        date_obs=hierarch(newhead,'DATE-OBS',count=count,errmes=errmes)
        date_obs=strsplit(date_obs,'-',/EXTRACT)
        ut=hierarch(newhead,'UTSTART',count=count,errmes=errmes)
        ut=float(strsplit(ut,':',/EXTRACT))
        ut=ut(0)+ut(1)/60.d0+(ut(2)+time*0.5d0)/3600.d0
        jdcnv,date_obs[0],date_obs[1],date_obs[2],ut,jd  ;covert UT to JD
        jd=jd-2400000.d0
      endif
;observatory coordinates
      if n_elements(obslon) eq 0 then obslon = ten( 70,44,12.096) ;observatory coordinates, seconds are unknown
      if n_elements(obslat) eq 0 then obslat = ten(-30,14,26.700)
      if n_elements(obsalt) eq 0 then obsalt = 2722.
    endif else begin
      gain=1.
      readn=1.
      time=1.
      backg=0.      
    endelse

return
end

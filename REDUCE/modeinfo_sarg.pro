pro modeinfo_sarg,mode,newhead,object,orient,xlo,xhi,ylo,yhi,gain,readn,backg,time $
                       ,obslat,obslon,obsalt,ra2000,de2000,jd,NONLINEAR=nonlinear
;===========================================================================================
;The TNG at La Palma, SARG spectrometer, RED and BLUE modes, 2 CCDs, 1x1 binning
    if n_elements(orient) eq 0 then orient = 6
    if not keyword_set(xr) then begin
     if mode eq 'sarg_red1' or mode eq 'sarg_blue1' then xr = [50, 2097]
     if mode eq 'sarg_red2'  then xr = [2198, 3730]
     if mode eq 'sarg_blue2' then xr = [2198, 4245]
    endif
    xlo = xr(0)
    xhi = xr(1)

    if not keyword_set(yr) then begin
     if mode eq 'sarg_red1' then yr = [0, 3970]
     if mode eq 'sarg_red2' then yr = [0, 3805]
     if mode eq 'sarg_blue1' then yr = [0, 3790]
     if mode eq 'sarg_blue2' then yr = [0, 3965]
    endif
    ylo = yr(0)
    yhi = yr(1)

    if n_elements(gain) eq 0 then begin
      gain = 1.55
    endif
    if n_elements(readn) eq 0 then begin
      readn = 7.5
    endif
    if n_elements(time) eq 0 then begin
      time = hierarch(newhead, 'EXPTIME', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
    if n_elements(backg) eq 0 then begin
      drk = 6.*time/3600.           ; dark current 6e-/hour
      if drk lt 1. then drk = 0.    ; set zero for short exposures
      sky = 0.
      backg = drk + sky
    endif

; For object frames, prepare for heliocentric correction
    imtype=hierarch(newhead,'OBS-TYPE',count=count,errmes=errmes)
    if strpos(imtype,'object') ne -1 then begin
      if n_elements(ra2000) eq 0 then begin
        ra=hierarch(newhead,'RA-DEG',count=count, errmes=errmes)
      endif
      if n_elements(de2000) eq 0 then begin
        dec=hierarch(newhead,'DEC-DEG',count=count, errmes=errmes)
      endif
      if n_elements(jd) eq 0 then begin
        date_obs=hierarch(newhead,'DATE-OBS',count=count,errmes=errmes)
        date_obs=strsplit(date_obs,'-',/EXTRACT)
        ut=hierarch(newhead,'EXPSTART',count=count,errmes=errmes)
        ut = strsplit(ut, ':', /extract)
        tmid = ten(ut[0], ut[1], ut[2]) + time/2d0/3600d0
        jdcnv,date_obs[0],date_obs[1],date_obs[2],tmid,jd  ;covert UT to JD
        jd=jd-2400000.
      endif
      precess, ra, dec, float(date_obs[0])+float(date_obs[1])/12.+float(date_obs[2])/365., 2000.
      ra2000=ra/15.d0
      de2000=dec
      if n_elements(obslon) eq 0 then obslon = double(ten(17,53,37.9)) ;observatory coordinates
      if n_elements(obslat) eq 0 then obslat = double(ten(28,45,28.3))
      if n_elements(obsalt) eq 0 then obsalt = 2387.2
    endif

return
end

pro modeinfo_nes,mode,newhead,object,orient,xlo,xhi,ylo,yhi,gain,readn,backg,time $
                       ,obslat,obslon,obsalt,ra2000,de2000,jd,NONLINEAR=nonlinear
;===========================================================================================
;The Zelentchuk 6m NES spectrometer.
    if n_elements(orient) eq 0 then orient = 3
    if not keyword_set(xr) then xr = [3, 1991]
    xlo = xr(0)
    xhi = xr(1)
    if not keyword_set(yr) then yr = [0, 2019]
    ylo = yr(0)
    yhi = yr(1)
    if n_elements(gain) eq 0 then begin
      gain = 1.25
    endif
    if n_elements(readn) eq 0 then begin
      readn = hierarch(newhead, 'RDNOISE', count=count)
      if count eq 0 then readn = 7.7
    endif
    if n_elements(backg) eq 0 then begin
      drk = 1.
      sky = 0.
      backg = gain * (drk + sky)        ;convert ADU to electrons
    endif
    if n_elements(time) eq 0 then begin
      time = hierarch(newhead, 'EXPTIME', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif

    if n_elements(ra2000) eq 0 then begin
      ra=hierarch(newhead,'RA',count=count, errmes=errmes)
      ra=ra/3600.                        ; RA in hours
    endif

    if n_elements(de2000) eq 0 then begin
      dec=hierarch(newhead,'DEC',count=count, errmes=errmes)
      dec=dec/3600.                      ; DEC in degrees
    endif

    if n_elements(jd) eq 0 then begin
      date_obs=hierarch(newhead,'DATE-OBS',count=count,errmes=errmes)
      date_obs=strsplit(date_obs,'-',/EXTRACT)
      t1=hierarch(newhead,'TM_START',count=count,errmes=errmes)
      t2=hierarch(newhead,'TM_END'  ,count=count,errmes=errmes)
      if (t2 ge t1) then ut=(t1+t2)*0.5d0/3600.d0 $
      else begin
        ut = (86400L+t2+t1)*0.5/3600.d0
        if (ut ge 24) then begin
          date_obs=hierarch(newhead,'DATE',count=count,errmes=errmes)
          date_obs=strsplit(date_obs,'-',/EXTRACT)
          ut=ut-24
        endif
      endelse
      jdcnv,date_obs[0],date_obs[1],date_obs[2],ut,jd  ;covert UT to JD
      jd=jd-2400000.d0
    endif
    ra=ra*15.d0
    precess, ra, dec, float(date_obs[0])+float(date_obs[1])/12.+float(date_obs[2])/365., 2000.
    ra2000=ra/15.d0
    de2000=dec

    if n_elements(obslon) eq 0 then obslon = double(ten(-41,26.5)) ;observatory coordinates
    if n_elements(obslat) eq 0 then obslat = double(ten( 43,39.2))
    if n_elements(obsalt) eq 0 then obsalt = 2070.

return
end

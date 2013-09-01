pro modeinfo_boes,mode,newhead,object,orient,xlo,xhi,ylo,yhi,gain,readn,backg,time $
                       ,obslat,obslon,obsalt,ra2000,de2000,jd,NONLINEAR=nonlinear
;===========================================================================================
;BOES spectrometer in South Korea
    if n_elements(orient) eq 0 then orient = 1
    if not keyword_set(xr) then xr = [50, 2097]
    xlo = xr(0)
    xhi = xr(1)
    if not keyword_set(yr) then yr = [1, 4100]
    ylo = yr(0)
    yhi = yr(1)
    if n_elements(gain) eq 0 then begin
      gain = hierarch(newhead, 'GAIN', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
    if n_elements(readn) eq 0 then begin
      readn = hierarch(newhead, 'RDNOISE', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
    if n_elements(backg) eq 0 then begin
      drk = 0.
      sky = 0.
      backg = gain * (drk + sky)        ;convert ADU to electrons
    endif
    if n_elements(time) eq 0 then begin
      time = hierarch(newhead, 'EXPTIME', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif

    obj=hierarch(newhead,'OBJECT',count=count, errmes=errmes)
    obj=strtrim(obj,2)
    if(obj ne 'ThAr' and obj ne 'THL' and obj ne 'bias') then  begin
      if n_elements(ra2000) eq 0 then begin
        ra=hierarch(newhead,'RA',count=count, errmes=errmes)
        ra=float(strsplit(ra,':',/EXTRACT))
        ra=ten(ra[0],ra[1],ra[2]) ; RA in hours
        ra=ra/3600.
      endif

      if n_elements(de2000) eq 0 then begin
        dec=hierarch(newhead,'DEC',count=count, errmes=errmes)
        dec=float(strsplit(dec,':',/EXTRACT))
        dec=ten(dec[0],dec[1],dec[2])
        dec=dec/3600.
      endif

      if n_elements(jd) eq 0 then begin
        date_obs=hierarch(newhead,'DATE-OBS',count=count,errmes=errmes)
        date_obs=strsplit(date_obs,'-',/EXTRACT)
        ut=hierarch(newhead,'UT',count=count,errmes=errmes)
        ut=float(strsplit(ut,':',/EXTRACT))
        exp=hierarch(newhead,'EXPTIME',count=count,errmes=errmes)
        ut=ut(0)+ut(1)/60.d0+(ut(2)+exp*0.5d0)/3600.d0
        jdcnv,date_obs[0],date_obs[1],date_obs[2],ut,jd  ;covert UT to JD
        jd=jd-2400000.d0
      endif
      if n_elements(epoch) eq 0 then begin
        epoch=hierarch(newhead,'EPOCH',count=count,errmes=errmes)
      endif
      ra=ra*15.d0
      precess, ra, dec, epoch, 2000.
      ra2000=ra/15.d0
      de2000=dec

      if n_elements(obslon) eq 0 then obslon = double(ten(-128,58,35.7)) ;observatory coordinates
      if n_elements(obslat) eq 0 then obslat = double(ten(36,9,53.2))
      if n_elements(obsalt) eq 0 then obsalt = 1174.
    endif

return
end

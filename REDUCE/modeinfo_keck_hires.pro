pro modeinfo_keck_hires,mode,newhead,object,orient,xlo,xhi,ylo,yhi,gain,readn,backg,time $
                       ,obslat,obslon,obsalt,ra2000,de2000,jd,NONLINEAR=nonlinear
;===========================================================================================
;all modes of HIRES spectrometer at Keck.
;
; Find out the gainb mode
    gain=hierarch(newhead, 'CCDGAIN', count=count, errmes=errmes)
    if(strtrim(gain,2) eq 'low') then low_gain=1 else low_gain=0
;
; Get what they tell us about the usefull area
    datasec=hierarch(newhead, 'DATASEC', count=count, errmes=errmes) ; DATASEC = '[1:2048,1:4608]' format
    datasec=strmid(datasec,1,strlen(datasec)-2)
    datasec=byte(datasec)
    ii=where(datasec eq (byte(':'))[0], nii)
    if(nii gt 0) then datasec(ii)=byte(',')
    datasec=string(datasec)
    xr_ccd=intarr(2)
    yr_ccd=intarr(2)
    reads,datasec,xr_ccd,yr_ccd

    if mode eq 'keck_hires_blue' then begin ; blue chip
      reorient = 3
      if(low_gain) then begin
; Low gain
        gain  = 1.9   ;e-/ADU
        readn = 2.8   ;e-
        dark  = 3.8   ;e-/hour
      endif else begin
; High gain
        gain  = 0.78  ;e-/ADU
        readn = 2.6   ;e-
        dark  = 5.5   ;e-/hour
      endelse
; Extraction area
      if keyword_set(xr) then begin
        xlo = xr(0)
        xhi = xr(1)
      endif else begin
        xlo = xr_ccd[0]+29
        xhi = xr_ccd[1]-1
      endelse
      if keyword_set(yr) then begin
        ylo = yr(0)
        yhi = yr(1)
      endif else begin
        ylo = yr_ccd[0]-1
        yhi = yr_ccd[1]-51
      endelse
    endif
    if mode eq 'keck_hires_green' then begin ; middle chip
      reorient = 3
      if(low_gain) then begin
; Low gain
        gain  = 2.2   ;e-/ADU
        readn = 3.1   ;e-
        dark  = 4.4   ;e-/hour
      endif else begin
; High gain
        gain  = 0.86  ;e-/ADU
        readn = 2.6   ;e-
        dark  = 7.7   ;e-/hour
      endelse
; Extraction area
      if keyword_set(xr) then begin
        xlo = xr(0)
        xhi = xr(1)
      endif else begin
        xlo = xr_ccd[0]+29
        xhi = xr_ccd[1]-1
      endelse
      if keyword_set(yr) then begin
        ylo = yr(0)
        yhi = yr(1)
      endif else begin
        ylo = yr_ccd[0]-1
        yhi = yr_ccd[1]-76
      endelse
    endif
    if mode eq 'keck_hires_red' then begin  ; red chip
      reorient = 3
      if(low_gain) then begin
; Low gain
        gain  = 2.2   ;e-/ADU
        readn = 3.1   ;e-
        dark  = 2.2   ;e-/hour
      endif else begin
; High gain
        gain  = 0.84  ;e-/ADU
        readn = 2.9   ;e-
        dark  = 4.2   ;e-/hour
      endelse
; Extraction area
      if keyword_set(xr) then begin
        xlo = xr(0)
        xhi = xr(1)
      endif else begin
        xlo = xr_ccd[0]+29
        xhi = xr_ccd[1]-1
      endelse
      if keyword_set(yr) then begin
        ylo = yr(0)
        yhi = yr(1)
      endif else begin
        ylo = yr_ccd[0]-1
        yhi = yr_ccd[1]-51
      endelse
    endif
    if n_elements(orient) eq 0 then orient = reorient
    if n_elements(time) eq 0 then begin
      time = hierarch(newhead, 'EXPTIME', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
    backg = dark*time+readn
; For object frames, prepare for heliocentric correction
    lamp  =hierarch(newhead,'LAMPNAME',count=count,errmes=errmes)
    target=hierarch(newhead,'TARGNAME',count=count,errmes=errmes)
    if(strtrim(lamp,2) eq 'none' and strtrim(target,2) ne ' horizon lock') then begin
      obslat=ten( 19, 49, 33.40757d0)
      obslon=ten(155, 28, 28.98665d0)
      obsalt=4159.6d0
      if n_elements(ra2000) eq 0 then begin
        ra = hierarch(newhead,'RA',count=count, errmes=errmes)
        ra=float(strsplit(ra,':',/EXTRACT))
        ra=ten(ra[0],ra[1],ra[2])
        ra=ra/3600.
      endif
      if n_elements(de2000) eq 0 then begin
        dec = hierarch(newhead,'DEC',count=count, errmes=errmes)
        dec=float(strsplit(dec,':',/EXTRACT))
        dec=ten(dec[0],dec[1],dec[2])
        dec=dec/3600.
      endif
      if n_elements(epoch) eq 0 then begin
        epoch=hierarch(newhead,'EQUINOX',count=count,errmes=errmes)
      endif
      ra=ra*15.d0
      precess, ra, dec, epoch, 2000.d0
      ra2000=ra/15.d0
      de2000=dec
      if n_elements(jd) eq 0 then begin
;        b=hierarch(newhead,'DATE_BEG',count=count,errmes=errmes)
;        date_obs=strmid(b,0,strpos(b,'T'))
;        ut=strmid(b,strpos(b,'T')+1)
;        date_obs=strsplit(date_obs,'-',/EXTRACT)
;        ut=float(strsplit(ut,':',/EXTRACT))
;        ut=ut(0)+ut(1)/60.d0+(ut(2))/3600.d0
;        jdcnv,date_obs[0],date_obs[1],date_obs[2],ut,jd  ;covert UT to JD
;        jd_beg=jd-2400000.d0
;        e=hierarch(newhead,'DATE_END',count=count,errmes=errmes)
;        date_obs=strmid(e,0,strpos(e,'T'))
;        ut=strmid(e,strpos(e,'T')+1)
;        date_obs=strsplit(date_obs,'-',/EXTRACT)
;        ut=float(strsplit(ut,':',/EXTRACT))
;        ut=ut(0)+ut(1)/60.d0+(ut(2))/3600.d0
;        jdcnv,date_obs[0],date_obs[1],date_obs[2],ut,jd  ;covert UT to JD
;        jd_end=jd-2400000.d0
;        jd=0.5d0*(jd_beg+jd_end)
        mjd = hierarch(newhead,'MJD',count=count,errmes=errmes)+time/2./3600./24.+0.5d0
        jd = mjd
      endif
    endif
  return
end

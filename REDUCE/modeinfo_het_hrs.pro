pro modeinfo_het_hrs,mode,newhead,object,orient,xlo,xhi,ylo,yhi,gain,readn,backg,time $
                       ,obslat,obslon,obsalt,ra2000,de2000,jd,NONLINEAR=nonlinear
;===========================================================================================
;The HET HRS spectrometer (red CCD in the mosaic, extent=1). Tested for 2x5 binning only
    trimsec= hierarch(newhead, 'TRIMSEC', count=count, errmes=errmes) ; data rectangle in the frame
    if(count eq 1) then begin
      trimsec=strtrim(trimsec,2)
      trimsec=byte(trimsec)
      ii=where(trimsec eq (byte(':'))[0], nii)
      if(nii gt 0) then trimsec(ii)=byte(',')
      trimsec=string(trimsec[1:n_elements(trimsec)-2])
      xr_ccd=intarr(2)
      yr_ccd=intarr(2)
      reads,trimsec,xr_ccd,yr_ccd
    endif
    if mode eq 'het_hrs_blue' then begin ; blue CCD
      reorient = 3
    endif
    if mode eq 'het_hrs_red' then begin  ; red CCD
      reorient = 3
    endif
    if n_elements(orient) eq 0 then orient = reorient
    if not keyword_set(xr) then begin
      xlo = xr_ccd[0]
      xhi = xr_ccd[1]
    endif else begin
      xlo = xr[0]
      xhi = xr[1]
    endelse
    if not keyword_set(yr) then begin
      ylo = yr_ccd[0]
      yhi = yr_ccd[1]
    endif else begin
      ylo = yr[0]
      yhi = yr[1]
    endelse

    gain=1.
    backg=0.
    readn=0.
    time=0
; For object frames, read CCD settings and prepare for heliocentric correction
    imtype=hierarch(newhead,'OBJECT',count=count,errmes=errmes)
    if strmatch(imtype,'object*',/fold_case) then begin
      if n_elements(gain) eq 0 then begin
        if mode eq 'het_hrs_red' then begin            ; blue CCD
          gain = hierarch(newhead, 'GAIN1', count=count, errmes=errmes)
        endif else if mode eq 'het_hrs_blue' then begin ; red CCD
          gain = hierarch(newhead, 'GAIN3', count=count, errmes=errmes)
        endif
        if count eq 0 then message, errmes
      endif
      if n_elements(readn) eq 0 then begin
        if mode eq 'het_hrs_red' then begin            ; blue CCD
          readn = hierarch(newhead, 'RDNOISE1', count=count, errmes=errmes)
        endif else if mode eq 'het_hrs_blue' then begin ; red CCD
          readn = hierarch(newhead, 'RDNOISE3', count=count, errmes=errmes)
        endif
      endif
      backg = gain * 0.1        ;convert ADU to electrons
      if n_elements(time) eq 0 then begin
        time = hierarch(newhead, 'EXPTIME', count=count, errmes=errmes)
        if count eq 0 then message, errmes
      endif
      if n_elements(ra2000) eq 0 then begin
        ra = hierarch(newhead,'RA',count=count, errmes=errmes)
        ra=float(strsplit(ra,':',/EXTRACT))
        ra=ten(ra[0],ra[1],ra[2]) ; RA in hours
        ra=ra/3600.
      endif
      if n_elements(de2000) eq 0 then begin
        dec = hierarch(newhead,'DEC',count=count, errmes=errmes)
        dec=float(strsplit(dec,':',/EXTRACT))
        dec=ten(dec[0],dec[1],dec[2])
        dec=dec/3600.
      endif
      if n_elements(epoch) eq 0 then begin
        epoch=hierarch(newhead,'EPOCH',count=count,errmes=errmes)
      endif
      ra=ra*15.d0
      precess, ra, dec, epoch, 2000.
      ra2000=ra/15.d0
      de2000=dec
      if n_elements(jd) eq 0 then begin
        jd = hierarch(newhead,'MJD',count=count,errmes=errmes)
      endif
;observatory coordinates
      if n_elements(obslon) eq 0 then obslon = ten(104,01,30) ;observatory coordinates, seconds are unknown
      if n_elements(obslat) eq 0 then obslat = ten( 30,41,30)
      if n_elements(obsalt) eq 0 then obsalt = 2002.
    endif

return
end

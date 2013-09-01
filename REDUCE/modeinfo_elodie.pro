pro modeinfo_elodie,mode,newhead,object,orient,xlo,xhi,ylo,yhi,gain,readn,backg,time $
                       ,obslat,obslon,obsalt,ra2000,de2000,jd,NONLINEAR=nonlinear
;===========================================================================================
;Observatoire de Haute-Province (OHP) ELODIE.
    if n_elements(orient) eq 0 then orient = 0
    if not keyword_set(xr) then xr = [0, 1023]
    xlo = 0
;    xhi = hierarch(newhead, 'NAXIS1', count=count, errmes=errmes)-1
    xhi=1023
    if not keyword_set(yr) then yr = [0, 1023]
    ylo = 0
    yhi = hierarch(newhead, 'NAXIS2', count=count, errmes=errmes)-1
    if n_elements(gain) eq 0 then begin
      gain = hierarch(newhead, 'CCDGAIN', count=count, errmes=errmes)
      if count eq 0 then gain=2.65
    endif
    if n_elements(readn) eq 0 then begin
      readn = 8.5
    endif
    if n_elements(backg) eq 0 then begin
      drk = 1.
      sky = 0.
      backg = gain * (drk + sky)        ;convert ADU to electrons
    endif
    if n_elements(time) eq 0 then begin
      time = hierarch(newhead, 'EXPTIME', count=count, errmes=errmes)
      if count eq 0 then time=0.
    endif

    imtype=hierarch(newhead,'IMATYP',count=count,errmes=errmes) ;check if the frame contains
    if strmatch(imtype,'OBJ*') then begin                      ;a stellar spectrum

      if n_elements(ra2000) eq 0 then begin
        ra=hierarch(newhead,'ALPHA',count=count, errmes=errmes)
        count = 0                                              ;parse RA string
        value=0d
        while strpos(ra, ':') ge 0 do begin
          value=value+float(ra)/60.^count
          ra=strmid(ra,strpos(ra,':')+1)
          count = count + 1
        endwhile
        ra2000=value
      endif

      if n_elements(de2000) eq 0 then begin
        de=hierarch(newhead,'DELTA',count=count, errmes=errmes)
        count = 0                                              ;parse DE string
        sig=1.
        value=0d
        while strpos(de, ':') ge 0 do begin
          value=value+float(de)*sig/60.^count
          if count eq 0 and value lt 0 then sig=-1
          de=strmid(de,strpos(de,':')+1)
          count = count + 1
        endwhile
        de2000=value
      endif

      if n_elements(jd) eq 0 then begin
        year=hierarch(newhead,'TIME1',count=count,errmes=errmes)
        month=hierarch(newhead,'TIME2',count=count,errmes=errmes)
        day=hierarch(newhead,'TIME3',count=count,errmes=errmes)
        ut=hierarch(newhead,'TIME4',count=count,errmes=errmes)
        jdcnv,year,month,day,ut,jd                            ;covert UT to JD
        jd=jd-2400000.
      endif

      if n_elements(obslon) eq 0 then obslon = double(ten(-5,42.8)) ;observatory coordinates
      if n_elements(obslat) eq 0 then obslat = double(ten(43,55.9))
      if n_elements(obsalt) eq 0 then obsalt = 665.
    endif

return
end

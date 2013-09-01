pro modeinfo_mcd2d,mode,newhead,object,orient,xlo,xhi,ylo,yhi,gain,readn,backg,time $
                       ,obslat,obslon,obsalt,ra2000,de2000,jd,NONLINEAR=nonlinear
;===========================================================================================
;The McDonald 2dcoude at R=240,000 focus.
  if mode eq 'mcd2d_cs21' then begin
    if n_elements(orient) eq 0 then orient = 0
    if not keyword_set(xr) then xr = [2, 2046]
    xyr = hierarch(newhead, 'TRIMSEC', count=count, errmes=errmes)
    if count eq 0 then begin
      xlo =    2
      xhi = 2048
      ylo =    2
      yhi = 2047
    endif else begin
      xyr=strsplit(strmid(xyr,1),'[:,\,]',/REGEX)
      xyr=fix(strmid(xyr,strsplit(strmid(xyr,1),'[:,\,]',/REGEX)+1))
      xlo=xyr(0)+1
      xhi=xyr(1)
      ylo=xyr(2)
      yhi=xyr(3)
    endelse
    if n_elements(gain) eq 0 then begin ; Electrons per ADU
      gain = hierarch(newhead, 'GAIN', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
    if n_elements(readn) eq 0 then begin ; Electrons
      readn = hierarch(newhead, 'RDNOISE', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
    if n_elements(backg) eq 0 then begin
      drk = 1.
      sky = 0.
      backg = gain * (drk + sky)        ;convert ADU to electrons
    endif
    if n_elements(time) eq 0 then begin
      time = hierarch(newhead, 'DARKTIME', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
  endif

;===========================================================================================
;The McDonald 2dcoude at R=60,000 focus.
  if mode eq 'mcd2d_cs23' then begin
    if n_elements(orient) eq 0 then orient = 0
    if not keyword_set(xr) then xr = [1, 2046]
    xyr = hierarch(newhead, 'TRIMSEC', count=count, errmes=errmes)
    if count eq 0 then begin
      xlo =    2
      xhi = 2048
      ylo =    2
      yhi = 2047
    endif else begin
      xyr=fix(strmid(xyr,strsplit(strmid(xyr,1),'[:,\,]',/REGEX)+1))
      xlo=xyr(0)+1
      xhi=xyr(1)
      ylo=xyr(2)
      yhi=xyr(3)
    endelse
    if n_elements(gain) eq 0 then begin ; Electrons per ADU
      gain = hierarch(newhead, 'GAIN*', count=count, errmes=errmes)
      if count eq 0 then gain=1. $ ; message, errmes $
      else gain=gain[0].(0)
    endif
    if n_elements(readn) eq 0 then begin ; Electrons
      readn = hierarch(newhead, 'RDNOISE*', count=count, errmes=errmes)
      if count eq 0 then readn=0. $ ; message, errmes $
      else readn=readn[0].(0)
    endif
    if n_elements(backg) eq 0 then begin
      drk = 1.
      sky = 0.
      backg = gain * (drk + sky)        ;convert ADU to electrons
    endif
    if n_elements(time) eq 0 then begin
      time = hierarch(newhead, 'DARKTIME', count=count, errmes=errmes)
      if count eq 0 then time=0. ; message, errmes
    endif

; For object frames, prepare for heliocentric correction
    imtype=hierarch(newhead,'IMAGETYP',count=count,errmes=errmes)
    if strmatch(imtype,'object*') then begin
      if n_elements(ra2000) eq 0 then begin
        ra = hierarch(newhead,'RA',count=count, errmes=errmes)
        ra = strsplit(ra, ':', /extract)
        ra2000 = ten(ra[0], ra[1], ra[2])
      endif
      if n_elements(de2000) eq 0 then begin
        de = hierarch(newhead,'DEC',count=count, errmes=errmes)
        de = strsplit(de, ':', /extract)
        de2000 = ten(de[0], de[1], de[2])
      endif
      if n_elements(jd) eq 0 then begin
        dateobs = hierarch(newhead,'DATE-OBS',count=count,errmes=errmes)
        dateobs = strsplit(dateobs,'-',/extract)
        ut = hierarch(newhead,'UT',count=count,errmes=errmes)
        ut = strsplit(ut, ':', /extract)
        tmid = ten(ut[0], ut[1], ut[2]) + time/2d0
        jd=julday(dateobs[1],dateobs[2],dateobs[0],tmid)  ;MM,DD,YY,decimal HH to JD
        jd=jd-2400000.
      endif
      if n_elements(obslon) eq 0 then obslon = ten(104,1,18) ;observatory coordinates
      if n_elements(obslat) eq 0 then obslat = ten(30,40,18)
      if n_elements(obsalt) eq 0 then obsalt = 2075.
    endif
  endif

return
end

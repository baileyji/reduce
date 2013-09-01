pro modeinfo_not_fies,mode,newhead,object,reorient,xlo,xhi,ylo,yhi,gain,readn,backg,time $
                  ,obslat,obslon,obsalt,ra2000,de2000,jd,NONLINEAR=nonlinear
;===========================================================================================
;all modes of FIES spectrometer.
  if strpos(mode,'not_fies') ge 0 then begin
    reorient = 1
    if n_elements(orient) eq 0 then orient = reorient
    if not keyword_set(xr) then begin
      nxd = hierarch(newhead, 'NAXIS1')
      presc = 50
      ovrsc = 100
      xlo = presc
      xhi = nxd - ovrsc - 1
    endif else begin
      xlo = xr(0)
      xhi = xr(1)
    endelse
    if not keyword_set(yr) then begin
      nyd = hierarch(newhead, 'NAXIS2')
      presc = 0
      ovrsc = 1
      ylo = presc
      yhi = nyd - ovrsc - 1
    endif else begin
      ylo = yr(0)
      yhi = yr(1)
    endelse
    if n_elements(gain) eq 0 then begin
      gain = hierarch(newhead, 'GAIN', count=count, errmes=errmes)
      if count eq 0 then gain=0.765
    endif
    if n_elements(readn) eq 0 then begin
      readn = hierarch(newhead, 'RDNOISE', count=count, errmes=errmes)
      if count eq 0 then readn=3.0
    endif
    if n_elements(backg) eq 0 then begin
      drk = hierarch(newhead, 'DARK', count=count, errmes=errmes)
      if count eq 0 then drk=0.
      backg = gain * drk          ;convert ADU to electrons
    endif
    if n_elements(time) eq 0 then begin
      time = hierarch(newhead, 'EXPTIME', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
; For object frames, prepare for heliocentric correction
    imtype=hierarch(newhead,'OBJECT',count=count,errmes=errmes)
    if(strpos(imtype,'ThAr') eq -1 and $
       strpos(imtype,'bias') eq -1 and $
       strpos(imtype,'flat') eq -1 and $
       strpos(imtype,'fluxtest') eq -1) then begin
      ra   = hierarch(newhead,'RA')   ; Catalogue RA in decimal deg.
      dec  = hierarch(newhead,'DEC')  ; Catalogue DEC in decimal deg
      ut   = hierarch(newhead,'UT')   ; UT of the start
      date = hierarch(newhead,'DATE-OBS') ; Start of observation
      date = float(strsplit(date,'[-,T,:]',/REGEX,/EXTRACT))
      equinox = date[0] +(date[1]+(date[2]+(date[3]+(date[4]+date[5]/60.d0)/60.d0)/24.d0)/365.25d0)/12.d0
      juldate, date, jd
      ra = ra*15.d0
      precess, ra, dec, equinox, 2000.d0
      ra2000 = ra/15.d0
      de2000 = dec
;observatory coordinates
      obslon = double(ten(17,53,06.3)) ; GPS longitude
      obslat = double(ten(28,45,26.2)) ; GPS latitude
      obsalt = 2382.d0
    endif
  endif
  return
end

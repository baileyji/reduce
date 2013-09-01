pro modeinfo_ces,mode,newhead,object,orient,xlo,xhi,ylo,yhi,gain,readn,backg,time $
                       ,obslat,obslon,obsalt,ra2000,de2000,jd,NONLINEAR=nonlinear
;===========================================================================================
;The ESO 3.6m CES spectrometer.
    if n_elements(orient) eq 0 then orient = 3
    if not keyword_set(xr) then xr = [51, 1024+50]
    xlo = hierarch(newhead, 'HIERARCH ESO DET OUT2 PRSCX', count=count, errmes=errmes)
    if count eq 0 then xlo = xr(0)-1
    xhi = hierarch(newhead, 'NAXIS1', count=count, errmes=errmes)
    ovrsc = hierarch(newhead, 'HIERARCH ESO DET OUT2 OVSCX', count=count, errmes=errmes)
    if count ne 0 then xhi = xhi - ovrsc else xhi = xr(1)-1
    if not keyword_set(yr) then yr = [0, 4095]
    ylo = 0
    yhi = hierarch(newhead, 'NAXIS2', count=count, errmes=errmes)-1
;    ylo = yr(0)
;    yhi = yr(1)
    if n_elements(gain) eq 0 then begin
      gain = hierarch(newhead, 'HIERARCH ESO DET OUT2 CONAD', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
    if n_elements(readn) eq 0 then begin
      readn = hierarch(newhead, 'HIERARCH ESO DET OUT2 RON', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
    if n_elements(backg) eq 0 then begin
      drk = hierarch(newhead, 'HIERARCH ESO INS DET1 OFFDRK', count=count, errmes=errmes)
      if count eq 0 then message, errmes, /CONTINUE
      sky = hierarch(newhead, 'HIERARCH ESO INS DET1 OFFSKY', count=count, errmes=errmes)
      if count eq 0 then message, errmes, /CONTINUE
      backg = gain * (drk + sky)        ;convert ADU to electrons
    endif
    if n_elements(time) eq 0 then begin
      time = hierarch(newhead, 'EXPTIME', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif

return
end

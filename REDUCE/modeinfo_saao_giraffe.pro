pro modeinfo_saao_giraffe,mode,newhead,object,orient,xlo,xhi,ylo,yhi,gain,readn,backg,time $
                       ,obslat,obslon,obsalt,ra2000,de2000,jd,NONLINEAR=nonlinear
;===========================================================================================
;The SAAO Giraffe spectrometer (South Africa).
  if mode eq 'saao_giraffe_top' then begin
    if n_elements(orient) eq 0 then orient = 2
    if not keyword_set(xr) then xr = [52, 1073]
    xlo = xr(0)
    xhi = xr(1)
    if not keyword_set(yr) then yr = [512, 1023]
    ylo = yr(0)
    yhi = yr(1)
    if n_elements(gain) eq 0 then begin
      gain = 1./0.833
    endif
    if n_elements(readn) eq 0 then begin
      readn = 6.7
    endif
    if n_elements(backg) eq 0 then begin
      drk = 1.
      sky = 0.
      backg = gain * (drk + sky)        ;convert ADU to electrons
    endif
    if n_elements(time) eq 0 then begin
      time = hierarch(newhead, 'ITIME', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
  endif

  if mode eq 'saao_giraffe_bottom' then begin
    if n_elements(orient) eq 0 then orient = 2
    if not keyword_set(xr) then xr = [52, 1073]
    xlo = xr(0)
    xhi = xr(1)
    if not keyword_set(yr) then yr = [0, 511]
    ylo = yr(0)
    yhi = yr(1)
    if n_elements(gain) eq 0 then begin
      gain = 1./0.926
    endif
    if n_elements(readn) eq 0 then begin
      readn = 6.1
    endif
    if n_elements(backg) eq 0 then begin
      drk = 1.
      sky = 0.
      backg = gain * (drk + sky)        ;convert ADU to electrons
    endif
    if n_elements(time) eq 0 then begin
      time = hierarch(newhead, 'ITIME', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
  endif

return
end

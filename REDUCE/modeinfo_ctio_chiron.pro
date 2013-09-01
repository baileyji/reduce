function nonlinear_ctio_chiron,image,header,AMPLIFIER=amplifier,GAIN=gain

; Readout speed affects non-linearity correction
  rdspeed = strtrim(sxpar(header, 'SPEEDMOD', count=count), 2)

; Check the number of amplifiers used
  ver = sxpar(header, 'E_HVERS', count=count)
  if(ver gt 1.001) then begin
    n_amp = sxpar(header, 'E_AMPL', count=count)
  endif else begin
    n_amp = 1
  endelse

  if(n_amp eq 2) then nlc = (rdspeed eq 'fast')? [5.0d-6,4.3d-6]:[4.5d-6,4.0d-6]
  
; For image fragment user must specify the amplifier and the gain
; Only one set of values is allowed as we do not have geometrical
; reference. This type of call is made from inside sumfits.
  if(keyword_set(amplifier) and keyword_set(gain)) then begin
    return, image * (1.0 - nlc[amplifier-1] * image) * gain
  endif

; Alternative mode: process the whole image at once
  if(n_amp gt 1) then begin
    image=float(image)
    for amp=1,n_amp do begin
      sf = suffix(amp, 1)
      xlo  = sxpar(header, 'E_XLO' +sf, count=count)
      ylo  = sxpar(header, 'E_YLO' +sf, count=count)
      xhi  = sxpar(header, 'E_XHI' +sf, count=count)
      yhi  = sxpar(header, 'E_YHI' +sf, count=count)
      gain = sxpar(header, 'E_GAIN'+sf, count=count)
      image[xlo:xhi,ylo:yhi] = image[xlo:xhi,ylo:yhi] $
                             * (1.0 - nlc[amp-1] * image[xlo:xhi,ylo:yhi]) * gain
    endfor
  endif

  return,image
end

pro modeinfo_ctio_chiron,mode,newhead,object,orient,xlo,xhi,ylo,yhi,gain,readn,backg,time $
                        ,obslat,obslon,obsalt,ra2000,de2000,jd,NONLINEAR=nonlinear
;===========================================================================================
;The CTIO Chiron spectrometer (2-amplifier mode)

    left  = '21' ; Left section of the CCD
    right = '22' ; Right section of the CCD
    reorient = 6 ; Wavelength runs right to left and top to bottom
    if n_elements(orient) eq 0 then orient = reorient

    dsec_left = hierarch(newhead, 'DSEC'+left, count=count, errmes=errmes) ; left data section
    if count eq 0 then message, errmes
    i1 = strpos(dsec_left,'[') & i2 = strpos(dsec_left,']')
    dsec_left = strsplit(strmid(dsec_left, i1+1, i2-i1-1),'[:,]',/EXTRACT)
    x0_left = 0L & x1_left = 0L & y0_left = 0L & y1_left=0L 
    reads,dsec_left,x0_left,x1_left,y0_left,y1_left 

    dsec_right = hierarch(newhead, 'DSEC'+right, count=count, errmes=errmes) ; left data section
    if count eq 0 then message, errmes
    i1 = strpos(dsec_right,'[') & i2 = strpos(dsec_right,']')
    dsec_right = strsplit(strmid(dsec_right, i1+1, i2-i1-1),'[:,]',/EXTRACT)
    x0_right = 0L & x1_right = 0L & y0_right = 0L & y1_right=0L 
    reads,dsec_right,x0_right,x1_right,y0_right,y1_right 

    xlo = [x0_left, x0_right] - 1
    xhi = [x1_left, x1_right] - 1
    ylo = [y0_left+1, y0_right+1] - 1 ; Takes care of funny increase in bias
    yhi = [y1_left-4, y1_right-4] - 1 ; at the top and the bottom

    if keyword_set(xr) then begin
      xlo = xr(0)
      xhi = xr(1)
    endif
    if keyword_set(yr) then begin
      ylo = yr(0)
      yhi = yr(1)
    endif

    if n_elements(gain) eq 0 then begin
      gain_left  = hierarch(newhead, 'GAIN'+left,  count=count, errmes=errmes)
      gain_right = hierarch(newhead, 'GAIN'+right, count=count, errmes=errmes)
      gain = [gain_left, gain_right]
      if count eq 0 then message, errmes
    endif

    if n_elements(readn) eq 0 then begin
      readn_left  = hierarch(newhead, 'RON'+left,  count=count, errmes=errmes)
      readn_right = hierarch(newhead, 'RON'+right, count=count, errmes=errmes)
      readn = [readn_left, readn_right]
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

; For object frames, prepare for heliocentric correction
    imtype=hierarch(newhead,'OBJECT',count=count,errmes=errmes)
    if(imtype ne 'quartz' and imtype ne 'bias' and imtype ne 'ThAr') then begin
      if n_elements(ra2000) eq 0 then begin
        ra = hierarch(newhead,'RA',count=count, errmes=errmes)
        ra = float(strsplit(ra,':',/EXTRACT))
        ra = ten(ra[0],ra[1],ra[2])
        ra2000 = ra/3600.
      endif
      if n_elements(de2000) eq 0 then begin
        dec = hierarch(newhead,'DEC',count=count, errmes=errmes)
        dec = float(strsplit(dec,':',/EXTRACT))
        dec = ten(dec[0],dec[1],dec[2])
        de2000 = dec/3600.
      endif
      if n_elements(jd) eq 0 then begin
        b = hierarch(newhead, 'DATE-OBS', count=count, errmes=errmes)
        date_obs = strmid(b, 0, strpos(b,'T'))
        ut =strmid(b, strpos(b, 'T') + 1)
        date_obs = strsplit(date_obs, '-', /EXTRACT)
        ut = float(strsplit(ut, ':', /EXTRACT))
        ut = ut[0] + ut[1] / 60.d0 + ut[2] / 3600.d0
        jdcnv,date_obs[0], date_obs[1], date_obs[2], ut, jd  ;convert UT to JD
        jd = jd - 2400000.d0 + time * 0.5d0
      endif
;observatory coordinates
      if n_elements(obslon) eq 0 then obslon =  70.80627d0
      if n_elements(obslat) eq 0 then obslat = -30.16908d0
      if n_elements(obsalt) eq 0 then obsalt = 2216.d0
    endif
  nonlinear='nonlinear_ctio_chiron'
  return
end

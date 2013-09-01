pro modeinfo_subaru_hds,mode,newhead,object,orient,xlo,xhi,ylo,yhi,gain,readn,backg,time $
                       ,obslat,obslon,obsalt,ra2000,de2000,jd,NONLINEAR=nonlinear
;===========================================================================================
;The Subaru HDS spectrometer.
;Acoording to the manual as of June-2008 HDS contains two
;2k x 4k CCDs referred in instrument mode as subaru_hds_blue and
;subaru_hds_red. Each CCD is (normally) read using 2 amplifiers
;in two opposite ends of a short side. The read is performed with
;overscan creating what appears to be a gap between the orders.
;The two read-outs have different gain and bias thus are better
;treated as separate parts.
;
  case mode of
'subaru_hds_blue_left':
'subaru_hds_red_left' :  begin
      biny = hierarch(newhead, 'BIN-FCT1' $  ; binning factor in x-dispersion direction
                  , count=count, errmes=errmes)
      binx = hierarch(newhead, 'BIN-FCT2' $  ; binning factor in dispersion direction
                  , count=count, errmes=errmes)
      xlo = 0                               ; starting column
      xhi = hierarch(newhead, 'H_OSMIN1' $  ; ending column of the frame
                  , count=count, errmes=errmes)
      xhi = xhi - 2
      ylo = hierarch(newhead, 'H_OSMIN2' $  ; starting row of the frame
                  , count=count, errmes=errmes)
      ylo = ylo + 1 
      yhi = hierarch(newhead, 'H_OSMAX2' $  ; ending row of the frame
                  , count=count, errmes=errmes)
      yhi = yhi - 3 
      gain = hierarch(newhead, 'H_GAIN1' $  ; gain at the left part of the frame
                  , count=count, errmes=errmes)
    end
'subaru_hds_blue_right':
'subaru_hds_red_right' : begin
      biny = hierarch(newhead, 'BIN-FCT1' $  ; binning factor in x-dispersion direction
                  , count=count, errmes=errmes)
      binx = hierarch(newhead, 'BIN-FCT2' $  ; binning factor in dispersion direction
                  , count=count, errmes=errmes)
      xlo = hierarch(newhead, 'H_OSMAX2' $  ; starting column of the frame
                  , count=count, errmes=errmes)
      xhi = hierarch(newhead, 'NAXIS1' $  ; full X-size of the whole frame
                  , count=count, errmes=errmes)
      xhi = xhi - 1
      ylo = hierarch(newhead, 'H_OSMIN2' $  ; starting row of the frame
                  , count=count, errmes=errmes)
      ylo = ylo + 1 
      yhi = hierarch(newhead, 'H_OSMAX2' $  ; ending row of the frame
                  , count=count, errmes=errmes)
      yhi = yhi - 3 
      gain = hierarch(newhead, 'H_GAIN2' $  ; gain at the left part of the frame
                  , count=count, errmes=errmes)
     end
default: return
    endcase
;
    reorient = 3
    readn = 4.4 ; Read-out noise, why the hell it is not in the header?
    readn = readn / gain ; convert to ADU
    dark = 10. ; dark current is 10 e-/hour
    dark = dark / gain / 3600.  ; convert to ADU/s
;
    if n_elements(orient) eq 0 then orient = reorient
;
    if n_elements(time) eq 0 then begin
      time = hierarch(newhead, 'EXPTIME', count=count, errmes=errmes)
      if count eq 0 then message, errmes
    endif
    backg = readn + time * dark  ;total noise
; For object frames, prepare for heliocentric correction
    imtype=hierarch(newhead,'OBJECT',count=count,errmes=errmes)
    if strmatch(imtype,'object*',/fold_case) then begin
      if n_elements(ra2000) eq 0 then begin
        ra2000 = hierarch(newhead,'RA2000',count=count, errmes=errmes)
        ra2000 = ra2000/15.
      endif
      if n_elements(de2000) eq 0 then begin
        de2000 = hierarch(newhead,'DEC2000',count=count, errmes=errmes)
      endif
      if n_elements(jd) eq 0 then begin
        jd = hierarch(newhead,'MJD',count=count,errmes=errmes)
      endif
; Observatory coordinates
      if n_elements(obslon) eq 0 then obslon = ten(115,28,50) ;observatory coordinates
      if n_elements(obslat) eq 0 then obslat = ten(19,49,43)
      if n_elements(obsalt) eq 0 then obsalt = 4139.
    endif
;
  return
end

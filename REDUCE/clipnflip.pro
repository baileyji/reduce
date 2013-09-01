function clipnflip, image, header, xr=xr, yr=yr, orient=orient
;+
;Process an image and associated FITS header already in memory as follows:
; 1) Trim image to desired subregion: newimage = image(xlo:xhi,ylo:yhi)
; 2) Transform to standard orientation (red at top, orders run left to right)
;
;Input:
; image (array) Raw image to be processed.
;
;Input/Output:
; [header] (string array) FITS header associated with image. Keyword arguments
;   xr=, yr=, and orient= need not be specified if header contains the
;   corresponding cards E_XLO, E_XHI, E_YLO, E_YHI, and E_ORIENT.
; [xr=] (vector(2)) first and last column to retain in input image.
; [yr=] (vector(2)) first and last row to retain in input image.
; [orient=] (integer) flag indicating how trimmed image should be reoriented.
;
;Output:
; Return value of function is either the processed image (array, reoriented
;   and trimmed) or an error message (scalar string).
;
;History:
; 2000-Aug-25 Valenti  Adapted from getimage.pro.
;-

if n_params() lt 1 then begin
  print,'syntax: newimage = clipnflip(image [,header ,xr= ,yr= ,orient=])'
  return, 'Syntax error'
end

;Check that required arguments have proper structure.
  err = ''
  sz = size(image)
  if sz(0) ne 2 then err = 'Image must be a two-dimensional array'
  if n_elements(header) gt 0 then begin     ;header specified
    sz = size(header)
    if sz(0) ne 1 or sz(2) ne 7 then begin
      message, /info, 'Header must be a string vector'
      err = 'Bad argument passed to trimflip'
    endif
  endif else begin              ;keywords only
    if not keyword_set(xr) or not keyword_set(yr) $
                           or n_elements(orient) eq 0 then begin
      message, /info, 'Must specify either HEADER or {XR=, YR=, and ORIENT=}'
      err = 'Bad keyword argument passed to trimflip'
    endif else begin
      if n_elements(xr) ne 2 then begin
        message, /info, 'XR= keyword must be 2 element vector'
        err = 'Bad keyword argument passed to trimflip'
      endif
      if n_elements(yr) ne 2 then begin
        message, /info, 'YR= keyword must be 2 element vector'
        err = 'Bad keyword argument passed to trimflip'
      endif
    endelse
  endelse
  if keyword_set(err) then return, err

;;Check the number of regions
;  n_regions = sxpar(header, 'E_NREG', count=count)
;  if count eq 0 then n_regions=1 else n_regions=n_regions>1

;Make sure trim region is specificied by procedure or header keyword.
;This part depends on how many amplifiers were used for the readout
  bad = ''
  ver = sxpar(header, 'E_HVERS', count=count)
  if(ver gt 1.001) then begin ; Earlier versions did not support mutiple amplifiers
    n_amp = sxpar(header, 'E_AMPL', count=count)
  endif else begin
    n_amp = 1
  endelse

  if(n_amp gt 1) then begin   ; More than one amplifier
    xlo = sxpar(header, 'E_XLO*', count=count)
    if count eq 0 then bad = 'XLO'
    xhi = sxpar(header, 'E_XHI*', count=count)
    if count eq 0 then bad = 'XHI'
    ylo = sxpar(header, 'E_YLO*', count=count)
    if count eq 0 then bad = 'YLO'
    yhi = sxpar(header, 'E_YHI*', count=count)
    if count eq 0 then bad = 'YHI'

;Check we have encountered inconsistent ranges
    if keyword_set(bad) then begin
      message, /info, bad + ' not specified by argument or in header'
      return, 'Unable to trim image'
    endif

;Make sure trim region is a subset of actual image.
    sz = size(image)
    i1 = where(xlo lt 0 or xlo ge sz[1], n1)
    i2 = where(ylo lt 0 or ylo ge sz[2], n2)
    i3 = where(xhi lt 0 or xhi ge sz[1], n3)
    i4 = where(yhi lt 0 or yhi ge sz[2], n4)
    i5 = where(xlo ge xhi, n5)
    i6 = where(ylo ge yhi, n6)

    err = ''
    if n1 gt 0 then err = 'Error specifying X trim region:' $
         + ' XLO';+strtrim(xlo,2)
    if n2 gt 0 then err = 'Error specifying Y trim region:' $
         + ' YLO';+strtrim(ylo,2)
    if n3 gt 0 then err = 'Error specifying X trim region:' $
         + ' XHI';+strtrim(xhi,2)
    if n4 gt 0 then err = 'Error specifying Y trim region:' $
         + ' YHI';+strtrim(yhi,2)
    if n5 gt 0 then err = 'Error specifying X region boundaries:' $
         + ' XLO>=XHI';+strtrim(xlo[i5],2)+'>='+strtrim(xhi[i5],2)
    if n6 gt 0 then err = 'Error specifying Y region boundaries:' $
         + ' YLO>=YHI';+strtrim(ylo[i6],2)+'>='+strtrim(yhi[i6],2)
    if keyword_set(err) then begin
      message, /info, err
      return, 'Unable to trim image'
    endif

    linear = sxpar(header, 'E_LINEAR', count=count)
    if(count gt 0 and not linear) then begin
      pref = sxpar(header, 'E_PREFMO', count=count)
      image = call_function('nonlinear_' + pref, image, header)
      i = where(strpos(header, 'E_LINEAR') ge 0, ni)
      if(ni gt 0) then header = [header[0:i-1], header[i+1:*]]
      sxaddpar, header, 'E_LINEAR', 'T'  $
              , ' Image corrected of non-linearity'
      ii = where(strpos(header, 'E_GAIN*') ge 0, ni)
      if(ni gt 0) then begin
        for i=0,ni-1 do begin
          k=ii[i]
          header = [header[0:k-1], header[k+1:*]]
        endfor
      endif          
      sxaddpar, header, 'E_GAIN', 1  $
              , ' Image was converted to e-'
    endif

;Trim image to leave only the subimage containing valid image data.
;For two amplifiers we assume a single vertical or horizontal gap.
;With four amplifiers we can have a cross.

    if(n_amp eq 2) then begin
      if(xlo[0] eq xlo[1]) then begin
        xsize = xhi[0] - xlo[0] + 1
        ysize = yhi[0] - ylo[0] + 1 $
              + yhi[1] - ylo[1] + 1
        timage = make_array(xsize, ysize, TYPE = sz[3], /NOZERO)
        ysize = yhi[0] - ylo[0] + 1
        timage[0:xsize-1, 0:ysize-1] $
             = image[xlo[0]:xhi[0] $ ;fill trimmed image
                   , ylo[0]:yhi[0]]
        timage[0:xsize-1, ysize:ysize+yhi[1]-ylo[1]] $
             = image[xlo[1]:xhi[1] $ ;fill trimmed image
                   , ylo[1]:yhi[1]]
      endif else if(ylo[0] eq ylo[1]) then begin
        xsize = xhi[0] - xlo[0] + 1 $
              + xhi[1] - xlo[1] + 1
        ysize = yhi[0] - ylo[0] + 1
        timage = make_array(xsize, ysize, TYPE = sz[3], /NOZERO)
        xsize = xhi[0] - xlo[0] + 1
        timage[0:xsize-1, 0:ysize-1] $
             = image[xlo[0]:xhi[0] $ ;fill trimmed image
                   , ylo[0]:yhi[0]]
        timage[xsize:xsize+xhi[1]-xlo[1],0:ysize-1] $
             = image[xlo[1]:xhi[1] $ ;fill trimmed image
                   , ylo[1]:yhi[1]]
      endif else begin
        message,'The two CCD sections are aligned neither in X nor in Y'
      endelse
    endif else if(n_amp eq 4) then begin
      message,'4-amplifier section is not implemented yet'
    endif
  endif else begin
    if keyword_set(xr) then begin
      xlo = xr(0)
      xhi = xr(1)
    endif else begin
      xlo = sxpar(header, 'E_XLO', count=count)
      if count eq 0 then bad = 'XLO'
      xhi = sxpar(header, 'E_XHI', count=count)
      if count eq 0 then bad = 'XHI'
    endelse
    if keyword_set(yr) then begin
      ylo = yr(0)
      yhi = yr(1)
    endif else begin
      ylo = sxpar(header, 'E_YLO', count=count)
      if count eq 0 then bad = 'YLO'
      yhi = sxpar(header, 'E_YHI', count=count)
      if count eq 0 then bad = 'YHI'
    endelse

;Check we have encountered inconsistent ranges
    if keyword_set(bad) then begin
      message, /info, bad + ' not specified by argument or in header'
      return, 'Unable to trim image'
    endif

;Make sure trim region is a subset of actual image.
    sz = size(image)
    err = ''
    if xhi lt xlo then err = 'Error specifying X trim region:' $
         + ' XHI<XLO (' + strtrim(xhi, 2) + '<' + strtrim(xlo, 2) + ')'
    if yhi lt ylo then err = 'Error specifying Y trim region:' $
         + ' YHI<YLO (' + strtrim(yhi, 2) + '<' + strtrim(ylo, 2) + ')'
    if xlo lt 0 or xhi ge sz(1) then err = 'X trim region [' $
         + strtrim(xlo, 2) + ',' + strtrim(xhi, 2) + '] not contained' $
         + ' in image (0<X<' + strtrim(sz(1)-1, 2) + ')'
    if ylo lt 0 or yhi ge sz(2) then err = 'Y trim region [' $
         + strtrim(ylo, 2) + ',' + strtrim(yhi, 2) + '] not contained' $
         + ' in image (0<Y<' + strtrim(sz(2)-1, 2) + ')'
    if keyword_set(err) then begin
      message, /info, err
      return, 'Unable to trim image'
    endif

;Trim image to leave only the subimage containing valid image data.
    timage = image(xlo:xhi, ylo:yhi)      ;trimmed image
  endelse

;Flip image (if necessary) to achieve standard image orientation.
  if n_elements(orient) eq 0 then begin
    orient = sxpar(header, 'E_ORIENT', count=count)
    if count eq 0 then begin
      message, /info, 'ORIENT not specified by argument or in header'
      return, 'Unable to reorient image'
    endif
  endif
  return, rotate(temporary(timage), orient)

end

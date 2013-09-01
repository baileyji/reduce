pro getspec, im, orc, xwd, spec, sunc, gain, dark ,readn  $
           , NOOPT=noopt, COLRANGE=colrange, SF_SMOOTH=sf_smooth $
           , SP_SMOOTH=sp_smooth, OSAMPLE=osample, SWATH_WIDTH=swath_width $
           , MASK=mask, PLOT=iplot, ORDER_RANGE=order_range $
           , TELLURIC=telluric, FILENAME=filename $
           , SLIT_TILT_FILE=slit_tilt_file $
           , SWATH_BOUNDARIES=swath_boundaries $
           , WING_SMOOTH_FACTOR=wing_smooth_factor

;Subroutine extracts spectra from all orders described in orc.
; im (input array (# columns , # rows)) image from which orc and back were
;   determined and from which spectrum is to be extracted.
; orc (input array (# coeffs , # orders)) polynomial coefficients (from FORDS)
;   that describe the location of complete orders on the image.
; xwd (input array (2, # orders)) fractional extraction widths (from GETXWD)
;   xwd(0,*) is the useful fraction of the order below the central line
;   and xwd(1,*) is the corresponding fraction above
; spec (output array (# columns , # orders)) final extracted, background
;   subtracted spectrum from im.
; sunc (output array (# columns , # orders)) final extracted, background
;   subtracted uncertainty spectrum from im.  If not using optimal extraction,
;   this array is set to unity everywhere.
;24-Oct-89 JAV	Create.
;01-Nov-89 GBB	Modified to allow no background subtraction.
;10-Nov-89 JAV  Cleaned up background subtraction logic.
;03-Dec-89 JAV	Added fractional extraction width to argument list.
;14-Dec-89 JAV	Fixed checks for swath off edge of spectrum.
;19-Jan-89 JAV	Fixed coefficient calculation in 'Arc Off Edge of Image' tests.
;23-Jan-89 JAV	Really fixed 'Arc Off Edge if Image' tests.
;06-Jun-90 JAV	Added argument to GETARC call so total counts are returned.
;04-Sep-90 JAV	Fixed background subtraction bug; backgd/pixel is subtracted
;		 from spectrum counts/pixel BEFORE conversion to total counts.
;13-Sep-90 JAV	Added user specified extraction width logic ($hamxwd stuff).
;18-Apr-92 JAV	Updated global variable list/interpretations. Changed xwd
;		 logic.
;05-May-98 CMJ  Added sunc (the uncertainty spectrum) to routine.
;16-Jul-98 CMJ  In implementation of slit function, have changed the usage
;                so that the central 60% of an order around the current center
;                is used to define the slit function, instead of the full
;                100%.  Ideally, the full 100% should be used, but due to
;                spike-like features in the flatfield (which come from the
;                sharp fall in intensity at the top and bottom of the slit),
;                the central regions of the troughs in a flatfielded image
;                are messed up and this can dominate the slit function in
;                weakly exposed spectra.
;17-Jul-98 CMJ  Changed the immediate above up to 66% as the flat-field
;                problem seems to be much better, though it could still
;                affect very weak spectra.
;20-Jul-98 CMJ  Added the keyword noopt to override optimal extraction.
;                This allows control of this parameter from hamspec and
;                helps facilitate the reduction of long slit day sky and
;                wide flat spectra.
;13-Aug-98 CMJ  Added special handling of order 38 in id=34 spectra.
;06-Jun-99 JAV  Added trace diagnostic message every 5 orders.
;26-Jan-00  NP  Removed common ham_aux, replaced with data from
;               inst_setup structure available in ham.common
;01-May-04  NP  Replaced single value of xwd with an array in order
;               to treat properly spactra with complex image slicers and
;               polarization analyzer.
;25-Feb-08  NP  Return swath boundaries if requested by the caller
;
;BUG: IDL will halt with an error if image is indexed with an invalid index.
;     ANA used to return the last element in the array. New logic is needed
;     here to handle orders that are partially off the chip.


if n_params() lt 4 then begin
  print,'syntax: getspec,im,orc,xwd,spec[,sunc...]'
  retall
end

  ham_optimal = 1

;Define useful quantities.
  ncol = n_elements(im(*,0))				;# columns in image
  nrow = n_elements(im(0,*))				;# rows in image

;Default order range span the whole image
  if(not keyword_set(colrange)) then begin
    colrange=intarr(2,nord)
    colrange(1,*)=ncol-1
  endif

  if((size(orc))(0) gt 1) then begin
    ncoef = n_elements(orc(*,0))			;# polyn. coeffs
    nord =  n_elements(orc(0,*))			;# full orders in orc
  endif else begin
    ncoef = n_elements(orc)
    nord = 1
  endelse
  ix = findgen(ncol)					;column indicies
  spec = fltarr(ncol,nord)				;init spectrum
  if keyword_set(noopt) then begin
     sunc = fltarr(ncol,nord)+1.                        ;thar unc = 1
  endif else begin
     sunc = fltarr(ncol,nord)				;init uncertainties
  endelse
  orcend = fltarr(ncoef,nord+2)				;init extended orcs

;GETARC needs order location coefficients (orc) on both sides of arc swath to
;  be extracted. In order to extract the lowest and highest orders specified
;  in orc, we need to extend orc one extra order on each end. We shall do so
;  by linearly extrapolating the last two orc on each end.
;Extend orc on the low end. Check that requested swath lies on image.
  if(nord gt 1) then begin
    ixx=ix(colrange(0,0):colrange(1,0))
    orclo = 2*orc(*,0) - orc(*,1)			;extrapolate orc
    if(xwd(0,0) gt 1.5) then begin                      ;Extraction width in pixels
      coeff=orc(*,0)
      coeff(0)=coeff(0)-xwd(0,0)
    endif else begin                                    ;Fraction extraction width
;     coeff = 0.5 * ((2+xwd)*orc(*,0) - xwd*orc(*,1))
      coeff = 0.5 * ((2+xwd(0,0))*orc(*,0) - xwd(0,0)*orc(*,1))
    endelse
    y = fltarr(ncol)
    y(colrange(0,0):colrange(1,0)) = poly(ixx,coeff)	;low edge of arc
    yoff = where(y lt 0,noff)				;pixels off low edge
  endif else begin
    noff = 0
    orclo = [0.,replicate(0.,ncoef-1)]
  endelse
  if noff gt 0 then begin				;check if on image
;   GETARC will reference im(j) where j<0. These array elements do not exist.
    print,'GETSPEC: Top order off image in columns [' $
      + strtrim(string(yoff(0)),2) + ',' $
      + strtrim(string(yoff(noff-1)),2) + '].'
  endif

;Extend orc on the high end. Check that requested swath lies on image.
  if(nord gt 1) then begin
    ix1 = colrange(0,nord-1)
    ix2 = colrange(1,nord-1)
    ixx=ix(ix1:ix2)
    orchi = 2*orc(*,nord-1) - orc(*,nord-2)             ;extrapolate orc
    if(xwd(1,nord-1) gt 1.5) then begin                 ;Extraction width in pixels
      coeff=orc(*,nord-1)
      coeff(0)=coeff(0)+xwd(1,nord-1)
    endif else begin                                    ;Fraction extraction width
;     coeff = orc(*,nord-1) + xwd*(orc(*,nord-1) - orc(*,nord-2))/2
      coeff = 0.5 * ((2+xwd(1,nord-1))*orc(*,nord-1) - xwd(1,nord-1)*orc(*,nord-2))
    endelse
    y = replicate(nrow-1., ncol)
    y(ix1:ix2) = poly(ixx,coeff)			;high edge of arc
    yoff = where(y gt nrow-1,noff)			;pixels off high edge
  endif else begin
    noff = 0
    orchi = [nrow-1.d0,replicate(0.,ncoef-1)]
  endelse
  if noff gt 0 then begin
;   GETARC will reference im(j) where j > ncol*nrow-1. These array elements do
;     not exist.
    print,'GETSPEC: Bottom order off image in columns [' $
      + strtrim(string(yoff(0)),2) + ',' $
      + strtrim(string(yoff(noff-1)),2) + '].'
  endif

;Define an order set (orcend) extended one extra order on either end.
  for n = 1,nord do orcend(*,n) = orc(*,n-1)
  orcend(*,0) = orclo
  orcend(*,nord+1) = orchi

;;Now loop through orders extracting spectrum and maybe subtracting background.
;  if subback then begin
;    print,'GETSPEC: Extracting spectrum and subtracting background.'
;  end else begin
;    print,'GETSPEC: Extracting spectrum - no background subtraction.'
;  end

  if not keyword_set(noopt) then begin                 ;check for optimal ext.
     if ham_optimal eq 1 then begin
        copt = 'y'
     endif else begin
        copt = 'n'
     endelse
  endif else begin
     copt = 'n'                                        ;set for thoriums
  endelse
  if keyword_set(noopt) then copt = 'n'            ;override if requested

  if copt eq 'n' then begin                        ;no optimal extraction
     if max(xwd) le 1 then begin	           ;fractional orders
       print,'GETSPEC: ' $
         + 'Interpreting extraction width as fraction of each order to mash.'
     end else begin				   ;fixed pixel widths
       print,'GETSPEC: ' $
         + 'Interpreting extraction width as number of pixels to mash.'
     endelse
  endif else begin                                 ;optimal extraction
     print,'GETSPEC: ' $
         + 'Using optimal extraction to produce spectrum.'
  endelse

  if copt eq 'y' then begin                        ;Optimal extraction of the spectrum

    if(keyword_set(slit_tilt_file)) then begin     ;2D slit function is requested
      if(file_test(slit_tilt_file)) then restore,slit_tilt_file
    endif

    if(not keyword_set(sf_smooth)) then sf_smooth=6.

    if(keyword_set(order_range)) then begin
      if(n_elements(order_range) eq 2) then begin
        ofirst = order_range(0)                       ;First order to extract
        olast  = order_range(1)                       ;Last order to extract
        if(ofirst lt 1 or olast gt nord or ofirst gt olast) then begin
          ofirst = 1                                  ;First order to extract
          olast  = nord                               ;Last order to extract
        endif
      endif else begin
       ofirst = 1                                     ;First order to extract
       olast  = nord                                  ;Last order to extract
      endelse
    endif else begin
      ofirst = 1                                      ;First order to extract
      olast  = nord                                   ;Last order to extract
    endelse

    swath_boundaries=intarr(2,1)
    for onum=ofirst,olast do begin                    ;loop thru orders

      ncole = colrange(1,onum-1)-colrange(0,onum-1)+1 ;number of columns to extract
      cole0 = colrange(0,onum-1)                      ;first column to extract
      cole1 = colrange(1,onum-1)                      ;last column to extract

;Background must be subtracted for slit function logic to work but kept
;as part of the FF signal during normalization

      scatter_below  = 0
      yscatter_below = 0
      scatter_above  = 0
      yscatter_above = 0
      
      if nord le 10 then begin
        print,'GETSPEC: extracting relative order ' $
              + strtrim(onum, 2) + ' out of ' + strtrim(nord,2)
      endif else begin
        if (onum-1) mod 5 eq 0 then begin
          print, 'GETSPEC: extracting relative orders ' $
                  + strtrim(onum, 2) + '-' $
                  + strtrim((nord)<(onum+4), 2) + ' out of ' + strtrim(nord,2)
        endif
      endelse

      ycen = poly(ix, orcend(*,onum))			        ;row at order center
;
;ycen=ycen-20 ;Can force offsets from peak (e.g. for binaries)
;
      x_left_lim  = cole0                               ;First column to extract
      x_right_lim = cole1                               ;Last column to extract
      ixx = ix(x_left_lim:x_right_lim)
      ycenn = ycen(x_left_lim:x_right_lim)
      if(xwd(0,onum-1) gt 1.5) then begin               ;Extraction width in pixels
        ymin = ycenn - xwd(0,onum-1)
      endif else begin                                  ;Fractional extraction width
        ymin = ycenn - xwd(0,onum-1) * (ycenn - poly(ixx, orcend(*,onum-1))) ;trough below
      endelse
      ymin = floor(ymin)
      if min(ymin) lt 1 then ymin = ymin - min(ymin) + 1          ;help for orders at edge
      if(xwd(1,onum-1) gt 1.5) then begin               ;Extraction width in pixels
        ymax = ycenn + xwd(1,onum-1)
      endif else begin                                  ;Fractional extraction width
        ymax = ycenn + xwd(1,onum-1) * (poly(ixx, orcend(*,onum+1)) - ycenn) ;trough above
      endelse
      ymax = ceil(ymax)
      if max(ymax) ge nrow-1 then ymax = ymax - max(ymax) + nrow - 2 ;helps at edge

;Define a fixed height area containing one spectral order
      y_lower_lim = fix(min(ycen(cole0:cole1)-ymin))              ;Pixels below center line
      y_upper_lim = fix(min(ymax-ycen(cole0:cole1)))              ;Pixels above center line
      if(y_lower_lim eq 0 or y_upper_lim eq 0) then stop

;
;if onum lt nord then begin ;Can force use of previous slitf at edge of CCD
;
      unc=replicate(-2.,ncol)
      mkslitf, im, scatter_below, yscatter_below $
             , scatter_above, yscatter_above $
             , ycen, y_lower_lim, y_upper_lim $
             , yslitf, slitf, binc, onum $
             , x_left_lim=cole0, x_right_lim=cole1, PLOT=iplot $
             , LAMBDA_SF=sf_smooth, LAMBDA_SP=lambda_sp, blz=blz $
             , OSAMPLE=osample, SWATH_WIDTH=swath_width, MASK=mask $
             , GAIN=gain, READN=readn, /NO_SCATTER $
             , TELLURIC=telluric, FILENAME=filename $
             , SLIT_TILT=shear_x $
             , SWATH_BOUNDARIES=sw_b $
             , WING_SMOOTH_FACTOR=wing_smooth_factor $
             , UNCERTAINTY=unc

;Handle binary in slit. Determine yslitf range explicitly for each binary.
;iwhr=where(yslitf lt -7.0)        ;EQ Peg B (-7.0 is an example)
;yslitf=yslitf(iwhr)
;slitf_trim = slitf
;oldnorm = median(total(slitf, 1))
;slitf = slitf(iwhr,*)
;for j=0, nord-1 do slitf(*,j) = oldnorm * slitf(*,j) / total(slitf(*,j))
;

;
;endif ;End of logic to force use of previous slitf at edge of CCD
;

;      hamopt, im, bg, ycen, yslitf $
;            , slitf, binc, s, u, gain, dark, readn $
;            , x_left_lim=x_left_lim,x_right_lim=x_right_lim  ;optimally extract
;      spec(x_left_lim:x_right_lim,onum-1) = s		 ;save in output array
;      sunc(x_left_lim:x_right_lim,onum-1) = u		 ;save uncertainties
      spec(x_left_lim:x_right_lim,onum-1) = blz(x_left_lim:x_right_lim)	 ;save in output array
      swath_boundaries=[[swath_boundaries],[sw_b]]
      if(max(unc) lt -1) then begin
        sunc(x_left_lim:x_right_lim,onum-1) = sqrt(abs(blz(x_left_lim:x_right_lim)))
      endif else begin
        sunc(x_left_lim:x_right_lim,onum-1) = unc(x_left_lim:x_right_lim)
      endelse
    endfor
    swath_boundaries=swath_boundaries(*,1:*)

  endif else begin					;else: getarc style
    for onum=1,nord do begin				;loop thru orders
      x_left_lim  = colrange(0,onum-1)                  ;First column to extract
      x_right_lim = colrange(1,onum-1)                  ;Last column to extract
      awid = xwd(0,onum-1) + xwd(1,onum-1)
      getarc,im,orcend,onum,awid,arc,pix $
      	    ,x_left_lim=x_left_lim,x_right_lim=x_right_lim	;extract counts/pixel
;      if subback then arc = arc - back(*,onum-1)	        ;subtract backgd/pixel
      spec(x_left_lim:x_right_lim,onum-1) = arc * pix		;store total counts
      sunc(x_left_lim:x_right_lim,onum-1) = sqrt(abs(arc*pix*gain + dark + $
                            pix*readn^2)) / gain          ;estimate uncertainty
    endfor
  endelse

return
end

pro hamspec, im ,head, orc, dxwd, dsxwd, obase, spec, ThAr=thar $
           , LONG=long, SIG=sunc, COLRANGE=colrange, SF_SMOOTH=sf_smooth $
           , SP_SMOOTH=sp_smooth, OSAMPLE=osample, SWATH_WIDTH=swath_width $
           , PLOT=iplot, MASK=mask , ORDER_RANGE=order_range $
           , TELLURIC=telluric, FILENAME=filename $
           , SWATH_BOUNDARIES=swath_boundaries $
           , WING_SMOOTH_FACTOR=wing_smooth_factor $
           , SLIT_TILT_FILE=slit_tilt_file

;HAMSPEC extracts spectrum and background from each order
; [spec (optional output array (# columns , # full orders))] extracted and
;   background subtracted spectrum.
; [bkgc (optional output array (# coeffs , # full orders))] coefficients of fit
;   to background values subtracted from spectrum to give spec.
;23-Oct-89 JAV	Create.
;18-Apr-92 JAV	Updated global variable list/interpretations. Reenabled maskim
;		 when thar=1.
;29-Apr-92 JAV	Now clip baseline column of FITS images.
;19-Sep-92 JAV	obase read from .fmt file mandatory. Write obase to .isp file.
;05-Jul-94 CMJ  Moved backgound determination code to routine GETBKG in order
;                 to facilitate interchanging Jeff's method with Gibor's by
;                 changing what routine to use.  GETBKG is Jeff's trough
;                 method.  BCKGRC is Gibor's method using contf to do a more
;                 global fit.  I prefer GETBKG and use of the 2.5arcsec decker
;                 correction for Hamilton images.
;06-Dec-94 CMJ  RDFITS call now has ClipBase set to 32 for McD Cass Echelle.
;15-Dec-94 JAV	Value of clipbase now depends on ham_id. Currently 0, 1, or 32
;		 for KPNO, Lick/Keck, and McDonald, respectively. Added logic
;		 to handle KPNO specific requirements (clipbase, transpose).
;                Add bias frame subtraction for KPNO (ham_id=40). Extracted
;		 image preprocessing to getimage.pro.  Note: (CMJ) for
;                McDonald data; the 2.1m requires a ClipBase of 32 while the
;                AIS chip on the 2.7m requires ClipBase=0.
;06-Feb-95 CMJ  Before saving extracted spectrum, test for McDonald dewars
;                and call cformat.pro if present.  Both dewars need to have
;                the orders reversed and dewar 30 to have the ordering of the
;                orders reversed.
;21-Nov-97 JAV	Compute fraction of pixels between orders containing background
;		 and pass that information to gettrof via getbkg.
;05-May-98 CMJ  Added sunc (uncertainty spectrum) to getspec call.
;06-May-98 CMJ  Put a version number comment into writing of isp files.
;                 Current version number is 3.0.
;16-Jul-98 CMJ  Changed bgmul to 2.2 (old value was 1.4).  While working on
;                 McDonald high res data, the background was going too far
;                 into the wings of the real spectra and artificially raising
;                 the background level.  We may want to consider having two
;                 values of bgmul for high and low res data.
;17-Jul-98 CMJ  Made the size of bgmul depend on ham_id.  It =2.5 for id=34,
;                 and =2.0 for other id's.  Also changed logical handling of
;                 ham_dord so that the value =2 can be properly dealt with.
;                 This choice uses the default orders but tries to determine
;                 a global shift for them.
;20-Jul-98 CMJ  Added keyword long, for long slit spectra (typically wide
;                 flats and day skys from McDonald).  This keyword forces
;                 the use of the default order locations and the very narrow
;                 use of trough pixels to estimate backgrounds.
;17-Aug-98 CMJ  Put in a call to getxwd.pro to determine the extraction width
;                 for each spectrum individually.
;01-Sep-98 CMJ  Given improved background logic, have changed bgmul back
;                 to 1.7 to get best background estimate.
;02-Jun-99 JAV  Call cformat for KPNO feed, Camera 6, F3KB to reorder spectra.
;		 Pass biasfile='' for ham_id which do not have bias images.
;07-Jun-99 JAV  Allow for optimal extraction of longslit spectra.
;26-Jan-00  NP  Removed common ham_aux, replaced with data from
;                inst_setup structure available in ham.common
;25-Feb-08  NP  Return swath boundaries if requested by the caller

  if n_params() lt 7 then begin
    print,'syntax: hamspec,im,head,orc,dxwd,dsxwd,obase,spec[,MASK=mask $'
    print,'               [,SIG=sig[,/THAR[,COLRANGE=colrange $'
    print,'               [,SF_SMOOTH=sf_smooth[,SP_SMOOTH=sp_smooth $'
    print,'               [,OSAMPLE=osample[,SWATH_WIDTH=swath_width $'
    print,'               [,/PLOT[,ORDER_RANGE=order_range $'
    print,'               [,TELLURIC=telluric[,FILENAME=filename $
    print,'               [,SWATH_BOUNDARIES=swath_boundaries $
    print,'               [,SLIT_TILT_FILE=slit_tilt_file]]]]]]]]]]]]]]'
    retall
  end

;Extract frame information from the header
  readn = sxpar(head, 'E_READN')
  dark  = sxpar(head, 'E_BACKG')
  gain  = sxpar(head, 'E_GAIN')

;For long slit, heavily smooth to improve location of order centroids.
  if keyword_set(long) then begin
    im = median(im, 3)                 ;fix bad pixels
    for ismit=0, nsmit-1 do im = smooth(temporary(im), 7)
  endif

;Number of orders and columns
  nord = 1
  if((size(orc))(0) gt 1) then nord = (size(orc))(2);number of orders
  ncol = (size(im))(1)			            ;number of columns

;Determine fractional extraction width.
  if not keyword_set(thar) then begin            ;if we are not going thar
    if(keyword_set(dxwd)) then begin
      if(n_elements(dxwd) eq 1)           then xwd = replicate(0.5*dxwd,2,nord) $
      else if(n_elements(dxwd) eq 2*nord) then xwd = dxwd $
      else begin
        print,'Inconsistent size of dxwd ',size(dxwd),'. Should be 1 or [2,nord]'
        stop
      endelse
      print, 'HAMSPEC: Using default extraction width(s)'
    endif else begin
      print,'HAMSPEC: Using GETXWD to find extraction width.'
      getxwd,im,orc,xwd,sxwd,colrange=colrange,/GAUSS     ;find xwd
    endelse
  endif else begin
    if(keyword_set(dxwd)) then begin
      if(n_elements(dxwd) eq 1)           then xwd = replicate(0.5*dxwd,2,nord) $
      else if(n_elements(dxwd) eq 2*nord) then xwd = dxwd $
      else begin
        print,'Inconsistent size of dxwd ',size(dxwd),'. Should be 1 or [2,nord]'
        stop
      endelse
      print, 'HAMSPEC: Using default extraction width '+strtrim(mean(xwd),2)
    endif
  endelse

;Default order range span the whole image
  if(not keyword_set(colrange)) then begin
    colrange=intarr(2,nord)
    colrange(1,*)=ncol-1
  endif

;Extract spectrum.
  if keyword_set(thar) then begin		;true: extracting thar lamp
    getspec,im,orc,xwd,spec,0,gain,dark,readn $
           ,colrange=colrange, /NOOPT           ;extract w/o bg subtraction
    sunc=spec*0.                                ;set uncertainty spec. for thar
  endif else begin				;else: extracting spectrum
    getspec, im, orc, xwd, spec, sunc $         ;extract spectrum
           , gain, dark, readn $
           , colrange=colrange, sf_smooth=sf_smooth, sp_smooth=sp_smooth $
           , osample=osample, SWATH_WIDTH=swath_width $
           , MASK=mask, PLOT=iplot, ORDER_RANGE=order_range $
           , TELLURIC=telluric, FILENAME=filename $
           , SLIT_TILT_FILE=slit_tilt_file $
           , WING_SMOOTH_FACTOR=wing_smooth_factor $
           , SWATH_BOUNDARIES=swath_boundaries
  endelse

return
end

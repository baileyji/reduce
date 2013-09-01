pro hamflat, im, head, orders, blzcof, COLRANGE=colrange $
           , OSAMPLE=osample, FXWD=fxwd, SF_SMOOTH=sf_smooth $
           , SP_SMOOTH=sp_smooth, SWATH_WIDTH=swath_width, MASK=mask $
           , POLY=pol, PLOT=iplot, ORDER_HEIGHT=order_height $
           , NOSCATTER=noscatter, THRESHOLD=thres $
           , SWATH_BOUNDARIES=swath_boundaries $
           , WING_SMOOTH_FACTOR=wing_smooth_factor $
           , UNCERTAINTY=bunc $
           , POLARIZATION=polariz

;HAMFLAT reads in a master flat, determines scattered light inside inter-order
;  gaps, constructs models for each sp. order and each swath using
;  slit decomposition algorithm and replaces the original image with the
;  ratio of the original to the model where the model has signal above
;  threshold and 1 elsewhere. Sacttered light is subtracted when constructing
;  the model and added back to the model when computing the ratio. Polarization
;  requires special treatment of the scattered light. What flats are exposed
;  through a polarizing beam-splitter space between two polarizations of the
;  same order should not be used for estimating scattered light. HAMFLAT assumes
;  that the oprder location array points at a complete set of pairs so that for
;  each spectral order we have two polarizations. The scattered light is then
;  estimated between spectral orders and not between polarizations. 
;  As a by-product HAMFLAT returns the blaze functions for each order.
;
;Inputs:
; setfile (string) root of input filenames containing information common to
;   all observations made with a particular spectrograph setting. A complete
;   filename will be constructed by appending a standardized file extension.
;   The following file is expected to exist:
;   * setfile.ord - default order location coefficients (from HAMDORD)
;   The following files will be created by HAMFLAT:
;   * setfile.flt - normalized flat field image (assoc. var.).
;   * setfile.blz - blaze function polynomial fit coefficients (WDSKed).
; obsfile (string) root of input/output filenames containing information
;   specific to a given observation. A complete filename will be constructed
;   by appending a standardized file extension.
;   One of the following files is expected to exist (depending on ham_image):
;   * obsfile.fits - image in fits format (ham_image=0)
;   * obsfile.dsk - image in wdsk/rdsk format (ham_image=1)
;   The following file may be created by HAMFLAT:
;   * obsfile.ord - order location coefficients
;Outputs:
; im (array (ncol,nrow)) normalized version of the image in input OBSFILE.
;   Copy of im is written in setfile.flt.
;Notes:
; Previous version of hamflat, which fit along individual arcs, has been
;  renamed hamflat_arc.pro. The routines FITARCS, EXAFC, SMENDS, QUADEX,
;  NORMFLAT, and COLFIT are no longer used. MKSLITF is now used here, as
;  well as for optimal extraction of spectra.
;20-Oct-89 JAV  Create.
;18-Apr-92 JAV  Updated global variable list/interpretations. Fixed fits format
;        logic error. Implemented trace level. Fixed default order
;        location logic. Changed file extensions to lower case.
;29-Apr-92 JAV  Now clip baseline column of FITS images. Moved normalized flat
;        bad pixel logic from flatim.pro to here.
;06-Dec-94 CMJ  RDFITS call now has ClipBase set to 32 for McD Cass Echelle
;15-Dec=94 JAV  Value of clipbase now depends on ham_id. Currently 0, 1, or 32
;                for KPNO, Lick/Keck, and McDonald, respectively. Added logic
;                to handle KPNO specific requirements (clipbase, transpose).
;        Add bias frame subtraction for KPNO (ham_id=40). Extracted
;        image preprocessing to getimage.pro.
;03-Jun-99 JAV  Always pass biasfile='', even for KPNO flats. Modify offset
;        and trim parameters for KPNO flats, which are assumed to be
;        pretrimmed.
;06-Jun-99 JAV  Changed normalization algorithm from fits along individual
;        arcs to division by the slit function. Most of the changes
;        involve calling different routines (sfnorm, rather than
;        exorc, exafc, fitarcs, and normflat).
;08-Feb-2008 NP, added explicit parameter for the minimum flat signal to be
;            in normalization.
;25-Feb-2008 NP, return swath boundaries if requested.
;04-Mar-2008 NP, return uncertainties in the blaze functions.
;29-Jul-2009 NP, add handling of flats taken with polarizing beam splitter. 
;
  if n_params() lt 3 then begin
    print,'syntax: hamflat,im,head,orders,blzcoef[,COLRANGE=colrange[,/PLOT]]'
    retall
  end

;Force pixels to be positive to prevent negative or infinite flattened image
;  pixels.
  badp=where(im le 0,nbadp)             ;find pixels le 0
  if nbadp gt 0 then begin              ;true: fix pixels
    im(badp)=1.0                    ;and fix them
    print,'HAMFLAT:' + strcompress(nbadp) $
      + ' zero or negative pixels set to one.'
  endif

;Extract frame information from the header
  readn = sxpar(head, 'E_READN')
  dark  = sxpar(head, 'E_BACKG')
  gain  = sxpar(head, 'E_GAIN')
  
;Default minimum signal for the flat to be used
  if(keyword_set(thres)) then threshold=thres>10000L else threshold=90000L

;Check if the flat has reasonable amount of pixels with signal above the threshold
  ii=where(im gt threshold, nii)
  nn=double(n_elements(im))
  if(nii*10L lt nn) then begin
    print,'HAMFLAT: Your flat has only '+strtrim(round(nii/nn*100.),2) $
         +'% pixels with signal above '+strtrim(threshold,2)
    print,'HAMFLAT: Do you want to continue? [Y/n]'
    answer=get_kbrd(1)
    if(answer eq 'n' or answer eq 'N') then begin
      im[*,*]=1.
      return
    endif
  endif

;Use slit functions to normalize flat (and fetch blaze function).
  sfnorm, im, orders, blzcof=blzcof, dark=dark $
        , readn=readn, gain=gain, COLRANGE=colrange $
        , OSAMPLE=osample, FXWD=fxwd, SF_SMOOTH=sf_smooth $
        , SP_SMOOTH=sp_smooth, SWATH_WIDTH=swath_width, MASK=mask $
        , PLOT=iplot, POLY=pol, ORDER_HEIGHT=order_height $
        , NOSCATTER=noscatter $
        , THRESHOLD=threshold $
        , SWATH_BOUNDARIES=swath_boundaries $
        , WING_SMOOTH_FACTOR=wing_smooth_factor $
        , UNCERTAINTY=bunc $
        , POLARIZATION=polariz

;Force positive pixels in normalized flat.
  iwhr = where(im le 0, nwhr)				;look for bad pixels
  if nwhr gt 0 then begin				;true: bad pixels exist
    print,'HAMFLAT: ' + strtrim(string(nwhr),2) $
      + ' bad pixels in normalized flat set to unity.'
    im(iwhr) = 1.0					;adjust bad pixels
  endif

  print,'HAMFLAT: Flat field image processed - returning to caller.'
  return
end

pro sfnorm, im, orc, blzcof=blazecurve, dark=dark, readn=readn $
          , gain=gain, COLRANGE=colrange, OSAMPLE=osample $
          , SF_SMOOTH=sf_smooth, SP_SMOOTH=sp_smooth, SWATH_WIDTH=swath_width $
          , MASK=mask, POLY=pol, PLOT=iplot, FXWD=fxwd $
          , ORDER_HEIGHT=order_height, NOSCATTER=noscatter $
          , THRESHOLD=thres $
          , SWATH_BOUNDARIES=swath_boundaries $
          , WING_SMOOTH_FACTOR=wing_smooth_factor $
          , UNCERTAINTY=bunc $
          , POLARIZATION=polariz

;Use slit functions to normalize an echelle image of a flat field lamp.
;Inputs:
; im (array(ncol,nrow)) image from which orc and back were determined and
;   from which spectrum is to be extracted.
; orc (array(orcdeg,nord)) polynomial coefficients (from FORDS) that describe
;   the location of complete orders on the image.
; dxw float scalar, fractional width of the order to be normalized
;
;Outputs:
; blzcof (array(blzdeg,nord)) coefficients of (perhaps broken) polynomial
;   fit to extracted spectrum of flat lamp.
;Notes:
; Various diagnostic plots can be enabled by changing "if 0" to "if 1" below.
;History:
;05-Jun-1999 JAV, Adapted from getspec.pro
;26-Jan-2000 NP, removed common ham_aux, replaced with data from
;           inst_setup structure available in ham.common
;09-Apr-2001 NP, added parameters COLRANGE to handle partial orders,
;            OSAMPLE to control oversampling in MKSLITF, SWATH_WIDTH -
;            swath width used to determine the slit function in MKSLITF
;            SF_SMOOTH - to control the smoothness of the slit function
;            in MKSLITF. Also added the logic to handle MASK of bad pixels
;            The mask is supposed to have on ly two values: 1 for good and
;            0 for bad pixels. The flat field in the masked pixels is set to
;            1.
;08-Feb-2008 NP, added explicit parameter for the minimum flat signal to be
;            in normalization.
;25-Feb-2008 NP, Return swath boundaries if requested by the caller
;04-Mar-2008 NP, return uncertainties in the blaze functions
; 
  if n_params() lt 2 then begin
    print,'syntax: sfnorm,im,orc[,blzcof=blaze[,dark=dark[,readn=readn'
    print,'             [,gain=gain[,COLRANGE=colrange[,/PLOT]]]]]]'
    return
  end

;Internal parameters.
  nloop = 40                              ;number of rejection loops
  bsiz = 0.001                            ;binsize for sigma histogram

;Define useful quantities.
  sz = size(im)
  ncol = sz(1)                            ;number of columns in image
  nrow = sz(2)                            ;number of rows in image
  sz = size(orc)
  ncoef = sz(1)                           ;number of polynomial coefficients
  if(sz(0) gt 1) then nord =  sz(2) $     ;number of full orders in orc
  else                nord = 1

;Default order range span the whole image
  if(not keyword_set(colrange)) then begin
    colrange=intarr(2,nord)
    colrange(1,*)=ncol-1
  endif

;Default normalization width is the whole order 
  if(not keyword_set(fxwd))                 then exwd = replicate(0.5,2,nord) $
  else if(min(fxwd) gt 0) then exwd = fxwd  else exwd = replicate(0.5,2,nord)

;Construct xwd if only one walue is given
  if(n_elements(exwd) eq 1) then begin
    exwd = replicate(exwd*0.5,2,nord)                   ;extraction widths
  endif

;Initializations.
  ix = findgen(ncol)                      ;column indicies
  jx = lindgen(ncol)                      ;column indicies (as long integer)
  orcend = dblarr(ncoef,nord+2)           ;init extended orcs
  blazecurve = dblarr(ncol, nord)         ;optimally filtered blaze functions
  bunc       = dblarr(ncol, nord)         ;uncertainties in the blaze functions

;Getarc needs order location coefficients (orc) on both sides of arc swath to
;  be extracted. In order to extract the lowest and highest orders specified
;  in orc, we need to extend orc one extra order on each end. We shall do so
;  by linearly extrapolating the last two orc on each end.
;Extend orc on the low end. Check that requested swath lies on image.

  noff = 0
  if(nord gt 1) then begin
    orclo = 2*orc(*,0) - orc(*,1)         ;extrapolate orc
    y = poly(ix, 0.5*(orclo+orc(*,0)))    ;low edge of arc
;    yoff = where(y lt 0,noff)             ;pixels off low edge
  endif else begin
    orclo = [2,replicate(0.d0,n_elements(orc)-1)]
  endelse

;Check whether extraction is likely to result in invalid indexing of im.
  if noff gt 0 then begin                 ;check if on image
    print,'GETSPEC: Top order off image in columns [' $
      + strtrim(string(yoff(0)),2) + ',' $
      + strtrim(string(yoff(noff-1)),2) + '].'
  endif

;Extend orc on the high end. Check that requested swath lies on image.

  noff = 0
  if(nord gt 1) then begin
    orchi = 2*orc(*,nord-1) - orc(*,nord-2)  ;extrapolate orc
    y = poly(ix, 0.5*(orchi+orc(*,nord-1)))  ;high edge of arc
;    yoff = where(y gt nrow-1,noff)          ;pixels off high edge
  endif else begin
    orchi = [(size(im))(2)-2,replicate(0.d0,n_elements(orc)-1)]
  endelse

;Check whether extraction is likely to result in invalid indexing of im.
  if noff gt 0 then begin
    print,'GETSPEC: Bottom order off image in columns [' $
      + strtrim(string(yoff(0)),2) + ',' $
      + strtrim(string(yoff(noff-1)),2) + '].'
  endif

;Define an order set (orcend) extended one extra order on either end.
  for i = 1,nord do orcend(*,i) = orc(*,i-1)
  orcend(*,0) = orclo
  orcend(*,nord+1) = orchi

;Initializations.
  x = dindgen(ncol) / (ncol - 1)        ;normalized abscissa
  yprev = replicate(-1, ncol)           ;last column normalized

;Loop through the orders, fitting a slit function, optimally extracting
; the spectrum, fitting the extracted spectrum, and then using the fit
; and slit function to normalize the original flat. Regions not covered
; by the slit function are set to exactly 1.0.

  if(not keyword_set(sf_smooth)) then sf_smooth=8.  ;default slit function smoothness

; Fit the scattered light. The approximation is returned in 2D array bg for each
; inter-order troff

;Status report.
  print,'SFNORM: using mkscatter to evaluate the model for background.'
  if(keyword_set(pol)) then begin
    print,'SFNORM: using '+strtirm(pol,2 )+'-order polynomial to fit ' $
        + 'the background.'
  endif
  if(keyword_set(noscatter)) then begin
    bg=dblarr(ncol,nord+1)
    ybg=fltarr(ncol,nord+1)
    xx=dindgen(ncol)
    for onum=1,nord-1 do begin
      ybg(*,onum)=(poly(xx,orc(*,onum-1))+poly(xx,orc(*,onum)))*0.5
    endfor
    ybg(*,   0)=2.*poly(xx,orc(*,     0))-ybg(*,     1)
    ybg(*,nord)=2.*poly(xx,orc(*,nord-1))-ybg(*,nord-1)
  endif else begin
    mkscatter, im, orc, bg, ybg, colrange=colrange, LAMBDA_SF=sf_smooth $
             , SWATH_WIDTH=swath_width, OSAMPLE=osample, MASK=mask $
             , GAIN=gain, READN=readn, DATA=bg_data, POLY=pol $
             , ORDER_HEIGHT=order_height, LAMBDA_sp=400.
  endelse
  
;Status report.
  print,'SFNORM: using slit functions to normalize flat.'
  im_norm = replicate(1., ncol, nrow);              ;prepare a normalized array 
  im_ordr = fltarr(ncol, nrow);                     ;order shape array 

  ofirst = 1                                        ;starting order
  for onum=ofirst,nord do begin                     ;loop thru orders
    ncole = colrange(1,onum-1)-colrange(0,onum-1)+1 ;number of columns to extract
    cole0 = colrange(0,onum-1)                      ;first column to extract
    cole1 = colrange(1,onum-1)                      ;last column to extract

;Background must be subtracted for slit function logic to work but kept
;as part of the FF signal during normalization

    if(keyword_set(polariz)) then begin
      oo = ((onum-1)/2)*2+1        ;skip inter-polarization gaps
       scatter_below =  bg(*,oo-1)
      yscatter_below = ybg(*,oo-1)
       scatter_above =  bg(*,oo+1)
      yscatter_above = ybg(*,oo+1)
    endif else begin
       scatter_below =  bg(*,onum-1)
      yscatter_below = ybg(*,onum-1)
       scatter_above =  bg(*,onum)
      yscatter_above = ybg(*,onum)
    endelse

    if nord le 10 then begin
      print, 'SFNORM: processing relative order ' $
              + strtrim(onum, 2) + '  of ' + strtrim(nord,2)
    endif else begin
      if (onum-1) mod 5 eq 0 then begin
        if(onum lt nord) then $
          print, 'SFNORM: processing relative orders ' $
                + strtrim(onum, 2) + '-' $
                + strtrim((nord)<(onum+4), 2) + '  of ' + strtrim(nord,2) $
        else $
          print, 'SFNORM: processing relative order ' $
                + strtrim(onum, 2) + ' (the last one)'
      endif
    endelse

;Use order fits to define range of pixels to consider.
    ixx = ix(cole0:cole1)
    ycen = poly(ix, orcend(*,onum))                                ;row at order center
    if(exwd(0,onum-1) gt 1.) then  ymin=ycen(cole0:cole1)-exwd(0,onum-1) $;trough below
    else ymin = ycen(cole0:cole1)*(1.-exwd(0,onum-1)) $
              + poly(ixx, orcend(*,onum-1))*exwd(0,onum-1)
    ymin = ceil(ymin)
    if min(ymin) lt 1 then ymin = ymin - min(ymin) + 1             ;helps the bottom order
    if(exwd(1,onum-1) gt 1.) then  ymax=ycen(cole0:cole1)+exwd(1,onum-1) $;trough above
    else ymax = ycen(cole0:cole1)*(1.-exwd(1,onum-1)) $
              + poly(ixx, orcend(*,onum+1))*exwd(1,onum-1)
    ymax = floor(ymax)
    if max(ymax) ge nrow-1 then ymax = ymax - max(ymax) + nrow - 2 ;helps the top order

;Define a fixed height area containing one spectral order
;    if(onum eq ofirst) then $
;      y_lower_lim = fix(min(ycen(cole0:cole1)-ymin)) $             ;Pixels below center line
;    else $
;      y_lower_lim = fix(min(ycen(cole0:cole1)-yprev(cole0:cole1))) + 1
    y_lower_lim = ceil(min(ycen(cole0:cole1)-ymin))                ;Pixels below center line
    y_upper_lim = ceil(min(ymax-ycen(cole0:cole1)))+1              ;Pixels above center line
    yc = fix(ycen)                                                 ;Integer pixel offsets
                                                                   ;for the center line

;;Overlay order cross-cut profiles with subpixel shifts predicted from order
;; location fits. Smooth the resulting oversampled slit function. For all but
;; the first order, temporarily replace the original data values near the
;; normalization boundary between modified and unmodified pixels. This is
;; necessary because slit function construction extends back into parts of
;; the image that have already been modified. As soon as the slit function
;; has been calculated, restore the normalized data values.
;    if onum gt ofirst then begin
;      saven = im(isave)                               ;save normalized values
;      im(isave) = saver                               ;restore raw values
;    endif

;
;if max(ymax) lt nrow then begin ;logic to handle orders that stray off edge
;
;Optimally extract the current order.
    mkslitf, im, scatter_below, yscatter_below $
           , scatter_above, yscatter_above $ 
           , ycen, y_lower_lim, y_upper_lim $
           , yslitf, slitf, binc, onum $
           , x_left_lim=cole0, x_right_lim=cole1, PLOT=iplot $
           , LAMBDA_SF=sf_smooth, LAMBDA_SP=sp_smooth, blz=blz $
           , OSAMPLE=osample $
           , SWATH_WIDTH=swath_width, MASK=mask, GAIN=gain, READN=readn $
           , /NORMALIZE, NORM_IMAGE=im_norm, ORDER_IMAGE=im_ordr $
           , THRESHOLD=thres, SWATH_BOUNDARIES=swath_boundaries $
           , WING_SMOOTH_FACTOR=wing_smooth_factor $
           , UNCERTAINTY=bunc
;    if onum gt ofirst then im(isave) = saven          ;normalized values again
    nbinc = n_elements(binc)                          ;number of slitf bins
    if keyword_set(iplot) then begin                  ;plot diagnostic
      plot, yslitf, slitf(*,0)>slitf(*,nbinc-1), /xsty, /yno
      for i=1, nbinc-1 do begin                       ;loop thru other bins
        oplot, yslitf, slitf(*,i), co=i+1             ;plot local slit funct
      endfor
;      junk = get_kbrd(1)
    endif
;    save,onum,yslitf,slitf,file='flat_sf_'+suffix(onum,2)+'.sav'

;Trim end pixels of slitf, since they are sometimes bad.
    ysfmin = min(yslitf, max=ysfmax)                  ;range of yslitf
    nysf = n_elements(yslitf)                         ;number of subpixels
    dysf = (ysfmax - ysfmin) / (nysf - 1)             ;subpixel size
    ntrim = ceil(1.0 / dysf)                          ;subpixels to trim
    yslitf = yslitf(ntrim:nysf-ntrim-1)               ;trim subpixel indexes
    slitf = slitf(ntrim:nysf-ntrim-1,*)               ;trim slit functions
    ysfmin = min(yslitf, max=ysfmax)                  ;range of yslitf
    nysf = n_elements(yslitf)                         ;number of subpixels

;Note that mkslitf returns trimmed arrays s(ncole), u(ncole) etc.

    s=blz(cole0:cole1)

;Fit the extracted spectrum with a (possibly broken) polynomial, iteratively
; rejecting any anomolous features in the spectrum.
;    fit=median(s,5)                                    ;Reject bad columns
;    fit=middle(fit,8,ITER=nloop,EPS=readn/gain,/poly)  ;Use smoothed blaze
;    fit=middle(fit,3.d0,EPS=readn/gain,ITER=nloop)     ;Use smoothed blaze
    fit = s                                           ;Use the spectrum from slit_func

;Diagnostic plot of fit to extracted flat field order.
    if keyword_set(iplot) then begin
      plot, s, /yno, xsty=3, ysty=3, chars=1.4 $
          , xtit='Column Number' $
          , tit='Order '+strtrim(onum, 2)
      oplot, fit, co=3
      empty                                   ;flush graphics buffer
;      junk = get_kbrd(1)
    endif

    blazecurve(cole0:cole1,onum-1) = fit      ;Save the fit

;stop

;;Extract unnormalized data values near edge of new normalization boundary.
;;These values will be needed temporarily, when the slit function for the
;; next order is computed. Otherwise, the slit function for the next order
;; will include already normalized pixels, which causes the low end of the
;; slit function to be near zero.
;    ysave = fix(ycen) + y_upper_lim                   ;last pixel to modify
;    isave = [ jx + ncol * (ysave-2) $                 ;Indices of a few rows
;            , jx + ncol * (ysave-1) $                 ;is a precaution against
;            , jx + ncol *  ysave    ]                 ;double normalization
;    saver = im(isave)                                 ;save raw data
;
;;Use slit function and fit to order to normalize flat.
;    if(cole0 gt 0) then begin                         ;skip the beginning of the order
;      for icol=0,cole0-1 do begin                     ;loop thru image columns
;        if(onum eq ofirst) then y0 = 0 $              ;starting row to change
;        else                    y0 = ((yprev(icol)+1)>0)<(nrow-1)
;        if(onum eq nord) then   y1 = nrow-1 $         ;ending row to change
;        else                    y1 = ((yc(icol) + y_upper_lim)<(nrow-1))>0
;        im(icol,y0:y1) = 1.0
;      endfor
;    endif
;
;    for icol=cole0,cole1 do begin                     ;loop thru image columns
;      if(exwd(0,onum-1) gt 1.5) then begin            ;starting row to change
;        y0 = round(yc(icol) - exwd(0,onum-1))         ;exwd is in pixels
;      endif else begin                                ;or a fraction
;        y0 = yc(icol) - round(y_lower_lim*exwd(0,onum-1))
;      endelse
;      y0=y0>0
;      if(exwd(1,onum-1) gt 1.5) then begin            ;ending row to change
;        y1 = round(yc(icol) + exwd(1,onum-1))         ;exwd is in pixels
;      endif else begin                                ;or a fraction
;        y1 = yc(icol) + round(y_lower_lim*exwd(1,onum-1))
;      endelse
;      y1=y1<(nrow-1)
;;      iy = y0 + indgen(y1 - y0 + 1)                   ;list of rows to change
;;      yint = (iy - ycen(icol) - ysfmin) / dysf        ;interpolation grid
;
;;Linearly interpolate slit functions onto desired rows in current column.
;      isort = sort(abs(binc - icol))                  ;sort by distance
;      if(n_elements(isort) le 1) then begin           ;For short partial orders
;        ic0=isort                                     ;that may fit in one swath
;        ic1=isort
;        bc0=cole0
;        bc1=cole1
;      endif else begin
;        ic0 = isort(0)                                ;index of closest bin
;        ic1 = isort(1)                                ;index of next closest bin
;        bc0 = binc(ic0)                               ;closest bin center
;        bc1 = binc(ic1)                               ;next closest center
;      endelse
;
;      sf0 = interpolate(slitf(*,ic0), yint)           ;interpolate onto rows
;      sf1 = interpolate(slitf(*,ic1), yint)           ;interpolate onto rows
;      sf = sf0 + (sf1-sf0)/(bc1-bc0) * (icol-bc0)     ;extra/interpolate icol
;      bgg =  ( scatter_above(icol) -  scatter_below(icol)) $
;           / (yscatter_above(icol) - yscatter_below(icol)) $
;           * (findgen(y1 - y0 + 1) + y0 - yscatter_below(icol)) $
;           + scatter_below(icol)
;
;;Diagnostic plot.
;      if 0 and (icol mod 50 eq 0 or icol eq ncol-1) then begin
;        plot, iy, im(icol,y0:y1), xsty=3, ysty=3 $
;            , chars=1.4, ps=2 $
;            , xtit='Row Number', ytit='ADU' $
;            , tit='Order ' + strtrim(onum, 2) $
;                 +',  Column ' + strtrim(icol, 2)
;        oplot, iy, sf * fit(icol-cole0) + bgg, co=2
;        oplot, iy, bgg, co=4
;        junk = get_kbrd(1)
;      endif
;
;;Normalize the current segment of the flat.
;      if(keyword_set(mask)) then msk = mask(icol,y0:y1) $
;      else                       msk = replicate(1B,y1-y0+1)
;
;;If FF signal is less than 10000 (S/N<100) than there is no point in flatfielding
;      low_signal = where(im(icol,y0:y1)*GAIN lt threshold, n_low_signal)
;      if(n_low_signal gt 0) then msk(low_signal)=0B
;      im(icol,y0:y1) = msk * im(icol,y0:y1) / (sf * fit(icol - cole0) + bgg) $
;                     + (1B - msk)
;
;;Fill any part we missed with 1.0, especially before initial order.
;      if y0 gt (yprev(icol)+1) then im(icol,yprev(icol)+1:y0) = 1.0
;
;;If done with last order, fill to final edge with 1.0.
;      if onum eq nord and y1 lt nrow-1 then im(icol,y1+1:nrow-1) = 1.0
;
;;Update border of modifications so far.
;      yprev(icol) = y1
;
;;End of normalization loop through columns.
;    endfor
;
;    if(cole1 lt ncol-1) then begin                    ;skip the end of the order
;      for icol=cole1+1,ncol-1 do begin                ;loop thru image columns
;        yc = fix(ycen(icol))                          ;order center in column
;        if(onum eq ofirst) then y0 = 0 $              ;starting row to change
;        else                    y0 = ((yprev(icol)+1)>0)<(nrow-1)
;        if(onum eq nord) then   y1 = nrow-1 $         ;ending row to change
;        else                    y1 = ((yc + y_upper_lim)<(nrow-1))>0
;;        y0 = ((yprev(icol)+1)>0)<(nrow-1)             ;starting row to change
;;        y1 = ((yc + y_upper_lim)<(nrow-1))>0          ;ending row to change
;        im(icol,y0:y1) = 1.0
;      endfor
;    endif
;;stop

;End of loop through orders.
  endfor

  im = im_norm                          ;copy normalized FF to the output
  blazecurve = blazecurve > 1.          ;avoid division by zero in the future
  return
end

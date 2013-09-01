pro mkslitf, im, scatter_below, yscatter_below, scatter_above, yscatter_above $
           , ycen, y_lower_lim, y_upper_lim, yslitf $
           , slitf, bincen, ord_num, PLOT=iplot $
           , X_LEFT_LIM=x_left_lim, X_RIGHT_LIM=x_right_lim $
           , LAMBDA_SF=lam_sf, LAMBDA_SP=lam_sp, SWATH_WIDTH=swath_width $
           , BLZ=blz, OSAMPLE=osample, MASK=mask, GAIN=gain, READN=readn $
           , NO_SCATTER=no_scatter, TELLURIC=telluric, FILENAME=filename $
           , SLIT_TILT=shear_x, NORMALIZE=normalize, NORM_IMAGE=im_norm  $
           , ORDER_IMAGE=im_ordr, THRESHOLD=threshold $
           , SWATH_BOUNDARIES=swath_boundaries $
           , WING_SMOOTH_FACTOR=wing_smooth_factor $
           , UNCERTAINTY=sunc

;Determines slit function along echelle order
;Input:
; im (array(ncol,nrow)) image containing echelle spectrum
; back (vector(ncol) vector containing background under the current order
; ymin (vector(ncol)) row numbers along bottom of region to map
; ycen (vector(ncol)) row numbers of zero point for slit function
; ymax (vector(ncol)) row numbers along top of region to map
;Output:
; yslitf (vector(nslitf)) subpixel row offsets for slitf
; sflit (array(nslitf,nbin)) subpixel slit functions
; bincen (vector(nbin)) column of bin centers
;History:
;18-Nov-96 Valenti  Wrote.
;21-Nov-97 Valenti  Adapted phx_slitf for use in echelle reduction package
;30-Mar-98 Valenti  Don't use polynomials to fit slit function. Instead
;                    median filter, then bin, then Gaussian smooth.
;05-May-98 CMJ      Back to using polynomials to fit slit function.  Seems
;                    to work OK with 18th order on the binned, median filtered
;                    slit function.
;13-Aug-98 CMJ      Put in logic to calculate both the smoothed medianed
;                    slit function and the polynomial one and choose between
;                    them based on standard deviations, giving the fit
;                    a little extra room for messiness.  I also force
;                    whatever method is chosen on the first time through to
;                    use throughout.  I have also increased the smoothing
;                    size for high resolution data.
;29-Nov-98 JAV      Check for rspec exactly zero and fudge it to be one,
;                    so as to avoid divide by zero and subsequent badness.
;08-Dec-98 JAV      Added an inital filtering of bad pixels in the oversampled
;                    slit function *before* binning. The rejection threshold
;                    is 3 times the mean absolute value of the difference
;                    between the oversampled (sf) and median filtered (medsf)
;                    slit functions. Indices of the good pixels are contained
;                    in "igd".
;09-Dec-98 CMJ      Added the conditional on binning the slit function back
;                    in to aviod inappropriate referencing of the variables.
;                    The result is some SF bins have no points in them, so
;                    added a later check which interpolates over these bins.
;06-Jun-99 JAV      Increased trace level of slit function type message from
;                    10 to 20 (suppressing the messages by default).
;30-Mar-01 NP       I will write down what I have done. I promise.
;10-May-06 NP       Time to keep my promise.
;                   Changes to the algorithm:
;                   - slit decomposition function is used to split the selected
;                     rectangular subset of an order to the spectrum (sp) and
;                     oversampled (OSAMPLE) slit illumination function (SFSM)
;                     The initial rectangle is compressed by shifting each column with
;                     +-1 when the center line ycen crosses pixel boundary and stored
;                     in sf. The model (compressed) rectangle is returned from slit_func
;                     as sfbin. Spectra from various swaths are combine in 1D array blz.
;                   - mask is used to hangle severe cosmetics problems of the detector
;
;                   New parameters and keywords:
;                   - scatter_above, scatter_below, yscatter_above, yscatter_below.
;                     Scatter light estimate between sp. orders produced by mkscatter.
;                     It is subtracted here when building sf in case FF normalization.
;                     For science spectra extraction it is done by mkscatter.
;                   - ord_num is the relative order number used for infomation purpose only
;                   - x_left_lim, x_right_lim. Starting/ending columns of the order.
;                   - lambda_sf. (Default 0). Smoothing parameter for the slit illumination
;                     function.
;                   - lambda_sp. (Default 0). Smoothing parameter for the spectrum (should
;                     be zero for the science spectra extraction).
;                   - osample. Integer! (Default 1). Oversampling for the accurate reconstruction
;                     of the slit illumination function.
;                   - mask. Mask for bad pixels (byte array, 1B - good, 0B - bad).
;                   - no_scatter. When set, mkslitf assumes that the scattered light was
;                     taken care of by mkscatter.
;                   - telluric. If set to specific integer larger than 1 is used as the
;                     offset from the order center line. The sky is then estimated by computing
;                     median signal between this offset and the upper/lower limit of the
;                     extraction window.
;                   - filename. Processed filename, if set is used for window identification.
;25-Feb-08 NP       Return swath boundaries if requested, compute the uncertainty vector inside slit_func
;01-Jul-08 NP
;

if n_params() lt 8 then begin
  print, 'syntax: mkslitf, im, scatter_below, yscatter_below, scatter_above, yscatter_above $'
  print,'               [, ycen, y_lower_lim, y_upper_lim, yslitf $'
  print,'               [, slitf, bincen, ord_num, PLOT=iplot $'
  print,'               [, X_LEFT_LIM=x_left_lim, X_RIGHT_LIM=x_right_lim $'
  print,'               [, LAMBDA_SF=lam_sf, LAMBDA_SP=lam_sp, SWATH_WIDTH=swath_width $'
  print,'               [, BLZ=blz, OSAMPLE=osample, MASK=mask, GAIN=gain, READN=readn $'
  print,'               [, NO_SCATTER=no_scatter, TELLURIC=telluric,FILENAME=filename $'
  print,'               [, SWATH_BOUNDARIES=swath_boundaries $'
  print,'               [, SLIT_TILT=shear_x]]]]]]]]]'
  retall
endif


  whenstop=1800

;Get image size.
  sz = size(im)                                                   ;variable info
  ncol = sz(1)                                                    ;number of columns
  nrow = sz(2)                                                    ;number of rows
  if(not keyword_set(x_left_lim) ) then x_left_lim  = 0           ;first column to extract
  if(not keyword_set(x_right_lim)) then x_right_lim = ncol-1L     ;last column to extract
  if(keyword_set(gain))  then CCD_gain  = gain  else CCD_gain=1.  ;CCD gain e- per 1 ADU
  if(keyword_set(readn)) then CCD_readn = readn else CCD_readn=0. ;Readout noise
  noise=CCD_readn/CCD_gain

  USE_2D_SLIT_FUNC=0B
  if(keyword_set(shear_x)) then begin
    shear=-poly(dindgen(ncol),shear_x[*,ord_num-1])
    delta_x=ceil(abs(max([y_lower_lim,y_upper_lim])*max(abs(shear))))
    USE_2D_SLIT_FUNC=1B
  endif

  msk = 0
  imask = 0
  if(keyword_set(mask)) then begin
    sz = size(mask)
    if(sz(0) eq 2 and sz(1) eq ncol and sz(2) eq nrow) then imask = 1 $
    else begin
      print,'MKSLITF: Your mask is defined but does not match your image'
      stop
    endelse
  endif

;Internal program parameters.
  if(not keyword_set(swath_width)) then begin   ;Estimate the
    i = uniq(long(ycen))                        ;Points of row crossing
    ni = n_elements(i)                          ;This is how many times this order
                                                ;crosses to the next row
    if(ni gt 1) then begin                      ;Curved order crosses rows
      i = total(i(1:ni-1)-i(0:ni-2))/(ni-1)
      nbin = ((round(ncol/i) / 3)  > 3) < 20    ;number of swaths along the order
    endif else begin                            ;Perfectly aligned orders
      nbin = (ncol/400) > 3                     ;Still follow the changes in PSF
    endelse
    nbin = nbin * (x_right_lim $                ;Adjust for the true order length
                 - x_left_lim + 1) / ncol
  endif else $                                  ;If swath width is preset
    nbin = round((x_right_lim - x_left_lim + 1.) / swath_width)>1

  yc = fix(ycen)
;Find columns that contain spectrum.
  imin = yc - y_lower_lim                       ;bottom row of order
  imax = yc + y_upper_lim                       ;top row of order

;Calculate boundaries of distinct slitf regions.
  ibound = (x_right_lim - x_left_lim) * $
           findgen(nbin+1) / nbin + x_left_lim  ;boundaries of bins
  ibeg = ceil(ibound(0:nbin-1))                 ;beginning of each bin
  ibeg(0) = x_left_lim
  iend = floor(ibound(1:nbin))                  ;end of each bin
  if(nbin gt 1) then begin
    iend = [ibeg(1:*)-1, x_right_lim]
    ibeg_half = (ibeg[1:nbin-1] + ibeg[0:nbin-2]) / 2 ;boundaries of overalapping steps
    iend_half = (iend[1:nbin-1] + iend[0:nbin-2]) / 2
    ibeg_half = [ibeg_half, ibeg] & ibeg_half = ibeg_half[sort(ibeg_half)]
    iend_half = [iend_half, iend] & iend_half = iend_half[sort(iend_half)]
  endif else begin
    iend=[x_right_lim]
    ibeg_half=[x_left_lim ,x_right_lim]
    iend_half=[x_right_lim,x_right_lim]
  endelse
  bincen = 0.5*(ibeg_half + iend_half)          ;center of each bin
  if(arg_present(swath_boundaries)) then begin  ;store swath bounaries
    swath_boundaries=intarr(2,2*nbin-1)         ;if requested
    swath_boundaries(0,*)=ibeg_half[0:2*nbin-2]
    swath_boundaries(1,*)=iend_half[0:2*nbin-2]
  endif

;Initialize default parameters and arrays.
  if(keyword_set(osample)) then osamp=osample $ ;slitf pixels / real pixel
  else osamp = 10
  if(not keyword_set(lam_sf)) then lambda_sf=1. else lambda_sf=lam_sf
  if(not keyword_set(lam_sp)) then lambda_sp=0  else lambda_sp=lam_sp

  irow = findgen(nrow)                          ;indices of all rows
  ycene = ycen(x_left_lim:x_right_lim)
  nysf = y_upper_lim + y_lower_lim + 1L         ;subpixel range required
  yslitf0 = -y_lower_lim                        ;minimum value for yslitf
  yslitf1 =  y_upper_lim                        ;maximum value for yslitf

  blz = fltarr(ncol)

;/////////////////////////////////
;sfsm1=dblarr((nysf+1)*osample+1,nbin)
;sp1=dblarr(iend(0)-ibeg(0)+2,nbin)  

;Perform slit decomposition within each swath stepping through the order with
;half swath width. Spectra for each decomposition are combined with lenar weights.
  for ihalf=0, 2*nbin-2 do begin                             ;loop thru swaths
    ib = ibeg_half(ihalf)                                    ;left column
    if(USE_2D_SLIT_FUNC and ihalf gt 0)        then ib = ib - delta_x
    ie = iend_half(ihalf)                                    ;right column
    if(USE_2D_SLIT_FUNC and ihalf lt 2*nbin-2) then ie = ie + delta_x
    nc = ie - ib + 1                                         ;number of columns

;Load slit function data into vectors.
    nsf = nc * nysf                                          ;# slit func points
    j0 = lonarr(nc)
    j1 = lonarr(nc)
    sf = fltarr(nc,nysf)
    sfpnt = fltarr(nsf)
    ysfpnt = fltarr(nsf)
    if(imask) then msk = bytarr(nc,nysf)
    for j=0, nc-1 do begin                                   ;loop thru columns in region
      icen = yc(ib+j)                                        ;row closest to peak
      k0 = icen + yslitf0                                    ;lowest row to consider
      k1 = icen + yslitf1                                    ;highest row to consider
      j0(j) = j*nysf                                         ;begining of storage area
      j1(j) = j0(j) + k1 - k0                                ;new logic (works at edge)
      if(keyword_set(no_scatter)) then begin
        ssf = im(ib+j,k0:k1)
      endif else begin
        dy_scatter =   (dindgen(nysf) + k0 - yscatter_below(ib+j)) $
                     / (yscatter_above(ib+j) - yscatter_below(ib+j))
        scatter = (scatter_above(ib+j)-scatter_below(ib+j))* $ ;Interpolate background
                   dy_scatter+scatter_below(ib+j)
        ssf = (im(ib+j,k0:k1) - scatter)>0                     ;subtract scatter
      endelse
      sf(j,*) = ssf
      sfpnt(j0(j):j1(j)) = ssf
      ysfpnt(j0(j):j1(j)) = irow(k0:k1) - ycen(ib+j)
      if(imask) then msk(j,*) = mask(ib+j,k0:k1)
    endfor

;    j=where(sf gt 200, nj)
;    if(nj gt 0) then msk(j)=0B

;if(ord_num ge 65 and ib eq 0) then stop
    if(keyword_set(telluric)) then begin
      if(telluric gt 0 and telluric lt nysf/2) then tel_lim=telluric else tel_lim=5<(nysf/3)
      tel=total(sf,1)
      itel=indgen(nysf)
      itel=itel(where(abs(itel-nysf/2) ge tel_lim))
      tel=sf(*,itel)
      sc=dblarr(nc)
      for itel=0,nc-1 do sc(itel)=median(reform(tel(itel,*)))
;if(keyword_set(istop)) then stop
;stop
      sf = sf - sc#replicate(1,yslitf1-yslitf0+1)
      sfpnt = reform(transpose(sf),n_elements(sfpnt))
    endif

    jgood=0
    if(USE_2D_SLIT_FUNC) then begin
      slit_func_2d,sf,ycen(ib:ie)-yc(ib:ie),sp,sfsm $
               ,delta_x,shear,NOISE=noise $
               ,OVERSAMPLE=osamp,IM_OUT=sfbin,LAMBDA_SF=lambda_sf $
               ,LAMBDA_SP=lambda_sp,USE_COL=jgood,BAD=jbad,MASK=msk $
               ,WING_SMOOTH_FACTOR=wing_smooth_factor $
               ,UNCERTAINTY=unc
    endif else begin
      slit_func,sf,ycen(ib:ie)-yc(ib:ie),sp,sfsm,NOISE=noise $
               ,OVERSAMPLE=osamp,IM_OUT=sfbin,LAMBDA_SF=lambda_sf $
               ,LAMBDA_SP=lambda_sp,USE_COL=jgood,BAD=jbad,MASK=msk $
               ,WING_SMOOTH_FACTOR=wing_smooth_factor $
               ,UNCERTAINTY=unc
    endelse
;if(ib gt 1000) then stop
;if(ib gt 1500 and ord_num ge 26) then stop

;/////////////////////////////////
;sfsm1(*,i)=sfsm
;sp1(0:nc-1,i)=sp
;if(ib gt whenstop) then begin
;  whenstop=10000L
;  save,sfsm,file=strmid(filename,0,7)+'_lsf.sav'
;;  stop
;endif

    weight = replicate(1.,nc)          ;weight vector for combining the overalaping part
    if(ihalf gt 0) then weight(0:nc/2) = findgen(nc/2+1)/(nc/2)
    oweight = 1. - weight

;In case we do FF normalization replace the original image by the
;ratio of sf/sfbin where number of counts is larger than threshold
;and with 1 elsewhere
    scale=1
    if(keyword_set(normalize)) then begin
      if(imask) then ii = where(sfbin gt threshold/gain and msk eq 1B, nii) $
      else           ii = where(sfbin gt threshold/gain, nii)
      sss = replicate(1.,nc,nysf)                 ;1     - by default
      ddd = fltarr(nc,nysf)                       ;save orders, helps analysing ghosts
      if(nii gt 0) then sss[ii] = sf(ii)/sfbin(ii);ratio - when high signal
      ddd = sfbin
      if(ihalf gt 0) then begin
;
;Here is the layout to understand the lines below
;
;        1st swath    3rd swath    5th swath      ...
;     /============|============|============|============|============|
;
;               2nd swath    4th swath    6th swath
;            |------------|------------|------------|------------|
;            |.....|
;            overlap
;
;            +     ******* 1
;             +   *
;              + *
;               *            weights (+) previous swath, (*) current swath
;              * +
;             *   +
;            *     +++++++ 0


        overlap=iend_half[ihalf-1] - ibeg_half[ihalf] + 1
;        scale=median((blz[ibeg_half[ihalf]:iend_half[ihalf-1]]>1)/(sp[0:overlap-1]>1))
;print,scale
;stop
        scale=1
        if(nii gt 0) then sss[ii]=sss[ii]/scale
        sp=sp*scale
      endif else begin
        nc_old = nc                               ;previous swath width
        sss_old = fltarr(nc,nysf)                 ;dummy normalized ff from previous swath
        ddd_old = fltarr(nc,nysf)                 ;dummy order shape ff from previous swath
        overlap = nc_old + 1
      endelse
;
;This loop is complicated because swaths do overlap to ensure continuity of the
;spectrum.
      if(ihalf eq 0)             then ncc=ibeg_half[1] - ibeg_half[0] $
      else                            ncc=overlap
;if(ord_num eq 24 and ihalf ge 27) then stop
      for j=0, ncc-1 do begin                     ;loop thru columns in region
        icen = yc(ib+j)                           ;row closest to peak
        k0 = icen + yslitf0                       ;lowest row to consider
        k1 = icen + yslitf1                       ;highest row to consider
        j0(j) = j * nysf                          ;begining of storage area
        j1(j) = j0(j) + k1 - k0                   ;new logic (works at edge)
        jj = nc_old - ncc + j             ;column number in the previous swath
        im_norm(ib+j,k0:k1) = sss_old(jj,*) * oweight(j) $;replace the original
                            + sss    (j ,*) *  weight(j)
        im_ordr(ib+j,k0:k1) = ddd_old(jj,*) * oweight(j) $;replace the original
                            + ddd    (j ,*) *  weight(j)
      endfor
;if(ord_num eq 24 and ihalf ge 27) then stop
      if(ihalf eq 2*nbin-2) then begin            ;finish the very last swath   
        for j=ncc,nc-1 do begin                   ;loop thru columns in region
          icen = yc(ib+j)                         ;row closest to peak
          k0 = icen + yslitf0                     ;lowest row to consider
          k1 = icen + yslitf1                     ;highest row to consider
          j0(j) = j * nysf                        ;begining of storage area
          j1(j) = j0(j) + k1 - k0                 ;new logic (works at edge)
          im_norm(ib+j,k0:k1) = sss(j ,*)         ;replace the original
          im_ordr(ib+j,k0:k1) = ddd(j ,*)         ;replace the original
        endfor
      endif 
      nc_old  = nc
      sss_old = sss
      ddd_old = ddd
    endif ; End of FF normalization part

;if(keyword_set(istop)) then stop

;if ord_num eq 6 then stop

    nbad=n_elements(jbad)
    if(nbad eq 1) then if(jbad eq -1) then nbad=0

;    sp=total(sf,2)

    for j=0,nc-1 do begin          ; Normalize sfpnt for the plot
      sfpnt(j0(j):j1(j)) = sfpnt(j0(j):j1(j)) / (sp(j)>1.)
      sfbin(j,*) = sfbin(j,*) / (sp(j)>1.)
    endfor
    
    if(USE_2D_SLIT_FUNC) then begin
      if(ihalf gt 0 and ihalf lt 2*nbin-2) then begin   ; Copy the extracted spectrum
        blz(ib+delta_x:ie-delta_x) = blz(ib+delta_x:ie-delta_x) * oweight(delta_x:nc-delta_x-1) $
                                   + sp(delta_x:nc-delta_x-1)   *  weight(delta_x:nc-delta_x-1)
      endif else if(ihalf eq 0) then begin
        blz(ib        :ie-delta_x) = sp(0:nc-delta_x-1)
      endif else if(ihalf eq 2*nbin-2) then begin
        blz(ib+delta_x:ie)         = blz(ib+delta_x:ie) * oweight(delta_x:nc-1) $
                                   + sp(delta_x:nc-1)   *  weight(delta_x:nc-1)
      endif
    endif else begin
      blz(ib:ie) = blz(ib:ie) * oweight + sp * weight
    endelse
    if(arg_present(sunc)) then sunc(ib:ie) = sunc(ib:ie) * oweight + unc * weight

    if(ihalf eq 0) then begin
      nslitf = n_elements(sfsm)
      yslitf = yslitf0 + (findgen(nslitf)-0.5)/osamp-1.5 ;final subpixel scale
      slitf = fltarr(nslitf, 2*nbin-1)                   ;init final slit function
    endif

    sfsm2 = reform(transpose(sfbin),nsf)
    j = sort(ysfpnt)
;    if(keyword_set(istop)) then stop
;    if(nbad gt 0) then for k=0L,nbad-1 do jbad(k) = (where(j eq jbad(k)))(0)
; The line below does exactly the same as the line above. I have no clue why. NP
    if(nbad gt 0) then jbad = (sort(j))(jbad)
    ysfpnt = ysfpnt(j)
    sfpnt = sfpnt(j)
    sfsm2 = sfsm2(j)
;    j=uniq(ysfpnt)
;    jj=where(yslitf gt min(ysfpnt(j)) and yslitf le max(ysfpnt(j)))
;    sfsm(jj) = spl_interp(ysfpnt(j), sfsm2(j), spl_init(ysfpnt(j), $
;                          sfsm2(j)), yslitf(jj))
    slitf(*,ihalf) = sfsm/total(sfsm)*osamp         ;save slit function

;    plot,blz&oplot,indgen(ie-ib+1)+ib,sp,col=2
;    stop
;    sf0=sf & sfbin0=sfbin & sfsm0=sfsm

;Plot diagnostics.
    if keyword_set(iplot) then begin
      colors
      !p.multi=[0,2,2]

      pscale = mean(sp)
      sfplot=smooth(sfsm*pscale,(osamp-1)>1)/scale
      plot, ysfpnt, sfpnt*pscale, ps=3, /xsty, ysty=3, yr=minmax(sfpnt*pscale) $
          , /NODATA, title='Order '+strtrim(ord_num,2)+ $
          ', Columns '+strtrim(ib,2)+' through '+strtrim(ie,2)
      oplot, ysfpnt, sfpnt*pscale, ps=3
      if(nbad gt 0) then oplot, [ysfpnt(jbad)], [sfpnt(jbad)*pscale], ps=1, co=3
      ymiddle=0.5*(!y.crange(0)+!y.crange(1))
      oplot, yslitf, sfplot, thick=2, col=2
      oplot, !x.crange, [0.,0.], co=4

      plot, ysfpnt, sfpnt*pscale, ps=3, /xsty, ysty=3, yr=minmax(sfsm*pscale) $
          , /NODATA, title='Order '+strtrim(ord_num,2)+ $
          ', Columns '+strtrim(ib,2)+' through '+strtrim(ie,2)
      oplot, ysfpnt, sfpnt*pscale, ps=3
      if(nbad gt 0) then oplot, [ysfpnt(jbad)], [sfpnt(jbad)*pscale], ps=1, co=3
      ymiddle=0.5*(!y.crange(0)+!y.crange(1))
      oplot, yslitf, sfplot, thick=2, col=2
      oplot, !x.crange, [0.,0.], co=4

      plot, ysfpnt,(sfpnt - sfsm2)*pscale, ps=3, xs=1, title='Data - Fit'
      if(nbad gt 0) then oplot, [ysfpnt(jbad)], [(sfpnt(jbad) - sfsm2(jbad))*pscale] $
                              , ps=1, col=3
      oplot, !x.crange,  [0.,0.], co=4
      if(keyword_set(scatter_above) and keyword_set(scatter_below)) then $
        poffset = mean(scatter_below(ib:ie) + scatter_above(ib:ie))*0.5 $
      else poffset = 0
      oplot, yslitf,  sqrt(((sfsm*pscale/scale + poffset + CCD_readn^2)/CCD_gain)>0),col=2,thick=2
      oplot, yslitf, -sqrt(((sfsm*pscale/scale + poffset + CCD_readn^2)/CCD_gain)>0),col=2,thick=2
      oplot,  (!x.crange(1)-!x.crange(0))*[0.03,0.1]+!x.crange(0) $
           ,  (!y.crange(1)-!y.crange(0))*[0.06,0.06]+!y.crange(0), col=2,thick=2
      xyouts, (!x.crange(1)-!x.crange(0))*0.1+!x.crange(0) $
            , (!y.crange(1)-!y.crange(0))*[0.05,0.05]+!y.crange(0) $
            , '1!7r!3 + scatter + read noise'

      yr = sqrt(max(sfsm*pscale + poffset + CCD_readn^2)/CCD_gain)*2.
      yr = [-yr,yr]
      plot, ysfpnt,(sfpnt - sfsm2)*pscale, ps=3, xs=1, ys=3 $
;          , yr= minmax(sfpnt(jgood) - sfsm2(jgood))*pscale, title='Fit - Data'
          , yr = yr, title='Data - Fit'
      if(nbad gt 0) then oplot, [ysfpnt(jbad)], [(sfpnt(jbad) - sfsm2(jbad))*pscale] $
                              , ps=1, col=3
      oplot, !x.crange,  [0.,0.], co=4
      oplot, yslitf,  sqrt(((sfsm*pscale/scale + poffset + CCD_readn^2)/CCD_gain)>0),col=2,thick=2
      oplot, yslitf, -sqrt(((sfsm*pscale/scale + poffset + CCD_readn^2)/CCD_gain)>0),col=2,thick=2
      oplot,  (!x.crange(1)-!x.crange(0))*[0.03,0.1]+!x.crange(0) $
           ,  (!y.crange(1)-!y.crange(0))*[0.06,0.06]+!y.crange(0), col=2,thick=2
      xyouts, (!x.crange(1)-!x.crange(0))*0.1+!x.crange(0) $
            , (!y.crange(1)-!y.crange(0))*[0.05,0.05]+!y.crange(0) $
            , '1!7r!3 + scatter + read noise'

      if(keyword_set(filename)) then xyouts,0.5,0.5,filename,align=0.5,/NORM
      empty
      wshow,0
      !p.multi=[0,0,0]
    endif
  endfor

; Combine fitting errors and the Poisson shot noise
  sunc=sqrt((sunc+blz)>1.)

;if(not keyword_set(normalize) and ord_num ge 2) then begin
;  rdsk, sold, 'Reduced_Green/rj13.126', 1
;  sum=blz&for i=0,ncol-1 do sum(i)=total(im(i,(ycen(i)-10)>0:(ycen(i)+10)<(nrow-1)))
;  plot,blz,xs=1,tit=ord_num,psym=10,xr=[1000,2400]
;  if(ord_num ge 2) then begin
;    oplot,sold(*,ord_num-2),col=2,psym=10
;    oplot,(blz-sold(*,ord_num-2))*20.+1.d4
;    oplot,(sum-sold(*,ord_num-2))*20.+1.d4,col=4
;    oplot,sunc*20.+1.d4,col=3
;  endif
;  for ibb=0,nbin-1 do oplot,ibeg(ibb)+[0,0],!y.crange,line=2
;  oplot,!x.crange,1.d4+[0,0],line=4
;
;  bin=50.
;  n_bin=(max(sold(*,ord_num-2)-blz)-min(sold(*,ord_num-2)-blz))/bin
;  xhold_new=dindgen(n_bin-1)*bin+min(sold(*,ord_num-2)-blz)
;  hold_new=histogram(sold(*,ord_num-2)-blz,bin=bin)
;;  plot,xhold_new,hold_new,psym=10,xr=[-2000,2000]
;
;  n_bin=(max(sum-blz)-min(sum-blz))/bin
;  xhsum_new=dindgen(n_bin-1)*bin+min(sum-blz)
;  hsum_new=histogram(sum-blz,bin=bin)
;;  oplot,xhold_new,hold_new,psym=10,col=2
;
;  ssss=get_kbrd(1)
;endif
end

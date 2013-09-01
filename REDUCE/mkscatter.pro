pro mkscatter, im, orders, back, yback, colrange=colrange, LAMBDA_SF=lam_sf $
             , LAMBDA_sp=lam_sp, SWATH_WIDTH=swath_width, OSAMPLE=osample $
             , MASK=mask, GAIN=gain, READN=rdnoise, SUBTRACT=subtract $
             , POLY=pol, DATA=back_data,ORDER_HEIGHT=order_height $
             , DEBUG=debug,POLARIZ=polariz
;
; This subroutine is still far from perfect. NP wrote it in 2002 and since it
; was a pain in the back. The problem is that that it works very smoothlu with
; a good data set but once we get to the low S/N, overlapping orders and bad
; cosmetics things starting to fail. Typical failures are connected to the
; the identification of the walls of adjucent orders and the bottom.
;
; History of modifications:
;  2006-09-26 (NP) the condition for local maxima ID is set to
;                  (positive left der and non-positive right der) or
;                  (non-negative left der and negative right der) 


;Get image size.
  sz = size(im)                 ;variable info
  ncol = sz(1)                  ;number of columns
  nrow = sz(2)                  ;number of rows
  nord = n_elements(orders(0,*));number of orders

  msk = 0
  imask = 0
  if(keyword_set(mask)) then begin
    sz = size(mask)
    if(sz(0) eq 2 and sz(1) eq ncol and sz(2) eq nrow) then imask = 1 $
    else begin
      print,'MKSCATTER: Your mask is defined but does not match your image'
      stop
    endelse
  endif

  if(not keyword_set(gain))    then CCD_gain = 1. $       ; CCD gain
                               else CCD_gain = gain
  if(not keyword_set(rdnoise)) then CCD_rdnoise = 0. $    ; CCD rdnoise
                               else CCD_rdnoise = rdnoise/CCD_gain
  if(keyword_set(osample)) then osamp=osample $ ;slitf pixels / real pixel
  else                          osamp = 10
  if(not keyword_set(lam_sf)) then lambda_sf=1.
  if(keyword_set(order_extra_width)) then begin
    if(order_height gt 0. and order_height lt 2.) then $
      extra_offset=order_height
  endif else extra_offset=1.

  nloop = 10                    ;smoothing loop for scattered light

;Initialize arrays.
  xcol = dindgen(ncol)                          ;indices of all columns
  back = fltarr(ncol,nord+1)                    ;fitted scattered light model
  back_data = back                              ;scattered light data
  yback = fltarr(ncol,nord+1)                   ;scattered light coordinates
  ycen1 = poly(xcol, orders(*,0))               ;the shape of the starting order
  dbg=0

;Loop through orders
  for ord=1,nord-1 do begin
    if(not keyword_set(swath_width)) then width = 400 $    ; CCD rdnoise
                                     else width = swath_width<400
    ycen0 = ycen1
    ycen1 = poly(xcol, orders(*,ord))           ;the shape of the next order
    ibeg = colrange(0,ord-1)>colrange(0,ord)    ;left boundary to consider
    iend = colrange(1,ord-1)<colrange(1,ord)    ;right boundary to consider
    width = ((iend-ibeg+1)*0.1)>width           ;width to consider (swath < width < iend-ibeg+1)
    width = floor(width*0.5)                   ;half-width to avoid rounding errors
    height = mean(ycen1(ibeg:iend)-ycen0(ibeg:iend)) ;mean height
    height = round(height*0.5*extra_offset)>3   ;half-height to avoid rounding errors

    x_left_lim  = ibeg                          ;store for future use
    x_right_lim = iend

    icen = round(0.5*(ibeg+iend))               ;central column
    ibeg = (icen-width)>ibeg                    ;starting column of the region of interest
    iend = (icen+width)<iend                    ;ending column of the region of interest
    ycen = 0.5*(ycen1+ycen0)                    ;central line of the troff
    ymin = fix(ycen)-height                     ;bottom boundary of the region of interest
    ymax = fix(ycen)+height                     ;top boundary of the region of interest

    nc = 2*width+1                              ;dimensions of the region of interest
    nr = 2*height+1

;Copy the region of interest to sf
    sf = fltarr(nc,nr)
    if(imask) then msk = bytarr(nc,nr)
    for j=0,nc-1 do begin
      sf(j,*) = im(ibeg+j,ymin(j+ibeg):ymax(j+ibeg))
      if(imask) then msk(j,*) = mask(ibeg+j,ymin(j+ibeg):ymax(j+ibeg))
    endfor
;    display,im,/log&oplot,[ibeg,ibeg,iend,iend,ibeg], $
;                          [ycen[ibeg]-height,ycen[ibeg]+height, $
;                           ycen[iend]+height,ycen[iend]-height, $
;                           ycen[ibeg]-height]&wait,0.1

;    noise=CCD_readn*2.

    sf=sf-min(sf)
    slit_func,sf,ycen(ibeg:iend)-fix(ycen(ibeg:iend)),sp,sfsm,NOISE=noise $
             ,OVERSAMPLE=osamp,LAMBDA_SF=lam_sf,LAMBDA_SP=lam_sp $
             ,MASK=msk,IM_OUT=sfbin

    nslitf = n_elements(sfsm)
    yslitf = (findgen(nslitf)-0.5)/osamp-1.5-height            ;final subpixel scale

    dev = CCD_rdnoise/CCD_gain/sqrt(total(sp)/(iend - ibeg + 1.)) ;Typical pixel scatter normalized to counts
    var = stddev(sfsm)
;
; This is an improved version which is near the final form
;
    jback = bottom(sfsm<median(sfsm), 1, eps=dev, /POLY)       ;Fit the bottom with a smooth curve
    jback = (sfsm - jback) > 0.                                ;Find all positive spikes

;    k = where(jback(1:nslitf-2) gt jback(0:nslitf-3) and $     ;above the middle
;              jback(1:nslitf-2) gt jback(2:nslitf-1) and $
;              sfsm(1:nslitf-2) gt median(sfsm)+0.5*var, nback) + 1
;              
;    m1=where(yslitf(k) lt 0, n1) ;Part of the gap which is below the central line
;    m2=where(yslitf(k) gt 0, n2) ;Part of the gap which is above the central line
;    if(n1 eq 0 or n2 eq 0) then begin ;This happens e.g. if the two adjacent
;                                      ;orders have dramatically different level
;                                      ;of signal
;      k = where(jback(1:nslitf-2) gt jback(0:nslitf-3) and $     ;above the middle
;                jback(1:nslitf-2) gt jback(2:nslitf-1) and $
;                sfsm(1:nslitf-2) gt middle(sfsm, 1, eps=dev, /POLY), nback) + 1
;      m1=where(yslitf(k) lt 0, n1)
;      m2=where(yslitf(k) gt 0, n2)
;    endif
;    if(n1 eq 0 or n2 eq 0) then begin ;This can still happen of the overall
;                                      ;signal level is way too low
;      k = where(jback(1:nslitf-2) gt jback(0:nslitf-3) and $     ;above the middle
;                jback(1:nslitf-2) gt jback(2:nslitf-1) and $
;                sfsm(1:nslitf-2) gt middle(sfsm, 1, eps=dev, /POLY), nback) + 1
;      m1=where(yslitf(k) lt 0, n1)
;      m2=where(yslitf(k) gt 0, n2)
;    endif

;===========================================================================
    k = where(sfsm gt median(sfsm)+0.5*var, nback)
    if(nback gt 0) then begin
      m1=where(yslitf(k) lt 0, n1) ;Part of the gap which is below the central line
      m2=where(yslitf(k) gt 0, n2) ;Part of the gap which is above the central line
    endif else begin
      n1=0
      n2=0
    endelse
    if(n1 eq 0 or n2 eq 0) then begin ;This happens e.g. if the two adjacent
                                      ;orders have dramatically different level
                                      ;of signal
      k = where(sfsm gt middle(sfsm, 1, eps=dev, /POLY), nback)
      m1=where(yslitf(k) lt 0, n1)
      m2=where(yslitf(k) gt 0, n2)
    endif
    if(n1 eq 0 or n2 eq 0) then begin ;This happens e.g. if the two adjacent
                                      ;orders have dramatically different level
                                      ;of signal
      k = where(sfsm gt middle(sfsm, 1.d3, eps=dev), nback)
      m1=where(yslitf(k) lt 0, n1)
      m2=where(yslitf(k) gt 0, n2)
    endif
    ss=[0,jback,0]
    kk = where((ss(k+1) ge ss(k) and ss(k+1) gt ss(k+2)) $     ;Select local maxima
            or (ss(k+1) gt ss(k) and ss(k+1) ge ss(k+2)))
    k = k(kk)
              
    m1=where(yslitf(k) lt 0, n1) ;Part of the gap which is below the central line
    m2=where(yslitf(k) gt 0, n2) ;Part of the gap which is above the central line
;===========================================================================


    if(n1 eq 0 or n2 eq 0) then begin ;Ultimate method
;      plot,yslitf,sfsm,xs=1,title='Order #'+strtrim(ord+1,2)
;      oplot,yslitf,middle(sfsm, 1, eps=dev, /POLY)
;      oplot,yslitf,middle(sfsm, 1.d3, eps=dev)
;      oplot,!x.crange,median(sfsm)+0.5*var+[0,0],line=2
;      oplot,yslitf(k),sfsm(k),psym=2
      print,'MKSCATTER: Failed finding interorder boundaries. ' $
           +'Using 5 pixels around central line.'
      k=(sort(abs(yslitf)))(0)+[-3,3]
      nback=n_elements(k)
      dbg=1
    endif

    if(nback lt 2 or (nback ge 2 and $           ; Check if we have only noise
      (min(yslitf(k)) gt  height*0.5 or $
       max(yslitf(k)) lt -height*0.5))) then begin
      n1=min(abs(yslitf+height*0.5),k1)          ; In this case just select
      n2=min(abs(yslitf-height*0.5),k2)          ; the central half
      k=[k1,k2]
      nback=2
    endif

    if(jback(0) ge jback(1)) then begin
      k = [0,k]
      nback = nback + 1
    endif
    if(jback(nslitf-1) ge jback(nslitf-2)) then begin
      k = [k,nslitf-1]
      nback = nback + 1
    endif

    if(keyword_set(debug)) then begin
      plot,yslitf,sfsm,xs=1,title='Order #'+strtrim(ord+1,2)
      oplot,yslitf,middle(sfsm, 1, eps=dev, /POLY)
      oplot,yslitf,middle(sfsm, 1.d3, eps=dev)
      oplot,!x.crange,median(sfsm)+0.5*var+[0,0],line=2
      oplot,yslitf(k),sfsm(k),psym=2
    endif

    jback = where(yslitf(k) lt 0, n1)
;    m1 = max(sfsm(k(jback)), imax1)
    m1 = max(yslitf(k(jback)), imax1)
    imax1 = k(jback(imax1))

    jback = where(yslitf(k) gt 0, n2)
;    m2 = max(sfsm(k(jback)), imax2)
    m2 = min(yslitf(k(jback)), imax2)
    imax2 = k(jback(imax2))

    k=[imax1,imax2]
    if(nback eq 0) then begin
      print,'MKSCATTER: Failed to find the boundaries of the inter-order troff'
      print,'           Using default'
      print,'Order: '+strtrim(ord,2)+'  X-range: ('+strtrim(ibeg,2)+','+strtrim(iend,2)+')'
      k = [0, n_elements(nslitf-1)]
    endif

;if(dbg eq 1) then stop
    jback = bottom(sfsm(k(0):k(1))<median(sfsm), 1, eps=dev, /poly)  ;Fit bottom with a straight line
    iback = where(sfsm(k(0):k(1)) le jback+CCD_rdnoise, nback)       ;Find all the points below
    if(nback le 5) then begin
      print,'MKSCATTER: Major error in the order format: could not detect inter-order troff'
      plot,yslitf,sfsm,xs=1,title='Order #'+strtrim(ord,2)
      oplot,yslitf,bottom(sfsm<median(sfsm), 1, eps=dev, /POLY)
      oplot,yslitf(k),sfsm(k),psym=2
      stop
    endif
    iback = iback + k(0)
    imax1 = min(iback)
    imax2 = max(iback)
    sf1 = sfsm(iback)
    iback = 0
    jback = 0

    step = (max(sf1)-min(sf1))/(((n_elements(sf1)/10)<200)>10)
    h = histogram(sf1,bin=step)
    nh = n_elements(h)
    hmax = max(h,ihmax)
   
;    plot,step*dindgen(nh)+min(sf1),h,xs=1,psym=10 $
;        ,xr=[-20.*step,20*step]+step*ihmax+min(sf1)
    i0 = where(h(0:ihmax) lt 0.1*hmax, n1)
    if(n1 gt 0) then i0 = max(i0) else i0=0
    i1 = where(h(ihmax:nh-1) lt 0.1*hmax, n2)
    if(n2 gt 0) then i1 = min(i1)+ihmax else i1=nh-1
    ii = where(sf1-min(sf1) ge step*i0 and sf1-min(sf1) lt step*i1, nii1)+imax1
    if(nii1 le 0) then begin
      print,'MKSCATTER: Could not detect background points between orders '+ $
            strtrim(ord-1,2)+' and '+strtrim(ord,2)
    endif else begin
      y1 =  ceil(min(yslitf(ii)))                        ;bottom boundary of the troff
      y2 = floor(max(yslitf(ii)))                        ;top boundary of the troff
    endelse

    if(keyword_set(debug)) then begin
      oplot,[y1,y1],!y.crange & oplot,[y2,y2],!y.crange
      wait,0.1
;      if(ord gt 10) then stop
    endif

    if(nii1 le 0) then begin
      print,'MKSCATTER: Could not detect any background points'
      stop
    endif

    yback(*,ord) = ycen                                  ;scattered light coordinates

;    display,im,/LOG,yr=[min(ycen0+y1),max(ycen1+y2)]
;    oplot,ycen+y1,col=0,line=2
;    oplot,ycen+y2,col=0,line=3
;    oplot,ycen0,col=0
;    oplot,ycen1,col=0
;    stop

    for j=0,ncol-1 do begin                              ;scattered light in col j
      yy1 = round(ycen(j)+y1)>0                          ;first row
      yy2 = round(ycen(j)+y2)<(nrow-1)                   ;last row

      if(yy2-yy1 gt (0.3*(y2-y1)+1)) then begin          ;check if inside the frame
        scatter = reform(im(j,yy1:yy2))                  ;scattered light profile
        if(imask) then begin
          k = where(mask(j,yy1:yy2) eq 1B, nk)           ;check for bad pixels
          if(nk lt (yy2-yy1+1) and nk ge 2) then begin   ;if found, interpolate
            scatter = interpol(scatter(k), k, findgen(yy2-yy1+1))
          endif
        endif
        back_data(j,ord) = median(scatter)               ;take median
      endif else begin
        back_data(j,ord) = -10000                        ;bad point
      endelse

      if(ord eq 1) then begin                            ;for the first order try
        yy = ycen0(j) - (ycen(j) - ycen0(j))             ;to find background below
        yy1 = round(yy+y1)>0                             ;first row
        yy2 = round(yy+y2)<(nrow-1)                      ;last row
        yback(j,ord-1) = yy                              ;scattered light coordinates

        if(yy2-yy1 gt (0.3*(y2-y1)+1)) then begin        ;check if inside the frame
          scatter = reform(im(j,yy1:yy2))                ;scattered light profile
          if(imask) then begin
            k = where(mask(j,yy1:yy2) eq 1B, nk)         ;check for bad pixels
            if(nk lt (yy2-yy1+1) and nk ge 2) then begin ;if found, interpolate
              scatter = interpol(scatter(k), k, findgen(yy2-yy1+1))
            endif
          endif
          back_data(j,ord-1) = median(scatter)           ;take median
        endif else begin
          back_data(j,ord-1) = -10000
        endelse
      endif else if(ord eq nord-1) then begin            ;for the last order try
        yy = ycen1(j) + (ycen1(j) - ycen(j))             ;to find background above
        yback(j,ord+1) = yy                              ;scattered light coordinates
        yy1 = round(yy+y1)>0                             ;first row
        yy2 = round(yy+y2)<(nrow-1)                      ;last row

        if(yy2-yy1 gt (0.3*(y2-y1)+1)) then begin        ;check if inside the frame
          scatter = reform(im(j,yy1:yy2))                ;scattered light profile
          if(imask) then begin
            k = where(mask(j,yy1:yy2) eq 1B, nk)         ;check for bad pixels
            if(nk lt (yy2-yy1+1) and nk ge 2) then begin ;if found, interpolate
              scatter = interpol(scatter(k), k, findgen(yy2-yy1+1))
            endif
          endif
          back_data(j,ord+1) = median(scatter)           ;take median
        endif else begin
          back_data(j,ord+1) = -10000
        endelse
      endif
    endfor

    if(((ord+1) mod 10) eq 5 or ((ord+1) mod 10) eq 0 and ord ne nord-1) then $
      print,'MKSCATTER: Order '+strtrim(ord+1,2)+' of total ' $
           +strtrim(nord,2)+' was processed'
  endfor
  print,'MKSCATTER: Order '+strtrim(nord,2)+' was processed'

;Interpolate missing data from the adjacent troffs
  for i=0,ncol-1 do begin
    y=yback(i,*)
    k=where(back_data(i,*) gt -1000)
    back_data(i,*) = interpol(back_data(i,k), yback(i,k), yback(i,*))
  endfor

;Filter out noise: 0th background troff
  range = [colrange(0,0),colrange(1,0)]
  x = indgen(range(1) - range(0) + 1) + range(0)
  b = middle(back_data(x,0),20.,eps=dev)
  if(keyword_set(pol)) then $
    back(x,0) = middle(b,11,/poly,eps=dev,min=0) $
  else $
    back(x,0) = middle(b,lam_sp,eps=dev,min=0)

;All intermediate background troffs
  if(nord gt 1) then begin
    for ord=1,nord-1 do begin
      range = [colrange(0,ord-1)<colrange(0,ord), $
               colrange(1,ord-1)>colrange(1,ord)]
      x = indgen(range(1) - range(0) + 1) + range(0)
      b = middle(back_data(x,ord),20.,eps=dev)
      if(keyword_set(pol)) then $
        back(x,ord) = middle(b,11,/poly,eps=dev,min=0,/DOUBLE) $
      else $
        back(x,ord) = middle(b,lam_sp,eps=dev,min=0)
    endfor
  endif

;The last background troff
  range = [colrange(0,nord-1),colrange(1,nord-1)]
  x = indgen(range(1) - range(0) + 1) + range(0)
  b = middle(back_data(x,nord),20.,eps=dev)
  if(keyword_set(pol)) then $
    back(x,nord) = bottom(b,11,/poly,eps=dev,min=0,/DOUBLE) $
  else $
    back(x,nord) = bottom(b,lam_sp,eps=dev,min=0)

  if(keyword_set(subtract)) then begin
; Type conversion to avoid problems with UINT arrays when subtracting
    if(size(im, /TYPE) lt 4) then im = float(im)
    ycen = findgen(nrow)
    if(keyword_set(polariz)) then begin
      ii=indgen(nord/2)*2  ; polarization: orders come in pairs
      ii=[ii,max(ii)+2]    ; do not use interpolarization space
    endif else begin
      ii=indgen(nord+1)    ; no polarization, count all orders
    endelse
    for j = 0,ncol-1 do begin
      b = reform(back(j,*))  & b = b(ii)
      y = reform(yback(j,*)) & y = Y(ii)
      im(j,*) = im(j,*) -  interpol(b, y, ycen)
    endfor
  endif

  return
end

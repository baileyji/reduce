pro hamdord,im,orc,or_err=ome,or_range=orr,or_column=or_column $
           ,cross_disperser=cross,power=power,filter=filter $
           ,PLOT=iplot,POLARIM=polarim,MASK=mask,MANUAL=MANUAL $
           ,NOISE=noise,COLOR=color,THRES=thres
;Determines theorder locations for a particular spectrograph setting.
; Defaults are determined using the image (usually a narrow flat) in im.
; im (input 2D array) image where the orders are to be found. It should have
;   no overscan areas with orders preferentially align with the rows. Normally
;   it is a narrow flat with bias subtracted. The assumption about orientaton
;   is not critical for order detection: it is only used in polynomial fitting
;   part of the mark_orders. containing information specific
;   to a given observation.
; orc (output array [# coeffs , # orders]) coefficients from the
;   polynomial fits to the order locations.
; or_err (optional output vector (# orders)) each entry gives the mean of the
;   absolute value of the difference between order locations and the polynomial
;   fit to these locations.
; or_range (optional output vector int [# 2]) the number of the first and the
;   last detected order assuming that spacific type of cross-disperser. If the
;   number of detected orders is less than 3, orders are simply numbered from 0.
; cross_disperser (optional input string) contains the type of cross-disperser.
;   Default value is 'grating', the alternative is 'prism'. Anything that is
;   not a prism is interpreted as a grating. Case-insensitive. Ignored if the
;   number of detected orders is smaller that 3.
; power (optional input int) the power of the polynomial fit to the orders. The
;   default is 2.
; filter (optional input float) the with of the filter for detection of pixels
;   that belong to local maxima. The default is 20. Larger value detects less
;   pixels.
;
;Calls MARK_ORDERS
;03-Aug-2000 NP	This version is based on the new concept (no I/O insided
;               reduction subroutines) for REDUCE and on the new order
;               detection algorithm based on the cluster marking algorithms.

;@ham.common					;get common block definition

  if(n_params() lt 2 or (size(im))[0] ne 2) then begin
    print,'HAMDORD syntax:'
    print,'hamdord,im,orc[,or_err=ome[,or_range=orr[,cross_disperser=cross $'
    print,'              [,power=power],filter=filter[,/PLOT[,/POLARIM $'
    print,'              [,/MANUAL[,NOISE=noise]]]]]]]'
    print,'where: im  - 2D image with orders oriented along the rows'
    print,'       orc - 2D array of polynomial coefficients for detected orders'
    print,'       ome - RMS for each order'
    print,'       orr - 2-element int array with estimated spectrometer order range'
    print,'       cross - type of cross-disperser ''grating'' or ''prism''. Default is ''grating'''
    print,'       power - power of the polynomial for order fitting. Default is 2'
    print,'       filter - filter width for detection of orders. Default: 20.'
    print,'       /PLOT  - plot the process, also switches to the interactive mode'
    print,'       /POLARIM - assumes polarimetric mode with a double copy of each order'
    print,'       /MANUAL - all cluster merging is done in interactive mode'
    print,'       noise   - readout noise of the CCD may be important if the S/N'
    print,'                 in the image used is too low'
    print,'       /COLOR  - use colors when marking clusters.'
  endif

;Determine order locations
  if(not keyword_set(filter)) then filter=20.
  mark_orders,im,orc,POWER=power,FILTER=filter,ERROR=ome $
             ,PLOT=iplot,MASK=mask,MANUAL=MANUAL,NOISE=noise $
             ,COLOR=color,CLUSTER_LIMITS=cl_limits,THRES=thres

  nord=n_elements(orc[0,*])
  ncol=n_elements(im[*,0])
  nrow=n_elements(im[0,*])

  if(nord lt 3) then begin
;We have only one or two orders.
    orr=[0,(nord-1)>0]
    return
  endif else begin
;Find the useful column range for each order
    x=dindgen(ncol)
    or_column=intarr(2,nord)
;
;Some of the orders can be partially cut by the edge
;of the detector. Now we will try to determine the extent of each order which
;is useful for extraction. We will do it in 2 steps:
    ibad=-1
    for iord=0,nord-1 do begin
; Step 1: find middle lines between the current and the two adjacent orders
      if(iord eq 0) then begin
; Extrapolate below the first order
        y1=poly(x,orc[*,0])
        y2=poly(x,orc[*,1])
	x1=cl_limits(0,0)>cl_limits(0,1)
        x2=cl_limits(1,0)<cl_limits(1,1)
        if(x1 gt x2) then x1=x2
        deltay=ceil(min(y2(x1:x2)-y1(x1:x2))*0.5)
;        deltay=ceil(min(y2(x1:x2)-y1(x1:x2))*0.75)
        ym1=y1-deltay
        ym2=y1+deltay
      endif else if(iord eq nord-1) then begin
; Extrapolate above the last order
        y1=poly(x,orc[*,nord-2])
        y2=poly(x,orc[*,nord-1])
	x1=cl_limits(0,nord-2)>cl_limits(0,nord-1)
        x2=cl_limits(1,nord-2)<cl_limits(1,nord-1)
        if(x1 gt x2) then x1=x2
        deltay=ceil(min(y2(x1:x2)-y1(x1:x2))*0.5)
;        deltay=ceil(min(y2(x1:x2)-y1(x1:x2))*0.75)
        ym1=y2-deltay
        ym2=y2+deltay
      endif else begin
; For intermediate orders take the middle points
        y1=poly(x,orc[*,iord-1])
        y2=poly(x,orc[*,iord  ])
	x1=cl_limits(0,iord)>cl_limits(0,iord-1)
        x2=cl_limits(1,iord)<cl_limits(1,iord-1)
        if(x1 gt x2) then x1=x2
        deltay=ceil(min(y2(x1:x2)-y1(x1:x2))*0.5)
;        deltay=ceil(min(y2(x1:x2)-y1(x1:x2))*0.75)
        ym1=y2-deltay
        y1=y2
        y2=poly(x,orc[*,iord+1])
	x1=cl_limits(0,iord)>cl_limits(0,iord+1)
        x2=cl_limits(1,iord)<cl_limits(1,iord+1)
        if(x1 gt x2) then x1=x2
        deltay=ceil(min(y2(x1:x2)-y1(x1:x2))*0.5)
;        deltay=ceil(min(y2(x1:x2)-y1(x1:x2))*0.75)
        ym2=y1+deltay
      endelse

; Step 2: Find the range of columns where both ym1 and ym2 are
; within the image
      i=where(ym1(x1:x2) gt 0 and ym2(x1:x2) lt nrow, ni)
      or_column[*,iord]=[min(i),max(i)]+x1
; We can get one of the five cases: 1) all in range
;                                   2) left end in range
;                                   3) right end in range
;                                   4) part in the middle is in range
;                                   5) both ends are in range, but middle part is outside
; We choose to reject orders where less than 20% of columns is in range
; and all orders that belong to case 5.
      if(ni lt 0.2*ncol or (max(i)-min(i)+1 gt ni)) then ibad=[ibad,iord]
    endfor

    nbad=n_elements(ibad)
    if(nbad gt 1) then begin
      nbad=nbad-1
      ibad=ibad[1:nbad]
      igood=indgen(nord)
      igood[ibad]=-1
      igood=where(igood ge 0, ngood)
      if(ngood lt 1) then begin
        print,'HAMDORD: Not a single order was found in the image'
        stop
      endif
      orc=orc[*,igood]
      ome=ome[igood]
;      if(orr[0] lt orr[1]) then begin
;        orr[0]=orr[0]+min(igood)
;        orr[1]=orr[1]-(nord-max(igood)-1)
;      endif else begin
;        orr[0]=orr[0]-min(igood)
;        orr[1]=orr[1]+(nord-max(igood)-1)
;      endelse
      or_column=or_column[*,igood]
      nord=ngood
    endif
;
;We have at least 3 useful orders, we can try to find where we are in
;the focal plane by comparing order spacings.
    if(not keyword_set(cross)) then cross='grating'
    cross=strlowcase(cross)
    if(cross eq 'prism') then begin
      print,'HAMDORD: Well, I would need a specific expression for linear'
      print,'         dispersion as a function of the wavelength before'
      print,'         I can estimate order numbers for a prism.'
      orr=[0,nord-1]
    endif else begin
;
;We have a grating cross-disperser. That means:
; delta_y * order_number * (order_number - 1) = const.
;Therefore, we will check all order numbers from 1 to 1000 and look
;for minimum relative deviation from the mean
;      y=dblarr(nord)
;      for i=0,nord-1 do y[i]=poly(ncol*0.5d0,orc[*,i])
;      dy=y[1:nord-1]-y[0:nord-2]
      x=dindgen(ncol)+10
      dy=dblarr(nord-1)
      for i=0,nord-2 do dy[i]=mean(poly(x,orc[*,i+1])-poly(x,orc[*,i]))

      max_sp_order=500
      x=dindgen(max_sp_order)+10
;
;We may have increasing spacing which indicates that
;order numbers are decreasing
      if(total(dy[1:nord-2]-dy[0:nord-3]) lt 0) then begin
        dx = dindgen(nord-1)-1
        dx0= 1
      endif else begin
        dx =nord-2-dindgen(nord-1)
        dx0=-1
      endelse

      y=dblarr(max_sp_order)
      for i=0,max_sp_order-1 do begin
        yy=dy*(x[i]+dx)*(x[i]+dx+dx0)
;        yyavr=total(yy)/(nord-1)
;        y[i]=sqrt(total((yy-yyavr)^2))/yyavr
        y[i]=total(abs(yy(1:nord-2)-yy(0:nord-3)))
      endfor

      i=min(y,m)
      m=m+10
      if(dx0 gt 0) then orr=[m,m+nord-1] $
      else              orr=[m,m-nord+1]
    endelse

    if(keyword_set(iplot)) then begin
      display,im,/LOG
      x=dindgen(ncol)
      for iord=0,nord-1 do begin
        i1=or_column[0,iord]
        i2=or_column[1,iord]
        oplot,x[i1:i2],poly(x[i1:i2], orc[*,iord]),col=0
      endfor
    endif

    return
  endelse
end

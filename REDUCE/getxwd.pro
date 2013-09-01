pro getxwd,im,orc,xwd,sig,COLRANGE=colrange,GAUSS=gauss,DEBUG=debug
;Determines the fraction of each order to be mashed during spectral extraction.
; im (input array (# columns , # rows)) image to use in determining spectral
;   extraction width.
; orc (input array (# coeffs , # orders)) coefficients of polynomials
;   describing order locations.
; xwd (output scalar) fractional extraction width, i.e. percent of each order
;   to be mashed during spectral extraction.
; sig (output scalar) standard deviation estimate for the value of xwd
;   derived from mutiple order analysis. sig=0 for a single order.
; colrange (optional input integer 2 x Nord array) has the first and the last
;   useful column number for each order.
;03-Dec-89 JAV	Create.
;18-Aug-98 CMJ  Added sxwd as a return variable.
;26-Nov-98 JAV  Use pmax-(pmin<0) to locate signal when pmin le 0.
;14-Oct-00 NP   Added handling logic for partial orders. The useful range
;               of columns is passed via optional parameter colrange.

if n_params() lt 3 then begin
  print,'syntax: getxwd,im,orc,xwd[,sxwd[,colrange=colrange[,/gauss]].'
  retall
end

;Easily changed program parameters.
  soff = 10					;offset to edge of swath
  pkfrac = 0.9					;allowable fraction of peak

;Define useful quantities.
  ncol = n_elements(im(*,0))			;# columns in image
  nrow = n_elements(im(0,*))			;# rows in image
  ndeg = n_elements(orc(*,0))-1			;degree of poly fit to orders
  nord = n_elements(orc(0,*))			;# orders in orc

  if(not keyword_set(colrange)) then begin
    colrange=intarr(2,nord)
    colrange[1,*]=ncol-1
  endif

;Calculate from orc the location of order peaks in center of image.
  if nord eq 1 then begin
    xwd=fltarr(2,1)
    x1=colrange[0,0]
    x2=colrange[1,0]
    x=0.5*(x1+x2)
    y1=0.d0
    y2=floor(poly(x,orc[*,0]))
    y3=nrow
    prof=total(im[x-soff:x+soff,0:nrow-1],1)
    yyy=dindgen(nrow)
    if(keyword_set(gauss)) then begin
      pg=disp_gaussfit(yyy,prof,ag,yyy)
      kgood=where(abs(pg-prof) lt 0.33*ag(0), ngood)
      if(ngood lt n_elements(prof)) then begin
        pg=disp_gaussfit(yyy(kgood),prof(kgood),ag,yyy)
      endif
      yym1=floor(ag(1)-ag(2)-2)>0
      yym2=ceil(ag(1)+ag(2)+2)<(nrow-1)
      xwd[0,0]=pkfrac*(y1-yym1+1.)/(y1-y2+1);fraction of order below central line
      xwd[1,0]=pkfrac*(yym2-y1+1.)/(y3-y1+1);fraction of order above central line
    endif else begin
      pmin=min(prof)				;background trough counts
      pmax=max(prof)				;order peak counts
      keep=prof(where(prof gt pmin+sqrt(pmax-(pmin<0)),nkeep))
      xwd[0,0]=pkfrac*(0.5+0.5*nkeep) / n_elements(prof)	;fraction of order to extract
      xwd[1,0]=pkfrac*(0.5+0.5*nkeep) / n_elements(prof)	;fraction of order to extract
    endelse
    sig=0.
  endif else begin
    xwd = fltarr(2,nord)			;extraction widths at x
    for iord = 0,nord-1 do begin
      x1=colrange[0,iord]
      x2=colrange[1,iord]
      x=0.5*(x1+x2)
      if(iord eq 0) then begin
        y1=poly(x,orc[*,0])
        y3=poly(x,orc[*,1])
	y2=y1-(y3-y1)
        ym1=ceil(y1-0.5*(y3-y1))>0
        ym2=ceil(0.5*(y1+y3))<(nrow-1)
      endif else if(iord eq nord-1) then begin
        y1=poly(x,orc[*,nord-1])
        y2=poly(x,orc[*,nord-2])
	y3=y1+(y1-y2)
        ym1=ceil(0.5*(y1+y2))>0
        ym2=ceil(y1+0.5*(y1-y2))<(nrow-1)
      endif else begin
        y1=round(poly(x,orc[*,iord  ]))
        y2=floor(poly(x,orc[*,iord-1]))
        y3=ceil(poly(x,orc[*,iord+1]))
        prof=total(im[x-soff:x+soff,y2:y3],1)
;	s=min(prof(0:y1-y2),ym1)       & ym1=(ym1+y2)>0
;	s=min(prof(y1-y2+1:y3-y2),ym2) & ym2=(ym2+y1)<(nrow-1)
        ym1=ceil(0.5*(y1+y2))>0
        ym2=ceil(0.5*(y1+y3))<(nrow-1)
      endelse
      yym1=ym1
      yym2=ym2
      yyy=dindgen(ym2-ym1+1)
      prof=total(im[x-soff:x+soff,ym1:ym2],1)
      if(keyword_set(gauss)) then begin
        pg=disp_gaussfit(yyy,prof,ag,yyy)
	kgood=where(abs(pg-prof) lt 0.33*ag(0), ngood)
	if(ngood lt n_elements(prof)) then begin
	  pg=disp_gaussfit(yyy(kgood),prof(kgood),ag,yyy)
	endif
        yym1=floor(ag(1)-ag(2)+ym1-2)>0
        yym2=ceil(ag(1)+ag(2)+ym1+2)<(nrow-1)
        xwd[0,iord]=pkfrac*(y1-yym1+1.)/(y1-y2+1);fraction of order below central line
        xwd[1,iord]=pkfrac*(yym2-y1+1.)/(y3-y1+1);fraction of order above central line
      endif else begin
        pmin=min(prof)				;background trough counts
        pmax=max(prof)				;order peak counts
        keep=prof(where(prof gt pmin+sqrt(pmax-(pmin<0)),nkeep))
        xwd[0,iord]=pkfrac*(0.5+0.5*nkeep) / n_elements(prof)	;fraction of order to extract
        xwd[1,iord]=pkfrac*(0.5+0.5*nkeep) / n_elements(prof)	;fraction of order to extract
;      if(y1-xwd[0,iord]*(ym2-ym1+1) lt ym1) then xwd[0,iord]=(y1-ym1)/(ym2-ym1+1)
;      if(y1+xwd[1,iord]*(ym2-ym1+1) gt ym2) then xwd[1,iord]=(ym2-y1)/(ym2-ym1+1)
      endelse
      if(keyword_set(debug)) then begin
        plot,yyy+ym1,prof,tit=strtrim(iord,2),xs=3,xr=[yym1<ym1,yym2>ym2]
        oplot,[ym1,ym1],!y.crange
        oplot,[ym2,ym2],!y.crange
        oplot,[y1,y1],!y.crange,line=2
        oplot,[yym1,yym1],!y.crange,line=1
        oplot,[yym2,yym2],!y.crange,line=1
        s=get_kbrd(1)
      endif
    endfor
;    sig = stdev(vxwd,xwd)			;standard deviation
     sig=0.1*max(xwd)
  endelse
  print,'GETXWD: Extraction width (min,max) = ' $
    + strtrim(string(min(xwd),max(xwd),form='(f10.3,'','',f10.3)'),2)

  print,'GETXWD: Sigma = ' $
    + strtrim(string(sig,form='(f10.3)'),2)

  return
end

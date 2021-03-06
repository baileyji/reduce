pro Gaussbox,x,a,f,pder  ; f=gaussbroad(xbox,box,a(3))+a(4)+a(5)*x
                         ;                0      x<a(1)-a(2)
                         ; where box(x) = a(0)  -a(2)>=x-a(1)<=a(2)
                         ;                0      x>a(1)+a(2)

  n=((n_elements(x)-1)*10L+1L)<1000L
  delta=double(max(x)-min(x))/(n-1.d0)
  xx=dindgen(n)*delta+min(x)
  yy=xx*0.d0
  jj=round(lindgen(n_elements(x))*(n-1.D0)/(n_elements(x)-1.D0))
  a(1)=(a(1)>(min(x)+1))<(max(x)-1)
  a(2)=abs(a(2))>delta
  ii=where(xx ge a(1)-a(2) and xx le a(1)+a(2), nii)
  if(nii le 0) then begin
    print,'GAUSSBOX: That should never happen 0'
    stop
  endif
  if(nii gt 0) then yy(ii)=1.d0

  a(3)=abs(a(3))
  g=(gaussbroad(xx,yy,a(3)))(jj)

  f=a(0)*g+a(4)+a(5)*x

  if(n_params(0) le 3) then return   ; Do not need partial derivatives?
  pder=dblarr(n_elements(x),6)       ; MAKE ARRAY.

  pder(*,0)=g                        ; COMPUTE PARTIALS

  yyy=xx*0.d0                        ; Right derivative over a(1)
  iii=where(xx ge a(1)-a(2)+delta and xx le a(1)+a(2)+delta, niii)
  if(niii le 0) then begin
    print,'GAUSSBOX: That should never happen 1'
    stop
  endif
  yyy(iii)=1.d0
  dgda=((gaussbroad(xx,yyy,a(3)))(jj)-g)*a(0)/delta
  yyy=xx*0.d0                        ; Left derivative over a(1)
  iii=where(xx ge a(1)-a(2)-delta and xx le a(1)+a(2)-delta, niii)
  if(niii le 0) then begin
    print,'GAUSSBOX: That should never happen 2'
    stop
  endif
  if(niii gt 0) then yyy(iii)=1.d0
  dgda=((g-(gaussbroad(xx,yyy,a(3)))(jj))*a(0)/delta+dgda)*0.5d0
  pder(*,1)=dgda                     ; Save derivative over a(1)

  yyy=xx*0.d0                        ; Right derivative over a(2)
  iii=where(xx ge a(1)-a(2)-delta and xx le a(1)+a(2)+delta, niii)
  if(niii le 0) then begin
    print,'GAUSSBOX: That should never happen 3'
    stop
  endif
  if(niii gt 0) then yyy(iii)=1.d0
  dgda=((gaussbroad(xx,yyy,a(3)))(jj)-g)*a(0)/delta
  if(a(2) gt 0.01d0) then begin
    yyy=xx*0.d0                      ; Left derivative over a(2) if possible
    iii=where(xx ge a(1)-a(2)+delta and xx le a(1)+a(2)-delta, niii)
    if(niii le 0) then begin
      print,'GAUSSBOX: That should never happen 4'
      stop
    endif
    if(niii gt 0) then yyy(iii)=1.d0
    dgda=((g-(gaussbroad(xx,yyy,a(3)))(jj))*a(0)/delta+dgda)*0.5d0
  endif
  pder(*,2)=dgda                      ; Save the derivative over a(2)

  dgda=((gaussbroad(xx,yy,a(3)*1.001d0))(jj)-g)*a(0)/(1.d-3*a(3))
  dgda=(g-(gaussbroad(xx,yy,a(3)*0.999d0))(jj))*a(0)/(1.d-3*a(3))+dgda
  pder(*,3)=0.5d0*dgda

  pder(*,4)=1.d0

  pder(*,5)=x
  return
end

Function Gaussboxfit,x,y,a,xx

  on_error,2                    ;Return to caller if an error occurs
  n = n_elements(y)             ;# of points.
  c = poly_fit(x,y,1,yf)        ;fit a straight line.
  yd = y-yf                     ;difference.
  del = sqrt(total(yd*yd))/n
  ii = where(abs(yd) lt 3*del, nii)
  if(nii gt 6) then begin
    c = poly_fit(x(ii),y(ii),1) & yf = poly(x,c)
    yd = y-yf
  endif

  ymax=max(yd,imax) & xmax=x(imax)   ;x,y and subscript of extrema
  ymin=min(yd,imin) & xmin=x(imin)
  ii=(reverse(sort(yd)))(0:4)
  xmax=total(x(ii)*yd(ii))/total(yd(ii))
  a=dblarr(6)                   ;coefficient vector
  i0=imax                       ;emission
  i0 = (i0 > 1) < (n-2)         ;never take edges
  dy=yd(i0)                     ;diff between extreme and mean
  ii=where(yd gt 0.5*ymax)
  del=(0.5*(max(x(ii))-min(x(ii))))>1.d0
  a = double([yd(i0), xmax, del, del*0.1d0, c(0), c(1)]) ;estimates
  yy=mpcurvefit(x,y,replicate(1.d0,n),a,sigmaa, $
              function_name = 'GAUSSBOX',ITMAX=20,NODERIVATIVE=0) ;call curvefit

  gaussbox,xx,a,yy               ;compute the best fit on XX-grid
  xmax=a(1)                      ;coordinate of the center
  imax=max(where(xx lt xmax))    ;roughly half of the profile is in [0:imax]
  hm=a(0)*0.5+a(4)+a(5)*xmax
  ii=(sort(abs(hm-yy(0:imax))))(0:1)
  xleft=(hm-yy(ii(0)))/(yy(ii(1))-yy(ii(0)))*(xx(ii(1))-xx(ii(0)))+$
         xx(ii(0))
  ii=(sort(abs(hm-yy(imax+1:n_elements(yy)-1))))(0:1)+imax+1
  xright=(hm-yy(ii(0)))/(yy(ii(1))-yy(ii(0)))*(xx(ii(1))-xx(ii(0)))+$
          xx(ii(0))
  a(2)=abs(xright-xleft)
  return,yy
end

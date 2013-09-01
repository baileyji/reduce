pro lorentz,x,a,f,pder ; f=a(0)/(a(2)^2+z^2)+a(3)+a(4)*x
                       ; where: z=(x-a(1))

  z=x-a(1)
  f=a(0)/(a(2)^2+z^2)+a(3)+a(4)*x
  if(n_params(0) le 3) then return   ; Do not need partial derivatives?
  pder=dblarr(n_elements(x),5)       ; MAKE ARRAY.
  pder(*,0)= 1/(a(2)^2+z^2)          ; COMPUTE PARTIALS
  pder(*,1)= a(0)/(a(2)^2+z^2)^2*2*z
  pder(*,2)=-a(0)/(a(2)^2+z^2)^2*2*a(2)
  pder(*,3)= 1.d0
  pder(*,4)= x
  return
end

Function Lorentzfit,x,y,a,xx
  on_error,2                    ;Return to caller if an error occurs
  n = n_elements(y)     ;# of points.
  c = poly_fit(x,y,1,yf)    ;fit a straight line.
  yd = y-yf         ;difference.
  del = sqrt(total(yd*yd))/n
  ii = where(abs(yd) lt 3*del, nii)
  if(nii gt 6) then begin
    c = poly_fit(x(ii),y(ii),1) & yf = poly(x,c)
    yd = y-yf
  endif

  ymax=max(yd) & xmax=x(!c) & imax=!c   ;x,y and subscript of extrema
  ymin=min(yd) & xmin=x(!c) & imin=!c
  a=fltarr(5)           ;coefficient vector
  i0=imax                   ;emission
  i0 = i0 > 1 < (n-2)       ;never take edges
  dy=yd(i0)         ;diff between extreme and mean
  del = dy/exp(1.)      ;1/e value
  i=0
  while((i0+i+1) lt n) and $    ;guess at 1/2 width.
    ((i0-i) gt 0) and $
    (abs(yd(i0+i)) gt abs(del)) and $
    (abs(yd(i0-i)) gt abs(del)) do i=i+1
  a = double([yd(i0), x(i0), abs(x(i0)-x(i0+i)), c(0), c(1)]) ;estimates
  !c=0              ;reset cursor for plotting
  yy=curvefit(x,y,replicate(1.d0,n),a,sigmaa, $
          function_name = 'LORENTZ',ITMAX=20,NODERIVATIVE=0) ;call curvefit
  lorentz,xx,a,yy
  a(2)=2.d0*a(2)
  return,yy
end

pro gauss,x,a,f,pder  ; f=a(0)*exp(-z^2/2)+a(3)+a(4)*x
                      ; where: z=(x-a(1))/a(2)

  if a(2) ne 0.d0 then Z = (X-A(1))/A(2) $  ;GET Z
  else z= 10.D0
  ez=exp(-z^2/2.d0)*(abs(z) le 7.d0) ; GAUSSIAN PART IGNORE SMALL TERMS
  f=a(0)*ez+a(3)+a(4)*x              ; FUNCTIONS.
  if(n_params(0) le 3) then return   ; Do not need partial derivatives?
  pder=dblarr(n_elements(x),5)       ; MAKE ARRAY.
  pder(*,0)=ez                   ; COMPUTE PARTIALS
  if(a(2) ne 0.d0) then pder(*,1)=a(0)*ez*z/a(2)
  pder(*,2)=pder(*,1)*z
  pder(*,3)=1.d0
  pder(*,4)=x
  return
end

Function disp_gaussfit,x,y,a,xx
  on_error,2                    ;Return to caller if an error occurs
  n = n_elements(y)             ;# of points.
  iy=(sort(y))(0:3)
  c = poly_fit(x(iy),y(iy),1,/DOUBLE);fit a straight line.
  yf = poly(x,c)
  yd = y-yf                     ;difference.
;  del = sqrt(total(yd*yd))/n
;  ii = where(abs(yd) lt 3*del, nii)
;  if(nii gt 6) then begin
;    c = poly_fit(x(ii),y(ii),1,/DOUBLE) & 
;    yf = poly(x,c)
;    yd = y-yf
;  endif

  ymax=max(yd) & xmax=x(!c) & imax=!c   ;x,y and subscript of extrema
  ymin=min(yd) & xmin=x(!c) & imin=!c
  a=dblarr(5)           ;coefficient vector
  i0=imax               ;emission
  i0 = i0 > 1 < (n-2)   ;never take edges
  dy=yd(i0)             ;diff between extreme and mean
  del = dy/exp(1.D0)    ;1/e value
  i=0
  while((i0+i+1) lt n) and $    ;guess at 1/2 width.
    ((i0-i) gt 0) and $
    (abs(yd(i0+i)) gt abs(del)) and $
    (abs(yd(i0-i)) gt abs(del)) do i=i+1
  a = double([yd(i0), x(i0), abs(x(i0)-x(i0+i)), c(0), c(1)]) ;estimates
;  !c=0                                 ;reset cursor for plotting
  yy=mpcurvefit(x,y,replicate(1.d0,n),a,sigmaa, $
          function_name = 'GAUSS',ITMAX=20,NODERIVATIVE=0,/QUIET) ;call curvefit
  gauss,xx,a,yy
  a(2)=2.d0*sqrt(2.d0*alog(2.d0))*abs(a(2))
  return,yy
end

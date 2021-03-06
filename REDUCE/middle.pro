Function middle,f,filter,ITER=iter,EPS=eps,POLY=pol $
               ,MIN=mn1,MAX=mx1,LAM2=lam2,WEIGHT=weight
; middle tries to fit a smooth curve that is located
; along the "middle" of 1D data array f. Filter size "filter"
; together with the total number of iterations determine
; the smoothness and the quality of the fit. The total
; number of iterations can be controlled by limiting the
; maximum number of iterations (iter) and/or by setting
; the convergence criterion for the fit (eps)
; 04-Nov-2000 N.Piskunov wrote.
; 09-Nov-2011 NP added weights and 2nd derivative constraint as LAM2
;
  if(n_params() lt 2) then begin
    print,'syntax: middle,f,{filter/order}[,ITER=iter[,EPS=eps $'
    print,'             [,MIN=mn[,[MAX=mx[,/POLY[,LAM2=lam2[,WEIGHT=wgt]]]]]]]]'
    print,'where f      is the function to fit,'
    print,'      filter is the smoothing parameter for the optimal filter.'
    print,'             If POLY is set, it is interpreted as the order'
    print,'             of the smoothing polynomial,'
    print,'      iter   is the maximum number of iterations [def: 40]'
    print,'      eps    is convergence level [def: 0.001]'
    print,'      mn     minimum function values to be considered [def: min(f)]'
    print,'      mx     maximum function values to be considered [def: max(f)]'
    print,'      lam2   constraint on 2nd derivative'
    print,'      wgt    vector of weights.'
    return,0
  endif
  if(not keyword_set(iter)) then iter=40
  if(not keyword_set(eps)) then eps=0.001
  if(not keyword_set(mn1)) then mn=min(f) else mn=mn1
  if(not keyword_set(mx1)) then mx=max(f) else mx=mx1
  if(keyword_set(pol)) then begin
    j=where(f ge mn and f le mx, n)
    if(n le round(filter)) then return,f
    xx=(2*dindgen(n_elements(f))/(n_elements(f)-1)-1)
    fmin=min(f(j))-1
    fmax=max(f(j))+1
    ff=(f(j)-fmin)/(fmax-fmin)
    ff_old=ff
  endif else begin
    fmin=min(f)-1
    fmax=max(f)+1
    ff=(f-fmin)/(fmax-fmin)
    ff_old=ff
    n=n_elements(f)
  endelse
  if(keyword_set(weight)) then wgt=weight else wgt=1

  i=0
next:i=i+1
  if(keyword_set(pol)) then begin
    order=round(filter)
    if(order gt 0) then begin ; This is a bug in RSI poly routine
      t=median(poly(xx,poly_fit(xx,ff,order,/DOUBLE)),3)
      dev=sqrt(poly(xx,poly_fit(xx,(t-ff)^2,round(filter),/DOUBLE))>0.)
    endif else begin
      t=replicate(poly_fit(xx,ff,order,/DOUBLE),n_elements(f))
      dev=sqrt(replicate(poly_fit(xx,(t-ff)^2,order,/DOUBLE),n_elements(f))>0.)
    endelse
  endif else begin
    t=median(opt_filter(ff,filter,WEIGHT=weight,LAM2=lam2),3)
    dev=sqrt(opt_filter(wgt*(t-ff)^2,filter,WEIGHT=weight,LAM2=lam2)>0.d0)
  endelse
  ff=((t-dev)>ff)<(t+dev)
  dev2=max(wgt*abs(ff-ff_old))
  ff_old=ff
  if(dev2 gt eps and i lt iter) then goto,next

  if(keyword_set(pol)) then begin
    x=(2*dindgen(n_elements(f))/(n_elements(f)-1)-1)
    if(order gt 0) then begin ; This is a bug in RSI poly routine
      t=median(poly(xx,poly_fit(xx,ff,order,/DOUBLE)),3)
    endif else begin
      t=replicate(poly_fit(xx,ff,order,/DOUBLE),n_elements(f))
    endelse
  endif

  return,t*(fmax-fmin)+fmin
end

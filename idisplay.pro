Pro idisplay,im,xx,yy,XRANGE=xxr,YRANGE=yyr,LOG=lg,HIST=hist $
            ,MIN=mn,MAX=mx,title=tit,ASPECT=aspect

  if(n_params() lt 1) then goto,syntax

  sz=size(im)
  if(sz(0) ne 2) then begin
    print,'IDISPLAY: Image must be a 2D array'
    goto,syntax
  endif

  nx=sz(1)
  ny=sz(2)
  if(keyword_set(xx)) then x=xx else x=findgen(nx)
  if(keyword_set(yy)) then y=yy else y=findgen(ny)

  if(nx ne n_elements(x)) then  begin
    print,'IDISPLAY: The column coordinates do not match the number of columns'
    goto,syntax
  endif
  if(ny ne n_elements(y)) then begin
    print,'IDISPLAY: The row coordinates do not match the number of rows'
    goto,syntax
  endif
  
  if(keyword_set(mn)) then min_val=mn else min_val=min(im)
  if(keyword_set(mx)) then max_val=mx>min_val else max_val=max(im)>min_val

  imm=(im>min_val)<max_val
  if(keyword_set(lg)) then begin
    step=(max_val-min_val)/256.
    h=histogram(imm,min=min_val,max=max_val,bin=step)
    imm=bytscl(alog10(imm-min_val+step))
  endif else if(keyword_set(hist)) then begin
    imm=hist_equal(imm,min=min_val,max=max_val,bin=(max_val-min_val)/256.)
  endif else begin
    imm=bytscl(imm,min=min_val,max=max_val)
  endelse

  if(keyword_set(xxr)) then begin
    xr=xxr
    i=where(x ge xr(0) and x le xr(1), n)
    if(n gt 0) then begin
      x=x(i)
      imm=imm(i,*)
      nx=n
    endif
  endif else xr=minmax(x)

  if(keyword_set(yyr)) then begin
    yr=yyr
    i=where(y ge yr(0) and y le yr(1), n)
    if(n gt 0) then begin
      y=y(i)
      imm=imm(*,i)
      ny=n
    endif
  endif else yr=minmax(y)

  pmulti0=!p.multi(0)
  plot,minmax(x),minmax(y),xs=1,ys=1,xr=xr,xtick_get=xtck $
                          ,yr=yr,ytick_get=ytck,/nodata $
                          ,title=tit,iso=aspect
  pmulti1=!p.multi(0)
  !p.multi(0)=pmulti0

  if((!d.flags mod 2) eq 0) then begin ; Display
    x_pix=!x.window*!d.x_vsize & x_pix=[ceil(x_pix(0)),floor(x_pix(1))]
    dx=0.5*(x(2:nx-1)+x(1:nx-2))-0.5*(x(1:nx-2)+x(0:nx-3))
    dx=[(x(1)-x(0)),dx,(x(nx-1)-x(nx-2))]
    x_scl=(max(x_pix)-min(x_pix))/(max(x)-min(x)+0.5*(dx(0)+dx(nx-1)))
    x0=[x(0)-0.5*dx(0),0.5*(x(1:nx-1)+x(0:nx-2)),x(nx-1)+0.5*dx(nx-1)]
    x0=(x0-x0(0))*x_scl+x_pix(0)
    dx=round(dx*x_scl+0.5)>1

    y_pix=!y.window*!d.y_vsize & y_pix=[ceil(y_pix(0)),floor(y_pix(1))]
    dy=0.5*(y(2:ny-1)+y(1:ny-2))-0.5*(y(1:ny-2)+y(0:ny-3))
    dy=[(y(1)-y(0)),dy,(y(ny-1)-y(ny-2))]
    y_scl=(max(y_pix)-min(y_pix))/(max(y)-min(y)+0.5*(dy(0)+dy(ny-1)))
    y0=[y(0)-0.5*dy(0),0.5*(y(1:ny-1)+y(0:ny-2)),y(ny-1)+0.5*dy(ny-1)]
    y0=(y0-y0(0))*y_scl+y_pix(0)
    dy=round(dy*y_scl+0.5)>1
    for j=0,ny-1 do begin
      for i=0,nx-1 do begin
        tv,replicate(imm(i,j),dx[i],dy[j]),x0(i),y0(j)
      endfor
    endfor
  endif else begin ; Scalable pixels
    x_siz=!d.x_vsize/!d.x_px_cm*(!x.crange(1)-!x.crange(0))*!x.s(1)
    x_org=!d.x_vsize/!d.x_px_cm*(!x.s(1)*!x.crange(0)+!x.s(0))
;    x=!x.s(1)*x+!x.s(0)
    x0=x_org
    dx=0.5*(x(2:nx-1)+x(1:nx-2))-0.5*(x(1:nx-2)+x(0:nx-3))
    dx=[(x(1)-x(0)),dx,(x(nx-1)-x(nx-2))]
    x_scl=x_siz/(max(x)-min(x)+0.5*(dx(0)+dx(nx-1)))
;    x0=[x(0)-0.5*dx(0),0.5*(x(1:nx-1)+x(0:nx-2)),x(nx-1)+0.5*dx(nx-1)]
;    x0=(x0-x0(0))*x_scl+x_org
    dx=dx*x_scl

    y_siz=!d.y_vsize/!d.y_px_cm*(!y.crange(1)-!y.crange(0))*!y.s(1)
    y_org=!d.y_vsize/!d.y_px_cm*(!y.s(1)*!y.crange(0)+!y.s(0))
;    y=!y.s(1)*y+!y.s(0)
    dy=0.5*(y(2:ny-1)+y(1:ny-2))-0.5*(y(1:ny-2)+y(0:ny-3))
    dy=[(y(1)-y(0)),dy,(y(ny-1)-y(ny-2))]
    y_scl=y_siz/(max(y)-min(y)+0.5*(dy(0)+dy(ny-1)))
;    y0=dblarr(ny)
;    for j=1,ny-1 do y0(j)=total(dy(0:j-1))
;    y0=[y(0)-0.5*dy(0),0.5*(y(1:ny-1)+y(0:ny-2)),y(ny-1)+0.5*dy(ny-1)]
;    y0=y0-y0(0)
;    y0=y0*y_scl+y_org
    dy=dy*y_scl
;    nxx=(((x(nx-1)-x(0))/min(dx))<10000)+1
;    xxx=dindgen(nx)/(nx-1)*(x(nx-1)-x(0))
    y0=y_org

; Rebinning scale rbn=10. - mm
;                    = 1. - cm
    rbn=10.

    xr=[x(0)-0.5*dx(0),x(nx-1)+0.5*dx(nx-1)]
    x(0)=x0;-0.5*dx(0)
    for i=1,nx-1 do x(i)=x0+total(dx(0:i-1))
    i=uniq(fix(x*rbn))
    nx=n_elements(i)

    i0=0
    xxx=dblarr(nx)
    dxx=dblarr(nx)
    immm=bytarr(nx,ny)
    for ii=0,nx-1 do begin
      xxx(ii)=x(i0)
;      dxx(ii)=x(i(ii)+1)-xxx(ii)
      dxx(ii)=total(dx(i0:i(ii)))
      if(i(ii) gt i0) then $
        immm(ii,*)=byte(total(imm(i0:i(ii),*),1)/(i(ii)-i0+1)) $
      else $
        immm(ii,*)=imm(i(ii),*)
      i0=i(ii)+1
    endfor
    x=xxx
    dx=dxx
    if(dx(0) gt 0) then dx=dx+0.001d0 else dx=dx-0.001d0

    yr=[y(0)-0.5*dy(0),y(ny-1)+0.5*dy(ny-1)]
    y(0)=y0;-0.5*dy(0)
    for j=1,ny-1 do y(j)=y0+total(dy(0:j-1))
    j=uniq(fix(y*rbn))
    ny=n_elements(j)

    j0=0
    yyy=dblarr(ny)
    dyy=dblarr(ny)
    imm=bytarr(nx,ny)
    for jj=0,ny-1 do begin
      yyy(jj)=y(j0)
;      dyy(jj)=y(j(jj)+1)-yyy(jj)
      dyy(jj)=total(dy(j0:j(jj)))
      if(j(jj) gt j0) then $
        imm(*,jj)=byte(total(immm(*,j0:j(jj)),2)/(j(jj)-j0+1)) $
      else $
        imm(*,jj)=immm(*,j(jj))
      j0=j(jj)+1
    endfor
    y=yyy
    dy=dyy
    if(dy(0) gt 0) then dy=dy+0.001d0 else dy=dy-0.001d0
    immm=0

    for j=0,ny-1 do begin
      for i=0,nx-1 do begin
        tv,replicate(imm(i,j),2,2),x(i),y(j),xs=dx[i],ys=dy[j],/centimeter
      endfor
    endfor
;    stop
  endelse

  plot,minmax(x),minmax(y),xs=1,ys=1,xr=xr,yr=yr,/nodata,/noerase,iso=aspect
  !p.multi(0)=pmulti1
  return
  
syntax:
  print,'Usage: idisplay,im[,x[,y[,XRANGE=xr[,YRANGE=yr[,/LOG'
  print,'                  [,/HIST[,MIN=mn[,MAX=mx[,title=tit]]]]]]]]]'
  print,'where: im    - 2D array to be displayed'
  print,'       x     - 1D array with coordinates of the of the columns'
  print,'       y     - 1D array with coordinates of the of the rows'
  print,'       xr    - 2-element array with the range of columns to display'
  print,'       yr    - 2-element array with the range of rows to display'
  print,'       /LOG  - display image in log scale (takes precedence over /HIST)'
  print,'       /HIST - display histogram equalized image'
  print,'       mn    - minimum value to display'
  print,'       mx    - maximum value to display'
  print,'       tit   - plot title'
  return
end

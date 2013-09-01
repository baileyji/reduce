Pro disp_2d,m_wave,m_pix,m_ord,m_flag,maxordx,maxordy,npix,nord,obase,wvc $
           ,resid_pix=resid_pix,resid_m_s=resid_m_s

common thid_common, thid_order, thid_pixel, thid_maxordx, thid_maxordy

;Prepare line list for 2D fit
  nm = n_elements(m_wave)
  m_w = double(m_wave)
  m_p = double(m_pix)
  m_o = double(m_ord)
  m_f = double(m_flag)

  igd = where(m_flag eq 0, ngd)

  m_w = m_w(igd)
  m_p = m_p(igd)
  m_o = m_o(igd)
  m_f = m_f(igd)
  m_ml = m_o * m_w

;List of sp. orders in the correct sequence
  olist = obase + indgen(nord)      ;list of orders
  if(m_ord(nm-1) lt m_ord(0)) then olist = rotate(olist, 5)

;Determine what order polynomials can be constrained.
  ordxok = intarr(maxordx+1)
  for i=0, n_elements(ordxok)-1 do begin
    bsiz = npix / float(i+1)
    h = histogram(m_p, min=0, max=npix-1, bin=bsiz)
    if min(h) gt 0 then ordxok(i) = 1
  endfor
  nordxok = total(ordxok)
  ordyok = intarr(maxordy+1)
  for i=0, n_elements(ordyok)-1 do begin
    bsiz = nord / float(i+1)
    h = histogram(m_o, min=obase, max=obase+nord-1, bin=bsiz)
    if min(h) gt 0 then ordyok(i) = 1
  endfor
  nordyok = total(ordyok)

;Choose polynomial order in X.
  if nordxok eq 1 then begin
    iwhr = where(ordxok eq 1)
    ordx = iwhr(0)
    print, 'Polynomial order in X: ', strtrim(ordx, 2)
  endif else begin
    list = ''
    for i=0, nordxok-1 do begin
      if ordxok(i) eq 1 then begin
        list = list + ',' + strtrim(i, 2)
      endif
    endfor
    ordx=max(where(ordxok eq 1))
;    print, 'Polynomial order in X: ', strtrim(ordx, 2)
  endelse

;Choose polynomial order in Y.
  if nordyok eq 1 then begin
    iwhr = where(ordyok eq 1)
    ordy = iwhr(0)
;    print, 'Polynomial order in Y: ', strtrim(ordy, 2)
  endif else begin
    list = ''
    for i=0, nordyok-1 do begin
      if ordyok(i) eq 1 then begin
        list = list + ',' + strtrim(i, 2)
      endif
    endfor
    ordy=max(where(ordyok eq 1))
;    print, 'Polynomial order in Y: ', strtrim(ordy, 2)
  endelse

;Decide on cross-terms (use all that are OK).
  maxcross=6
  xok = intarr(maxcross)
  xlist = ''
  nxok = 0
  if ordx ge 1 and ordy ge 1 and ngd ge ordx+ordy+nxok+2 then begin
    xok(0) = 1
    xlist = xlist + 'XY '
    nxok = nxok + 1
  endif
  if ordx ge 2 and ordy ge 1 and ngd ge ordx+ordy+nxok+2 then begin
    xok(1) = 1
    xlist = xlist + 'XXY '
    nxok = nxok + 1
  endif
  if ordx ge 1 and ordy ge 2 and ngd ge ordx+ordy+nxok+2 then begin
    xok(2) = 1
    xlist = xlist + 'XYY '
    nxok = nxok + 1
  endif
  if ordx ge 2 and ordy ge 2 and ngd ge ordx+ordy+nxok+2 then begin
    xok(3) = 1
    xlist = xlist + 'XXYY '
    nxok = nxok + 1
  endif
  if ordx ge 3 and ordy ge 1 and ngd ge ordx+ordy+nxok+2 then begin
    xok(4) = 1
    xlist = xlist + 'XXXY '
    nxok = nxok + 1
  endif
  if ordx ge 1 and ordy ge 3 and ngd ge ordx+ordy+nxok+2 then begin
    xok(5) = 1
    xlist = xlist + 'XYYY '
    nxok = nxok + 1
  endif
  nxok = total(xok)
;  if nxok gt 0 then print, 'Cross terms: ' + xlist

;Load common block used by fitting routine.
  thid_order = m_o / 1d2
  thid_pixel = m_p / 1d3
  thid_maxordx = maxordx
  thid_maxordy = maxordy

;Use Marquardt to fit marked lines.
  dummy = fltarr(ngd)
  sig = replicate(1.0, ngd)
  par = dblarr(maxordx + 1 + maxordy + maxcross)
  par(0:1) = poly_fit(m_p/1d3, m_ml, 1, /double)
  dpar = [ 1d0, replicate(1.d-2, maxordx) $
              , replicate(1.d-2, maxordy) $
              , replicate(1.d-2, maxcross) ]
  if ordx lt maxordx then dpar(ordx+1:maxordx) = 0.0
  if ordy lt maxordy then dpar(ordy+maxordx+1:maxordx+maxordy) = 0.0
  iwhr = where(xok eq 0, nwhr)
  if nwhr gt 0 then dpar(maxordx+maxordy+1+iwhr) = 0.0
  mlfit = marq('thid_func', dummy, m_ml, sig, par, dpar, trace=0)

;Construct wavelength coefficients.
  vers = 2.5
  fill = 0.0
  wvc = [ vers, npix, nord, obase, fill, fill, fill $
        , maxcross, maxordx, maxordy, par ]

;Calculate wavelengths for every pixel.
  if(arg_present(resid_pix)) then begin
    wave=dblarr(npix,nord)
    dummy = fltarr(npix)
    thid_pixel = dindgen(npix) / 1d3
    for i=0, nord-1 do begin
      thid_order = replicate(olist(i), npix) / 1d2
      wave(*,i) = thid_func(dummy, par) / olist(i)
    endfor
    resid_pix = fltarr(nm)
    resid_m_s = fltarr(nm)
    xpix = dindgen(npix)
    for i=0,n_elements(olist)-1 do begin
      order=olist[i]
      iwhr = where(m_ord eq order, nwhr)
      if nwhr gt 0 then begin
        resid_pix(iwhr) = m_pix(iwhr)  - interpol(xpix, wave(*,order-obase), m_wave(iwhr))
        resid_m_s(iwhr) = m_wave(iwhr) - interpol(wave(*,order-obase), xpix, m_pix(iwhr))
      endif
    endfor
    resid_m_s=resid_m_s/m_wave*2.9979246d8
  endif
  return
end

;  f=findfile('\nikolai\reduction\elodie_order*.sav')
;  f=f(sort(f))
;  m_wave=0.d0
;  m_pix=0.d0
;  m_flag=0
;  m_ord=0.
;  npix=1024
;  ord=0
;  for i=0,n_elements(f)-1 do begin
;    order=fix(strmid(f[i],strlen(f[i])-7,3))
;    restore,f[i]
;    nlines=n_elements(line)
;    ord=[ord,order]
;    m_wave=[m_wave,line.wll]
;    m_pix=[m_pix,line.posm]
;    m_flag=[m_flag,line.flag]
;    m_ord=[m_ord,replicate(order,nlines)]
;  endfor
;  ord=ord[1:*]
;  m_wave=m_wave[1:*]
;  m_pix=m_pix[1:*]
;  m_flag=m_flag[1:*]
;  m_ord=m_ord[1:*]
;  nm=n_elements(m_wave)
;
;  obase=min(m_ord)
;  nord=max(m_ord)-obase+1
;  maxordx=5
;  maxordy=5
;
;  SDtoFWHM=2.d0*sqrt(2.d0*alog(2.d0))
;  !p.region=0 & !p.position=0
;  !p.multi=0
;  device,decomposed=0 & colors
;  iter=0
;  xmin=-0.5>min(resid)
;  xmax= 0.5<max(resid)
;next:
;
;disp_2d,m_wave,m_pix,m_ord,m_flag,maxordx,maxordy,npix, $
;        nord,obase,wvc,residual=resid
;
;  dx=(xmax-xmin)*0.01d0
;  i=where(m_flag eq 0)
;  h=histogram(resid(i),bin=dx,min=xmin,max=xmax)
;  x=xmin+dindgen(n_elements(h))*dx
;  xx=rebin(x,n_elements(x)*10)
;  yy=gaussfit(x,h,a,xx)
;  plot,x,h,psym=10,xr=[-1,1]*0.5
;  oplot,xx,yy,col=2
;  i=where(abs(resid) gt a(2)*4./SDtoFWHM, ni)
;  if(ni gt 0) then m_flag(i)=1
;  xmin=-6.*a(2)/SDtoFWHM>min(resid)
;  xmax= 6.*a(2)/SDtoFWHM<max(resid)
;  iter=iter+1
;  if(iter lt 4) then goto,next
;
;  !p.multi=[0,1,2]
;  plot,m_ord,resid,psym=1,xs=1,yr=[-1,1]*6*a(2),ys=1, $
;       xtit='Orders',ytit='Residuals'
;  i=where(m_flag ne 0, n_flag)
;  if(n_flag gt 0) then oplot,m_ord(i),resid(i),psym=1,col=3
;  oplot,!x.crange, [a(2),a(2)]*3/SDtoFWHM,col=4
;  oplot,!x.crange,-[a(2),a(2)]*3/SDtoFWHM,col=4
;
;  plot,m_pix,resid,psym=1,xs=1,yr=[-1,1]*6*a(2),ys=1,xtit='Pixels',ytit='Residuals'
;  if(n_flag gt 0) then oplot,m_pix(i),resid(i),psym=1,col=3
;  oplot,!x.crange, [a(2),a(2)]*3/SDtoFWHM,col=4
;  oplot,!x.crange,-[a(2),a(2)]*3/SDtoFWHM,col=4
;
;  for i=0,n_elements(f)-1 do begin
;    restore,f[i]
;    ii=where(m_ord eq ord[i],nii)
;    if(ni eq n_elements(line)) then line.flag=m_flag[ii]
;    save,b,line,file=f[i]
;  endfor
;
;  !p.multi=0
;  wshow,0
;
;end

function kpno_thar,wrange=wrange,path=path,thar=thar,maxline=maxline
  common kpno_atlas,file

  if(keyword_set(path)) then file=path $
  else begin                       file='c:\IDL\kpno_thar'
    help,calls=a
    a=strmid(a(0),strpos(a(0),'<')+1,strlen(a(0))-strpos(a(0),'<')-1)
    if(strpos(a,'/nfs') eq 0) then a=strmid(a,4)
    if(strpos(a,'/export') eq 0) then a=strmid(a,7)
    if(!VERSION.OS eq 'Win32') then delimiter='\' else delimiter='/'
    prefix=strmid(a,0,strpos(a,delimiter,/REVERSE_SEARCH)); kpno_thar directory
    if(strpos(a,delimiter) lt 0) then cd,CURRENT=prefix
    file=prefix+delimiter+'kpno_thar'
  endelse

  restore,file+'.sav'

  thar=readfits(file+'.fits',h,/SILENT) ; Extract ThAr spectrum
  w0=sxpar(h,'CRVAL1')
  dw=sxpar(h,'CDELT1')
  wl=dindgen(n_elements(thar))*dw+w0

  if(keyword_set(wrange)) then begin ;Select wavelength range if specified
    i=where(wl ge wrange[0] and wl le wrange[1], npix)
    if(npix gt 0) then begin
      wl=wl[i] & thar=thar[i]
    endif else return,-1

    i=where(w_th ge wrange[0] and w_th le wrange[1], npix)
    if(npix gt 0) then w_th=w_th[i] else return,-2
  endif

  n=n_elements(thar)
  i=where(thar[1:n-2]-thar[0:n-3] gt 0 and $ ; Select peaks using
          thar[1:n-2]-thar[2:n-1] gt 0 and $ ; 2nd derivative
      thar[1:n-2] gt median(thar), npix)
  if(npix gt 0) then i=i+1 else return,-3

  ii=-1      ;Identify sp. lines
  wth=0.d0
  rth=0.d0
  for j=0,npix-1 do begin
    d=min(abs(wl[i[j]]-w_th),k)
    if(d lt dw) then begin
      ii=[ii,j] & wth=[wth,w_th[k]] & rth=[rth,thar[i[j]]]
    endif
  endfor
  ii=ii[1:*]
  wth=wth[1:*]
  rth=rth[1:*]
  thar=[[wl],[thar]]

  if(keyword_set(maxline)) then begin
    n=n_elements(wth)
    if(n gt (maxline>1)) then begin
      i=(reverse(sort(rth)))[0:maxline-1]
      wth=wth[i]
      rth=rth[i]
    endif
  endif

  return,[[wth],[rth]]
end

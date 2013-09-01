pro get_shear,thar_filename,save_file=save_file,PLOT=plot,inst_mode=inst_mode,exten=exten

  if(n_params() lt 1) then begin
    print,'Syntax: get_shear,thar_filename,[save_file=save_file[' $
         +',inst_mode=inst_mode[,/PLOT]]]'
    print,'   thar_filename - name of the ThAr FITS file'
    print,'   save-file     - REDUCE save file containing order locations'
    print,'   inst_mode     - REDUCE instrument setup'
    print,'   /PLOT         - a flag plotting intermediate fits'
    return
  endif

; Read ThAr file
  if(keyword_set(inst_mode)) then instr=inst_mode $
  else                            instr,'bHROS'
  fits_read,thar_filename,ThAr,hdr,exten=exten
  hdr=modeinfo(hdr,instr)
  ThAr=clipnflip(ThAr,hdr)
  ncol=n_elements(ThAr[*,0])
  nrow=n_elements(ThAr[0,*])

; Restore order location
  if(keyword_set(save_file)) then restore,save_file $
  else                            restore,'QSO_.ord.sav'
  nord=n_elements(orders[0,*])

; Make symmetric extraction window and reduce it a bit
  xwd=min(def_xwd)-2

; Slit shearing will be fit with parabola as function of pixel number
  shear_x=dblarr(3,nord)

  if(keyword_set(plot)) then colors
  for iord=0,nord-1 do begin ; Go through sp. orders

; Local strongest sp. lines
    sp=dblarr(ncol)
    ycen=poly(dindgen(ncol),orders[*,iord])
    for i=col_range[0,iord],col_range[1,iord] do $
      sp[i]=total(ThAr[i,ycen[i]-xwd:ycen[i]+xwd])
    if(keyword_set(plot)) then begin
      plot,sp,xs=1,ys=3
      oplot,!x.crange,mean(sp)+100+[0,0]
    endif
    locmax=where(sp[1:ncol-2] gt sp[0:ncol-3] and $
                 sp[1:ncol-2] gt sp[2:ncol-1] and $
                 sp[1:ncol-2] gt mean(sp)+100, nmax)+1

; Mark lines which are too close to each other to confuse Gauss fit
    ibad=where(locmax[1:nmax-1]-locmax[0:nmax-2] lt 20, nbad)
    if(nbad gt 0) then begin

; For bad pairs keep the strongest of the two lines
      for i=0,nbad-1 do $
        ibad[i]=(sp[locmax[ibad[i]+1]] gt sp[locmax[ibad[i]]])? $
                ibad[i]:ibad[i]+1
      locmax[ibad]=-1
      locmax=locmax[where(locmax ge col_range[0,iord]+13 and $
                          locmax lt col_range[1,iord]-14, nmax)]
    endif

; Keep not more than 31 strongest lines per order
    i=reverse(sort(sp[locmax]))
    i=i[0:40<(nmax-1)]
    nmax=n_elements(i)
    locmax=locmax[i]
    i=sort(locmax)
    locmax=locmax[i]
    if(keyword_set(plot)) then oplot,locmax,sp[locmax],psym=2


    xx=dindgen(19)-9.d0   ; Gaussian fit will be based on 19 pixels
    xcen=dblarr(2*xwd+1)  ; Line centers will be stored along the slit
    shear=dblarr(nmax)    ; Shearing will be determined for each line separately
    for iline=0,nmax-1 do begin ; Loop through the lines
      x=locmax[iline]+xx
      for irow=-xwd,xwd do begin
        s=ThAr[x,irow+ycen[locmax[iline]]] ; Extract short horizontal strip
        yy=disp_gaussfit(x,s,a,x)          ; Do Gauss fit
        xcen[irow+xwd]=a[1]                ; Store line center
      endfor
      if(keyword_set(plot)) then plot,xcen,ys=3,psym=2,tit='Order ' $
                                      +strtrim(iord+1,2)+' line ' $
                                      +strtrim(iline+1,2)+' out of ' $
                                      +strtrim(nmax,2) $
                                      ,xtit='Slit position' $
                                      ,ytit='ThAr line position'
      a=poly_fit(dindgen(2*xwd+1)-xwd,xcen,1,/double) ; Linear fit to slit image
      j=sort(abs(poly(dindgen(2*xwd+1)-xwd,a)-xcen))  ; Outlyers
      a=poly_fit(double(j[0:2*xwd*0.8])-xwd $         ; Fit again using best 80%
                  ,xcen[j[0:2*xwd*0.8]],1,/DOUBLE)
      j=sort(abs(poly(dindgen(2*xwd+1)-xwd,a)-xcen))  ; Outlyers again
      if(keyword_set(plot)) then oplot,j[0:2*xwd*0.8] $
                            ,xcen[j[0:2*xwd*0.8]],psym=2,col=3
      a=poly_fit(double(j[0:2*xwd*0.8])-xwd $         ; Final fit using best 80%
                  ,xcen[j[0:2*xwd*0.8]],1,/DOUBLE)
      if(keyword_set(plot)) then oplot,poly(dindgen(2*xwd+1)-xwd,a)
      if(keyword_set(plot)) then empty
      shear[iline]=a[1]                               ; Store line shear
    endfor
    if(keyword_set(plot)) then plot,locmax,shear,ys=3,xs=1,psym=2 $
                                   ,tit='Slit tilt for order ' $
                                   +strtrim(iord,2),xtit='Pixel #' $
                                   ,ytit='Slit tilt'
    a=poly_fit(locmax,shear,1,/double)
    j=sort(abs(poly(locmax,a)-shear))
    a=poly_fit(locmax[j[0:nmax*0.7]],shear[j[0:nmax*0.7]],2,/DOUBLE)
    j=sort(abs(poly(locmax,a)-shear))
    if(keyword_set(plot)) then oplot,locmax[j[0:nmax*0.7]] $
                                   ,shear[j[0:nmax*0.7]],psym=2,col=3
    a=poly_fit(locmax[j[0:nmax*0.8]],shear[j[0:nmax*0.8]],2,/DOUBLE)
    if(keyword_set(plot)) then oplot,locmax,poly(locmax,a)
    if(keyword_set(plot)) then empty
    shear_x[*,iord]=a
;    ss=get_kbrd(1)
  endfor

; Fix order 5
  i=[0,1,2,3,4,6]
  shear_x[0,5]=poly(5.d0,poly_fit(double(i),shear_x[0,i],1,/DOUBLE))
  shear_x[1,5]=poly(5.d0,poly_fit(double(i),shear_x[1,i],1,/DOUBLE))
  shear_x[2,5]=poly(5.d0,poly_fit(double(i),shear_x[2,i],1,/DOUBLE))
  
  delta_x=intarr(nord)
  for i=0,nord-1 do delta_x[i]=ceil(max(abs(poly([0,ncol] $
                                   ,shear_x[*,i]))*(xwd+1.)))
  f=strmid(thar_filename,0,strpos(thar_filename,'fits')+4) $
           +'.slit_tilt.sav'
  save,delta_x,shear_x,file=f
  print,'Slit tilt functions and the PSF windows are stored in '+f
  
  return
end

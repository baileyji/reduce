Function bandsol,aa,rr,DOUBLE=dbl

  n=n_elements(rr)
  sz=size(aa)
  if(n eq 0 or sz(0) ne 2 or sz(1) ne n) then begin
    print,'bandsol solve a sparse system of linear equations with band-diagonal matrix.'
    print,'Band is assumed to be symmetrix relative to the main diaginal. Usage:'
    print,'res=bandsol(a,r[,/DOUBLE])'
    print,'where a is 2D array [n,m] where n - is the number of equations and m'
    print,'        is the width of the band (3 for tri-diagonal system),'
    print,'        m is always an odd number. The main diagonal should be in a(*,m/2)'
    print,'        The first lower subdiagonal should be in a(1:n-1,m-2-1), the first'
    print,'        upper subdiagonal is in a(0:n-2,m/2+1) etc. For example:'
    print,'               / 0 0 X X X \'
    print,'               | 0 X X X X |'
    print,'               | X X X X X |'
    print,'               | X X X X X |'
    print,'           A = | X X X X X |'
    print,'               | X X X X X |'
    print,'               | X X X X X |'
    print,'               | X X X X 0 |'
    print,'               \ X X X 0 0 /'
    print,'      r is the array of RHS of size n.'
    print,'      /DOUBLE forces the calculations to double precision.'
    return,0
  endif
  nd=sz(2)

  if(keyword_set(dbl)) then begin
    a=double(aa)
    r=double(rr)
  endif else begin
    a=float(aa)
    r=double(rr)
  endelse

  for i=0L,n-2L do begin
    r(i)=r(i)/a(i,nd/2) & a(i,*)=a(i,*)/a(i,nd/2)
    for j=1L,(nd/2)<(n-i-1L) do begin
      r(i+j)=r(i+j)-r(i)*a(i+j,nd/2-j)
      a(i+j,0:nd-j-1)=a(i+j,0:nd-j-1)-a(i,j:nd-1)*a(i+j,nd/2-j)
    endfor
  endfor

  r(n-1L)=r(n-1L)/a(n-1L,nd/2)
  for i=n-1L,1L,-1L do begin
    for j=1L,(nd/2)<i do begin
      r(i-j)=r(i-j)-r(i)*a(i-j,nd/2+j)
    endfor
    r(i-1L)=r(i-1L)/a(i-1L,nd/2)
  endfor
  r(0L)=r(0L)/a(0L,nd/2)

  return,r
end

Pro slit_func_2d,im,ycen,sp,sf,delta_x,shear_x,OVERSAMPLE=oversample $
             ,LAMBDA_SF=lamb_sf,LAMBDA_SP=lamb_sp,IM_OUT=im_out $
             ,USE_COL=use_col,MASK=mask,NOISE=noise,BAD=jbad $
             ,MODEL_ONLY=model,WING_SMOOTH_FACTOR=wing_smooth_factor $
             ,UNCERTAINTY=unc
common bandsolv,band_solv_name

  if(not keyword_set(band_solv_name)) then begin
    help,calls=a
    delimiter=path_sep()
;    i1=strpos(a[0],delimiter)
    i1=strpos(a[0],'<')+1
    prefix=file_dirname(strmid(a[0],i1))
;    a=strmid(a(0),strpos(a(0),'<')+1,strlen(a(0))-strpos(a(0),'<')-1)
;    if(strpos(a,'/nfs') eq 0) then a=strmid(a,4)
;    if(strpos(a,'/export') eq 0) then a=strmid(a,7)
;    if(!VERSION.OS eq 'Win32') then delimiter='\' else delimiter='/'
;    prefix=strmid(a,0,strpos(a,delimiter,/REVERSE_SEARCH))   ;slit_func directory
;    if(strpos(a,delimiter) lt 0) then cd,CURRENT=prefix
    band_solv_name=prefix+delimiter+'bandsol.so.'   $
                                   +!version.os+'.' $
                                   +!version.arch+'.' $
                                   +strtrim(!version.MEMORY_BITS,2)
  endif
  if(not keyword_set(oversample)) then oversample=1
  if(oversample lt 1) then oversample=1
  if(not keyword_set(mask)) then begin
    mmsk=byte(im*0)+1B
  endif else begin
    if((size(mask))(0) ne (size(im))(0) or $
       (size(mask))(1) ne (size(im))(1) or $
       (size(mask))(2) ne (size(im))(2)) then begin
      print,'SLIT_FUNC: Mask must have the same size as the image'
      stop
    endif
    mmsk=mask
  endelse
  osample=long(oversample)
  oind=lindgen(osample+1L)*(osample+2L)
  weight=1./double(osample)
  if(not keyword_set(lamb_sf)) then lamb_sf=0.1

  for reject=1,1 do begin
    if(keyword_set(use_col)) then begin
      imm=im(use_col,*)
      yycen=ycen(use_col)
      msk=mmsk(use_col,*)
    endif else begin
      use_col=indgen(n_elements(im(*,0)))
      imm=im
      yycen=ycen
      msk=mmsk
    endelse

    sz=size(imm)
    ncol=sz(1)
    nrow=sz(2)
    n_sf=(nrow+1L)*osample+1L
    nx=2L*delta_x+1L

    norm=n_elements(msk)/total(long(msk))

    yslit=dindgen(n_sf)/float(osample)-1.

    Akl_ind_sf=(lindgen(osample+1L)*n_sf)#replicate(1L,osample+1L)
    for i=0L,osample do Akl_ind_sf[*,i]=Akl_ind_sf[*,i]+(osample-i)*n_sf+i
    Akl_ind_sf=reform(Akl_ind_sf,n_elements(Akl_ind_sf))

;    Akl_ind_sp=(lindgen(delta_x+1L)*ncol)#replicate(1L,delta_x+1L)
;    for i=0L,delta_x do Akl_ind_sp[*,i]=Akl_ind_sp[*,i]+(delta_x-i)*ncol+i
    Akl_ind_sp=(lindgen(nx+1L)*ncol)#replicate(1L,nx+1L)
    for i=0L,nx do Akl_ind_sp[*,i]=Akl_ind_sp[*,i]+(nx-i)*ncol+i

                                                         ; We assume that sp and sf are supplied
    if(not keyword_set(model)) then begin                ; if only a model is requested
      sf=total(imm*msk,1)
      if(osample gt 2 and n_elements(sf) gt 5) then sf=median(sf,5) ; the spectrum
      if(mean(total(imm,2)) lt 1.d3) then $              ; robust guess for sf
      sf=exp(-((dindgen(nrow)-nrow/2.)/(nrow/4.))^2)     ; in case of low S/N
      sf=sf/total(sf)                                    ; slit function
      sp=total((imm*msk)*(replicate(1.,ncol)#sf),2)*norm ; Initial guess for
      if(osample gt 2) then sp=median(sp,5)              ; the spectrum
      sp=sp/total(sp)*total(imm*msk)
    endif

;    outliers=0L
    if(keyword_set(noise)) then dev=noise $
    else                        dev=sqrt(total(msk*(imm-sp#sf)^2)/total(msk))
    j=where(abs(imm-sp#sf) gt 3.*dev,nj)
    if(nj gt 0) then begin
      msk(j)=0B
;      outliers=outliers+nj
;      stop
    endif

;
; Computing the omega
;
; S_xy = Sum_(ix=-delta_x,delta_x) P_(x+ix) Sum_iy omega(ix,iy,x,y) * L(iy)
;
; Define omega and usefull indexing arrays used for optimization below
; {x,y} is the CCD pixel where the signal is computed, iy is subpixel index
; and ix is the contribution from the PSF centered on pixel x+ix
;
    omega=dblarr(nx,osample+1,ncol,nrow)      ; Intermediate values of omega are identical
    omega[delta_x,*,*,*]=weight
;
; We start by computing the omega. Note, that omega only depends on the geometry and therefore,
; does not change throughout the iterations:
;
    ix=delta_x
    sf_ind=0L                                  ; Index of slit function falling to row Y in cols x-delta_x:x+delta_x
    sf_ind_ind=lonarr(ncol)                    ; Starting point of indeces corresponding to pixel x
    iybottom=intarr(ncol)
    for x=0L,ncol-1L do begin                  ; Omega depends on x through the shift of the central line
      ix1=-(x<delta_x)
      ix2=(ncol-x-1L)<delta_x
      sf_ind_ind[x]=n_elements(sf_ind)
      for ix=ix1,ix2 do begin
        jx=ix+delta_x
        weight_bottom=yycen[x+ix] mod weight   ; Bottom of a pixel weight is the amount by which a
                                               ; subpixel sticks out into the pixel above
        weight_top=weight-weight_bottom        ; Top of a pixel weight is the amount of a subpixel
                                               ; left in the current pixel
        omega[jx,0,x,*]=weight_bottom          ; Correct bottom weights
        omega[jx,osample,x,*]=weight_top       ; Correct top weights
        i1=fix((1.d0-yycen[x+ix])*osample)     ; The very bottom part of the PSF sticks below the last
        if(ix eq 0) then iybottom[x]=i1        ; The coordinate of the bottom subpixel is needed for setting up SLEs
        if(i1 gt 0) then omega[jx,0:i1-1,x,0]=0; CCD pixel. If yycen==0 this part is exactly osample
        omega[jx,0:i1,x,0]=weight_bottom       ; subpixels. Larger yycen shifts the PSF up leaving smaller tail
        sf_ind=[sf_ind,i1+lindgen(osample+1L)]
        i1=fix(yycen[x+ix]*osample)            ; The very top pixel behaves in the opposite way. Larger yycen
        if(i1 lt osample) then omega[jx,i1+1:*,x,0]=0 ; leaves larger tail sticking outside the data.
        omega[jx,i1:*,x,0]=weight_top
        if(ix ne 0) then omega[jx,*,x,*]=0.d0
      endfor
    endfor

    sf_ind=sf_ind[1:*]
    sf_ind_ind=[sf_ind_ind-1L,n_elements(sf_ind)]
;
; For a tilted slit we assume that the angle is known as function of x
; so we can modify the omega accordingly
;
    for x=0L,ncol-1L do begin
      ix1=-(x<delta_x)
      ix2=(ncol-x-1L)<delta_x
      yy=dblarr((osample+1L)*nrow)
      o=reform(omega[delta_x,*,x,*],(osample+1L)*nrow)
      for iy=0L,(osample+1L)*nrow-1L do yy[iy]=total(o[0:iy])-o[iy]*0.5d0-nrow*0.5d0-yycen[x]
      shear=shear_x[x]*yy
      ix1=fix(shear)
      ii=where(shear ge 0, nii, complement=jj, ncomplement=njj)
      ix2=ix1
      if(nii gt 0) then ix2[ii]=ix2[ii]+1
      if(njj gt 0) then ix2[jj]=ix2[jj]-1
      omega[delta_x,*,x,*]=0.d0
      for iy=0L,(osample+1L)*nrow-1L do begin
        y=iy/(osample+1L)
        iyy=iy mod (osample+1L)
        if(x+ix1[iy] ge 0L and x+ix1[iy] lt ncol) then $
          omega[ix1[iy]+delta_x,iyy,x+ix1[iy],y]=o[iy]*(1.d0-abs(shear[iy]-fix(shear[iy])))
        if(x+ix2[iy] ge 0L and x+ix2[iy] lt ncol) then $
          omega[ix2[iy]+delta_x,iyy,x+ix2[iy],y]=o[iy]*abs(shear[iy]-fix(shear[iy]))
      endfor
    endfor
    u=replicate(1,osample+1)
    if(keyword_set(model)) then  goto,model_only
;
; Iterations start here
;
    iter=0
next:
    iter=iter+1                                ; Iteration counter

;===============================================================================
; (A) Detailed version
;     Compute the new PSF
;
;time0=systime(1)
    Akl=dblarr(n_sf,2*osample+1)               ; Initialize compact matrix
;    Aij=dblarr(n_sf,n_sf)
    Bl=dblarr(n_sf)                            ; and RHS
    for x=0L,ncol-1L do begin
      ix1=-(x<delta_x)
      ix2=(ncol-x-1L)<delta_x
      jx1=ix1+delta_x
      jx2=ix2+delta_x
      iy1=iybottom[x]-osample
      iy2=iy1+osample
      for y=0L,nrow-1L do begin
        iy1=iy1+osample
        iy2=iy1+osample
        ind=Akl_ind_sf+iy1
        dummy=total((sp[x+ix1:x+ix2]#u) $
                    *reform(omega[jx1:jx2,*,x,y]),1)
;        Aij[iy1:iy2,iy1:iy2]=Aij[iy1:iy2,iy1:iy2]+(dummy#dummy)*msk[x,y]
        Akl[ind]=Akl[ind]+(dummy#dummy)*msk[x,y]
        Bl[iy1:iy2]=Bl[iy1:iy2]+dummy*msk[x,y]*imm[x,y]
      endfor
    endfor

;time1=systime(1)

    lambda = lamb_sf*total(Akl(*,osample))/n_sf
    if(keyword_set(wing_smooth_factor)) then $
      lambda=lambda*(1.+wing_smooth_factor*(2.d0*dindgen(n_sf)/(n_sf-1)-1.d0)^2) $
    else $
      lambda=replicate(lambda,n_sf)

; 1st order Tikhonov regularization (minimum 1st derivatives)
; Add the following 3-diagonal matrix * lambda:
;  1 -1  0  0  0  0
; -1  2 -1  0  0  0
;  0 -1  2 -1  0  0
;  0  0 -1  2 -1  0
;      .  .  .

;    main_diag=lindgen(n_sf)*(n_sf+1L)
;    uppX_diag=lindgen(n_sf-1L)*(n_sf+1L)+1L
;    lowX_diag=lindgen(n_sf-1L)*(n_sf+1L)+n_sf
;    Aij[0,0]=Aij[0,0]+lambda
;    Aij[n_sf-1,n_sf-1]=Aij[n_sf-1,n_sf-1]+lambda
;    Aij[main_diag[1:n_sf-2]]=Aij[main_diag[1:n_sf-2]]+2.d0*lambda
;    Aij[uppX_diag]=Aij[uppX_diag]-lambda
;    Aij[lowX_diag]=Aij[lowX_diag]-lambda

;    for i=0L,osample do Akl[osample-i:*,i]=Aij[lindgen(n_sf-osample+i)*(n_sf+1L)+n_sf*(osample-i)]
;    for i=1L,osample do Akl[0:n_sf-i-1,i+osample]=Aij[lindgen(n_sf-i)*(n_sf+1L)+i]

;    la_ludc,Aij,index,/DOUBLE
;    sf1=la_lusol(Aij,index,Bl,/DOUBLE)
;    sf1=sf1/total(sf1)*osample


;    i=CALL_EXTERNAL(band_solv_name, 'bandsol', $
;                    Akl, Bl, n_sf, 2L*osample+1L)
;    sf=Bl/total(Bl)*osample
;stop

;    Akl(        0L,osample)=Akl(           0L,osample)+lambda    ; +lambda to the upper-left element
;    Akl(   n_sf-1L,osample)=Akl(      n_sf-1L,osample)+lambda    ; and to the lower-right
;    Akl(1L:n_sf-2L,osample)=Akl(   1L:n_sf-2L,osample)+2.*lambda ; +2*lambda to the rest of the main diagonal
;    Akl(0L:n_sf-2L,osample+1L)=Akl(0L:n_sf-2L,osample+1L)-lambda ; -lambda to the upper sub-diagonal
;    Akl(1L:n_sf-1L,osample-1L)=Akl(1L:n_sf-1L,osample-1L)-lambda ; -lambda to the lower sub-diagonal
    Akl(  0,osample)=Akl(  0,osample)+lambda(  0); +lambda to the upper-left element
    Akl(n_sf-1,osample)=Akl(n_sf-1,osample)+lambda(n_sf-1L); and to the lower-right
    Akl(1L:n_sf-2L,osample)=Akl(1L:n_sf-2L,osample)+2.*lambda(1L:n_sf-2L); +2*lambda to the rest of the main diagonal
    Akl(0L:n_sf-2L,osample+1L)=Akl(0L:n_sf-2L,osample+1L)-lambda(0L:n_sf-2L); -lambda to the upper sub-diagonal
    Akl(1L:n_sf-1L,osample-1L)=Akl(1L:n_sf-1L,osample-1L)-lambda(1L:n_sf-1L); -lambda to the lower sub-diagonal

;
; 2nd order Tikhonov regularization (minimum 2nd derivative)
; Add the following 5-diagonal matrix * lambda:
;  1 -2  1  0  0  0
; -2  5 -4  1  0  0
;  1 -4  6 -4  1  0
;  0  1 -4  6 -4  1
;      .  .  .

;    lambda=0.1*lambda
;    Akl(        0L,osample)   =Akl(        0L,osample)+1.*lambda    ; Main diagonal
;    Akl(   n_sf-1L,osample)   =Akl(   n_sf-1L,osample)+1.*lambda
;    Akl(        1L,osample)   =Akl(        1L,osample)+5.*lambda
;    Akl(   n_sf-2L,osample)   =Akl(   n_sf-2L,osample)+5.*lambda
;    Akl(2L:n_sf-3L,osample)   =Akl(2L:n_sf-3L,osample)+6.*lambda
;    Akl(        0L,osample+1L)=Akl(        0L,osample+1L)-2.*lambda ; upper sub-diagonal
;    Akl(   n_sf-2L,osample+1L)=Akl(   n_sf-2L,osample+1L)-2.*lambda
;    Akl(1L:n_sf-3L,osample+1L)=Akl(1L:n_sf-3L,osample+1L)-4.*lambda
;    Akl(        1L,osample-1L)=Akl(        1L,osample-1L)-2.*lambda ; lower sub-diagonal
;    Akl(   n_sf-1L,osample-1L)=Akl(   n_sf-1L,osample-1L)-2.*lambda
;    Akl(2L:n_sf-2L,osample-1L)=Akl(2L:n_sf-2L,osample-1L)-4.*lambda

;    sf=bandsol(Akl,Bl,/DOUBLE)
;    sf=sf/total(sf)*osample
    i=CALL_EXTERNAL(band_solv_name, 'bandsol', $
                    Akl, Bl, n_sf, 2L*osample+1L)
    sf=Bl/total(Bl)*osample

;time2=systime(1)

    sp_old=sp
    r=sp

    Akl=dblarr(ncol,2*nx+1)
;    Aij=dblarr(ncol,ncol)
    Bl=dblarr(ncol)
    for x=0L,ncol-1L do begin
      x1=(x-delta_x)>0L
      x2=(x+delta_x)<(ncol-1L)
      ix1=-(x<delta_x)
      ix2=(ncol-x-1L)<delta_x
      jx1=ix1+delta_x
      jx2=ix2+delta_x
      ind=Akl_ind_sp[0:x2-x1,0:x2-x1]+x1
      for y=0L,nrow-1L do begin
        dummy =total(transpose(reform(omega[jx1:jx2,*,x,y])) $
                    *reform(sf[sf_ind[sf_ind_ind[x]:sf_ind_ind[x+1]-1L]+y*osample] $
                    ,osample+1,ix2-ix1+1),1)
;        Aij[x1:x2,x1:x2]=Aij[x1:x2,x1:x2]+(dummy#dummy)*msk[x,y]
        Akl[ind]=Akl[ind]+(dummy#dummy)*msk[x,y]
        Bl[x1:x2]=Bl[x1:x2]+dummy*msk[x,y]*imm[x,y]
      endfor
    endfor

;Smoothing in the spectral direction is requested
    if(keyword_set(lamb_sp)) then begin
      lambda = lamb_sp*total(sp)/ncol
      Akl(        0L,delta_x)=Akl(           0L,delta_x)+lambda    ; +lambda to the upper-left element
      Akl(   ncol-1L,delta_x)=Akl(      ncol-1L,delta_x)+lambda    ; and to the lower-right
      Akl(1L:ncol-2L,delta_x)=Akl(   1L:ncol-2L,delta_x)+2.*lambda ; +2*lambda to the rest of the main diagonal
      Akl(0L:ncol-2L,delta_x+1L)=Akl(0L:ncol-2L,delta_x+1L)-lambda ; -lambda to the upper sub-diagonal
      Akl(1L:ncol-1L,delta_x-1L)=Akl(1L:ncol-1L,delta_x-1L)-lambda ; -lambda to the lower sub-diagonal
    endif

;akl1=akl*0
;for i=0L,nx do Akl1[nx-i:*,i]=Aij[lindgen(ncol-nx+i)*(ncol+1L)+ncol*(nx-i)]
;for i=1L,nx do Akl1[0:ncol-i-1,i+nx]=Aij[lindgen(ncol-i)*(ncol+1L)+i]
;stop

;    la_ludc,Aij,index,/DOUBLE
;    sp1=la_lusol(Aij,index,Bl,/DOUBLE)


;Get the new spectrum
;    sp=bandsol(Akl,Bl,/DOUBLE)
    i=CALL_EXTERNAL(band_solv_name, 'bandsol', $
                    Akl, Bl, ncol, 2L*nx+1L)
    sp=Bl

;stop

;time3=systime(1)

;Reconstruct the model
    dev_new=0.d0 ; New deviation estimate
    model=imm*0
    for x=0L,ncol-1L do begin
      ix1=-(x<delta_x)
      ix2=(ncol-x-1L)<delta_x
      jx1=ix1+delta_x
      jx2=ix2+delta_x
      for y=0L,nrow-1L do begin
        model[x,y]=total((sp[x+ix1:x+ix2]#u)*reform(omega[jx1:jx2,*,x,y]) $
              *reform(sf[sf_ind[sf_ind_ind[x]:sf_ind_ind[x+1]-1L]+y*osample] $
              ,osample+1,ix2-ix1+1))
      endfor
      dev_new=dev_new+total((imm[i,*]-model[x,*])^2)
    endfor
    if(iter gt 1) then dev=sqrt(dev_new/total(msk))

;Locate and mask outliers
    for x=0,ncol-1L do begin
      if(iter gt 1) then begin
        j=where(abs((imm[x,*]-model[x,*])) gt 3.*dev,nj,COMPLEMENT=b)
        if(nj gt 0) then begin
          msk[x,j]=0B
;          outliers=outliers+nj
;          stop
        endif
;        if(nj lt nrow) then msk(x,b)=1B
        if(nj lt nrow) then msk[x,b]=1B*(mmsk[use_col,*])[x,b]
      endif
    endfor


;time4=systime(1)

;  !p.multi=[0,1,2]
;  x=dindgen(n_elements(im(0,*)))+0.5
;  plot,yslit,sf,xs=3,ys=3,title='Slit Function. Iteration='+strtrim(iter,2)
;  for i=0,n_elements(im(*,0))-1 do oplot,x-ycen(i)+0.5,im(i,*)/sp(i),psym=3
;  plot,sp,xs=3,ys=3,title='Spectrum'
;  wshow,0
;  !p.multi=0
;   print,'Maximum change in the spectrum is:',max(abs(sp-sp_old)/max(sp))
;  print,time1-time0,time2-time1,time3-time2,time4-time3
;  stop

;    if(iter lt 4 and (max(abs(sp-sp_old)/max(sp)) gt 1.d-4 or outliers gt 0L)) then goto,next
    if(iter lt 8 and max(abs(sp-sp_old)/max(sp)) gt 1.d-5) then goto,next

model_only:
    jbad=0L
    if(arg_present(im_out)) then begin
;Reconstruct the model
      unc=fltarr(ncol)
      im_out=im*0
      for x=0L,ncol-1L do begin
        ix1=-(x<delta_x)
        ix2=(ncol-x-1L)<delta_x
        jx1=ix1+delta_x
        jx2=ix2+delta_x
        for y=0L,nrow-1L do begin
          im_out[x,y]=total((sp[x+ix1:x+ix2]#u)*reform(omega[jx1:jx2,*,x,y]) $
                *reform(sf[sf_ind[sf_ind_ind[x]:sf_ind_ind[x+1]-1L]+y*osample] $
                ,osample+1,ix2-ix1+1))
        endfor
      endfor
;      stop

;Locate and mask outliers
      for x=0,ncol-1L do begin
        j=where(abs((im[x,*]-im_out[x,*])) gt 3.*dev,nj,COMPLEMENT=b)
        if(nj gt 0) then begin
          msk[x,j]=0B
          unc[x]=sqrt(total((im[x,j]-im_out[x,j])^2)/nj)
        endif else unc[x]=-1
        if(nj lt nrow) then msk[x,b]=1B*(mmsk[use_col,*])[x,b]
      endfor

    endif
  endfor
  if(n_elements(jbad) gt 1) then jbad=jbad(1:n_elements(jbad)-1) $
  else                           jbad=-1

  return
end

;restore,'slit_func.sav'
;jgood=0
;slit_func,sf,yc,sp,sfsm,OVERSAMPLE=10,IM_OUT=sfbin,LAMBDA_SF=10. $
;            ,USE_COL=jgood,BAD=jbad,MASK=msk
;end

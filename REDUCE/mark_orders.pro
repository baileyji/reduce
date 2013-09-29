Pro locate_clusters,im,n,x,y,nx,ny,filter,MASK=mask,NOISE=noise

  offset=keyword_set(noise)?noise:0.0
  imm=bytarr(nx,ny)
  for ix=0,nx-1 do begin
    mm=float(im[ix,*])-offset-opt_filter(float(im[ix,*]),filter)
;    h=median(mm>0)*0.5
    h=median(mm[where(mm gt 0)])*0.5
    i=where(mm gt h+1, ni)
    if(ni gt 0) then imm[ix,i]=1B
;    if(ix eq nx/2) then stop
  endfor
  if(keyword_set(mask)) then imm=imm*mask
  y=where(imm gt 0B, ni)
  if(ni le 0) then begin
    print,'MARK_ORDERS: This should never happen'
    stop
  endif

;  imm=float(im)
;  for ix=0,nx-1 do imm[ix,*]=float(im[ix,*])-opt_filter(float(im[ix,*]),filter)
;  mm=minmax(imm)
;  nn=0.001*n
;  bin=(mm(1)-mm(0))/nn
;  h=histogram(imm,bin=bin)
;  m=max(h,imax)
;  h=[h[0:imax],h[imax-1-indgen(imax)]]
;  xx=(dindgen(2*imax+1)-imax)*bin
;  z=gaussfit(xx,h,ag,nterms=3)
;  h=ag(2)
;  y=where(imm gt 2.d0*h)

  x=y mod nx
  y=y/nx

  return
end

Pro mark_orders, im,orders, POWER=power, FILTER=filter $
               , ERROR=ord_err, THRES=thres, PLOT=iplot $
               , POLARIM=polarim, MASK=mask, MANUAL=MANUAL $
               , NOISE=noise, COLOR=color, CLUSTER_LIMITS=cl_limits

  if(n_params() lt 2) then begin
    print,'Syntax: mark_orders,im,orders[,POWER=power[,FILTER=filter $'
    print,'                  [,THRES=thres[,/PLOT[,/POLARIM[,MASK=mask $'
    print,'                  [,/MANUAL[,NOISE=noise]]]]]]]]'
    print,'        im      - 2D array [nx,ny] with spectral orders oriented horizontally'
    print,'        orders  - 2D array [power,norders] with polynomial approximation to order shape'
    print,'        power   - power of the polynomial fit. Default is 2 (parabola).'
    print,'        filter  - parameter of the optimal filter used to detect high pixels. Default: 40.'
    print,'        ord_err - [norders] error of the polynomial fit to spectral orders'
    print,'        thres   - minimum group of pixels to be considered is part of spectral order. Default:400'
    print,'        /PLOT   - if set, execute mark_orders in interactive mode'
    print,'        /POLARIM- if set, assumes pairs of orders prodeced by beamsplitter'
    print,'        mask    - 2D byte array [nx,ny] marks good (1B) and bad (0B) pixels'
    print,'        /MANUAL - all cluster merging is done in interactive mode'
    print,'        noise   - readout noise of the CCD may be important if the S/N'
    print,'                  in the image used is too low'
    print,'        cl_limits - optional parameter containing starting/ending columns of each cluster'
    print,'       /COLOR  - use colors when marking clusters.' 
   return
  endif

  sz=size(im)
  nx=sz(1)
  ny=sz(2)
  n=n_elements(im)

;Get X and Y coordinates of all pixels sticking above the filtered image
  if(not keyword_set(filter)) then filter=20.d0
  locate_clusters,im,n,x,y,nx,ny,filter,MASK=mask,NOISE=noise

;Locate and mark clusters of pixels
  if(keyword_set(iplot)) then begin
    if(keyword_set(color)) then colors else loadct,0,/SILENT
  endif
  if(not keyword_set(thres)) then thres=400
  index=reduce_cluster(x,y,nx,ny,nregions=nregions,thres=thres,PLOT=iplot)
  ind=where(index gt 0, n)
  if(n gt 0) then begin
    x=x(ind)
    y=y(ind)
    index=index(ind)
  endif else begin
    print,'MARK_ORDERS: No clusters identified'
    stop
  endelse
  
  if nregions le 1 then begin
    print,'MARK_ORDERS: No clusters identified'
    stop
  endif

;Find all the clusters crossing the bottom boundary
  icross=where(y eq 0, ncross)
  if(ncross gt 0) then begin         ;Found some pixels
    icross=index(icross)             ;Get cluster numbers
    icross=icross(sort(icross))      ;Sort for unique to work
    icross=icross(uniq(icross))      ;Keep different cluster numbers only
    ncross=n_elements(icross)        ;Number of suspicious clusters
    for i=0,ncross-1 do begin        ;Loop through clusters
      if(icross(i) gt 0) then begin
        ind=where(icross(i) eq index)  ;Find all the pixels in the border cluster
        bind=ind(where(y[ind] eq 0,nbind));All the border pixels in this cluster
        for ii=0L,nbind-1 do begin
          j=ind[where(x[ind] eq x[bind[ii]])]
          index[j]=0
        endfor
;        xbmin=min(x[bind])             ;Smallest X
;        xbmax=max(x[bind])             ;Largest X
;        xbmean=0.5*(xbmin+xbmax)       ;Mean X for border pixels
;        xmean=mean(x[ind])             ;Mean X for all pixels in this cluster
;        if(xmean lt xbmean) then begin ;Most of the cluster is located to the left of xbmean
;          j=ind[where(x[ind] ge xbmin)];Find the ones that can be influenced by the edge
;          index[j]=0                   ;Whipe them out
;        endif else begin               ;Most of the cluster is located to the right of xbmean
;          j=ind[where(x[ind] le xbmax)];Find the ones that can be influenced by the edge
;          index[j]=0                   ;Whipe them out
;        endelse
      endif
    endfor
  endif

;Find all the clusters crossing the top boundary
  icross=where(y eq ny-1, ncross)
  if(ncross gt 0) then begin         ;Found some pixels
    icross=index[icross]             ;Get cluster numbers
    icross=icross[sort(icross)]      ;Sort for unique to work
    icross=icross[uniq(icross)]      ;Keep different cluster numbers only
    ncross=n_elements(icross)        ;Number of suspicious clusters
    for i=0,ncross-1 do begin        ;Loop through clusters
      if(icross[i] gt 0) then begin
        ind=where(icross[i] eq index)  ;Find all the pixels in the border cluster
        bind=ind(where(y[ind] eq ny-1,nbind));All the border pixels in this cluster
        for ii=0L,nbind-1 do begin
          j=ind[where(x[ind] eq x[bind[ii]])]
          index[j]=0
        endfor
;        xbmin=min(x[bind])             ;Smallest X
;        xbmax=max(x[bind])             ;Largest X
;        xbmean=0.5*(xbmin+xbmax)       ;Mean X for border pixels
;        xmean=mean(x[ind])             ;Mean X for all pixels in this cluster
;        if(xmean lt xbmean) then begin ;Most of the cluster is located to the left of xbmean
;          j=ind[where(x[ind] ge xbmin)];Find the ones that can be influenced by the edge
;          index[j]=0                   ;Whipe them out
;        endif else begin               ;Most of the cluster is located to the right of xbmean
;          j=ind[where(x[ind] le xbmax)];Find the ones that can be influenced by the edge
;          index[j]=0                   ;Whipe them out
;        endelse
      endif
    endfor
  endif

  for ireg=1,nregions-1 do begin
    j=where(index eq ireg, nj)
    if(nj lt thres and nj gt 0) then index[j]=0
  endfor

  xx=dindgen(nx)
  if(not keyword_set(power)) then power=2

  skip=[0,0]
  silent=0

  loadct,0,/SILENT
  if(keyword_set(iplot)) then begin
    window,1
    wset,0
  endif
next:
  if(keyword_set(iplot)) then if(silent ne 1) then display,im,/log
;  if(keyword_set(iplot)) then begin
;    display,im,/log
;    j=where(index gt 0)
;    oplot,x[j],y[j],psym=3,col=0
;  endif
  yy=fltarr(nx,nregions-1)
  ord_err=dblarr(nregions-1)
  orders=dblarr(power+1,nregions-1)
  xoffset=dblarr(nregions-1)
  xscale=dblarr(nregions-1)
  ind=0
  cl_limits=intarr(2,nregions-1)   ; starting and ending columns of clusters
  for ireg=1,nregions-1 do begin
    j=where(index eq ireg, nj)
    if(nj gt 0) then begin ; Check if cluster is not empty
      pwr=power<(n_elements(uniq(x[j],sort(x[j])))-1); restrict the power if not enough x's
      if(pwr gt 0) then begin                        ; Check if cluster includes
        index[j]=ind+1                               ; enough different pixels
        xoffset[ind]=double(mean(x[j]))              ; Offset argument in poly_fit
        xscale[ind]=1.d0/max(abs(x[j]-xoffset[ind]))
        coef=poly_fit((x[j]-xoffset[ind])*xscale[ind],y[j],pwr,/double)
                                                     ; Check if this cluster is
                                                     ; not a vertical artifact
        if(abs((poly_fit((x[j]-xoffset[ind]),y[j],1,/double))(1)) lt $
                               4.*ny/nx) then begin
          cl_limits[*,ind]=minmax(x[j])
          orders[0:pwr,ind]=reform(coef)
          yy[*,ind]=poly((xx-xoffset[ind])*xscale[ind],orders[*,ind])
          ord_err[ind]=sqrt(total((y[j]-(poly((x[j]-xoffset[ind])*xscale[ind] $
                                        ,orders[*,ind])))^2)/nj)
          if(keyword_set(iplot) and silent ne 1) then oplot,xx,yy[*,ind],col=0
          ind=ind+1
        endif
      endif else index[j]=0
    endif
  endfor
  nregions=ind+1

  if(nregions le 2) then begin
    mean_sep=1
  endif else begin   ; take the mean distance between orders but be restrictive
    yo=dblarr(nregions-1)
    for kkk=0,nregions-2 do yo[kkk]=poly((double(nx)-xoffset[kkk])*xscale[kkk] $
                                        , orders[*,kkk])
    kkk=where(orders[0,0:nregions-2] ge 0 and yo lt ny, nkkk)
    if(nkkk gt 2) then begin
      yo=minmax(middle(orders(0,sort(orders[0,kkk])),1,/poly))
      mean_sep=(yo[1]-yo[0])/nkkk
    endif else mean_sep=10
;    mean_sep=(total(orders[0,1:nregions-2]- $
;                    orders[0,0:nregions-3])/(nregions-2))<40
  endelse
;  print,mean_sep

; Cluster merging loop.
; We loop through all clusters and identify the overap
; The results are stored in array imerge where the imerge(0,*) contains the
; number of overlaping pixels while imerge(1,*) points to the overlapping
; cluster number. The limits for the overlap region are stored in imerge(2:3,*).
  imerge=intarr(4,nregions-1)
  for i=1,nregions-2 do begin
    j=intarr(nregions-1-i)
    lm=intarr(2,nregions-1-i)
    for ii=i+1,nregions-1 do begin
      iii=where(yy[*,i -1] ge 0 and yy[*,i -1] lt ny and $
                yy[*,ii-1] ge 0 and yy[*,ii-1] lt ny and $
                abs(yy[*,i-1]-yy[*,ii-1]) lt 0.25*mean_sep, niii)
      j[ii-i-1]=niii
      lm[*,ii-i-1]=minmax(iii)
    endfor
    imerge[0,i-1]=max(j,niii)
    imerge[1,i-1]=niii+i
    imerge[2:3,i-1]=lm[*,niii]
    ii=where(skip[0,*] eq i and skip[1,*] eq niii+i, miii)
    if(miii gt 0) then imerge[0,i-1]=0
  endfor

skip_pair:
  lmerge=max(imerge[0,*],nmerge)  ; Largest overlap of lmerge pixels is between
                                  ; clusters nmerge and imerge(1,nmerge)
  if(lmerge gt 30) then begin
    jmerge=imerge[1,nmerge]       ; Cluster nmerge has largest overlap with jmerge
    if(keyword_set(iplot)) then begin
;=============================================
      x_plot_old=!x & y_plot_old=!y & pold=!p
      wset,1
      j=where(index eq jmerge+1)
      cluster_j_x=x[j]
      cluster_j_y=y[j]
      j=where(index eq nmerge+1)
      cluster_n_x=x[j]
      cluster_n_y=y[j]
      xr=minmax([cluster_j_x,cluster_n_x])
      xr=[(xr[0]-3)>0,(xr[1]+3)<(nx-1)]
      yr=minmax([cluster_j_y,cluster_n_y])
      yr=[(yr[0]-3)>0,(yr[1]+3)<(ny-1)]
      display,im,/log,xr=xr,yr=yr,tit='Clusters '+strtrim(jmerge,2)+ $
                                              '+'+strtrim(nmerge,2)
      oplot,cluster_j_x,cluster_j_y,psym=1,col=0
      oplot,cluster_n_x,cluster_n_y,psym=4,col=!d.n_colors-1
      oplot,xx,yy[*,nmerge],col=128
      oplot,xx,yy[*,jmerge],col=128,line=2
;      s=get_kbrd(1)
      x_plot_new=!x & y_plot_new=!y & pnew=!p
      wset,0
      !x=x_plot_old & !y=y_plot_old & !p=pold
;=============================================
bad_answer:
;      oplot,xx,yy[*,nmerge] & oplot,xx,yy[*,jmerge],line=2
      if((imerge[3,nmerge]<cl_limits[1,nmerge]) $    ; If 90% overlap - merge
        -(imerge[2,nmerge]>cl_limits[0,nmerge]) gt $ ; without asking
          0.9*(cl_limits[1,nmerge]-cl_limits[0,nmerge]) or $
         (imerge[3,nmerge]<cl_limits[1,jmerge]) $
        -(imerge[2,nmerge]>cl_limits[0,jmerge]) gt $
          0.9*(cl_limits[1,jmerge]-cl_limits[0,jmerge]) and $
          not keyword_set(MANUAL)) then begin
          s='m'
          silent=1
      endif else begin
        silent=0
        oplot,xx,yy[*,nmerge] & oplot,xx,yy[*,jmerge],line=2
        print,'Clusters '+strtrim(nmerge,2)+'+'+strtrim(jmerge,2)+ $
              ' seem to belong to the same group. You can also select'+ $
              ' a different pair with (C)ursor.'
        print,'Shall they be (M)erged, (I)gnored or shall I discard (S)olid '+ $
              'line or (D)ashed line? '
        s=get_kbrd(1)
      endelse
      if(s eq 'm' or s eq 'M') then begin
        print,'Merging clusters '+strtrim(nmerge,2)+' and '+strtrim(jmerge,2) $
             +' of total '+strtrim(nregions,2)
        j=where(index eq nmerge+1)
        index[j]=jmerge+1
      endif else if(s eq 's' or s eq 'S') then begin
        print,'discard Solid'
        oplot,xx,yy[*,nmerge] & oplot,xx,yy[*,jmerge],col=0
        print,'I am about to discard the marked cluster. '+ $
              'Are you sure (Y/N)? [N]'
        s=get_kbrd(1)
        if(s ne 'y' and s ne 'Y') then goto,bad_answer
        j=where(index eq nmerge+1)
        index[j]=0
      endif else if(s eq 'D' or s eq 'd') then begin
        print,'discard Dashed'
        oplot,xx,yy[*,nmerge],col=0 & oplot,xx,yy[*,jmerge]
        print,'I am about to discard the marked cluster. '+ $
              'Are you sure (Y/N)? [N]'
        s=get_kbrd(1)
        if(s ne 'y' and s ne 'Y') then goto,bad_answer
        j=where(index eq jmerge+1)
        index[j]=0
      endif else if(s eq 'C' or s eq 'c') then begin
        oplot,xx,yy[*,nmerge],col=0 & oplot,xx,yy[*,jmerge]
        wset,1
        !x=x_plot_new & !y=y_plot_new & !p=pnew
        print,'Click mouse button on the first cluster you want to merge'
        cursor,xcur_n,ycur_n,/DOWN,/DATA
        yo=dblarr(nregions-1)
        for kkkn=0,nregions-2 do yo[kkkn]= $
          poly((double(xcur_n)-xoffset[kkkn])*xscale[kkkn], orders[*,kkkn])
        yo=min(abs(ycur_n-yo),kkkn)
        oplot,xx,yy[*,kkkn],col=0,line=0,thick=3
        print,'Click on the cluster with which you want to merge'
        cursor,xcur_j,ycur_j,/DOWN,/DATA
        yo=dblarr(nregions-1)
        for kkkj=0,nregions-2 do yo[kkkj]= $
          poly((double(xcur_j)-xoffset[kkkj])*xscale[kkkj], orders[*,kkkj])
        yo=min(abs(ycur_j-yo),kkkj)
        oplot,xx,yy[*,kkkj],col=0,line=1,thick=3
        print,nmerge,kkkn,jmerge,kkkj
        nmerge=kkkn
        jmerge=kkkj
        j=where(index eq nmerge+1)
        index[j]=jmerge+1
        wset,0
        !x=x_plot_old & !y=y_plot_old & !p=pold
      endif else if(s eq 'I' or s eq 'i') then begin
        print,'Ignore'
        oplot,xx,yy[*,nmerge],col=0
        oplot,xx,yy[*,jmerge],col=0
        imerge[0,nmerge]=0
        skip=[[skip],[nmerge,jmerge]]
        goto,skip_pair
      endif else goto,bad_answer
      goto,next
    endif else begin
      j=where(index eq nmerge+1)
      index[j]=jmerge+1
      goto,next
    endelse
  endif
;  endif else begin ; Automatically discard order with high curvature
;    if(abs(orders[power,nmerge]) ge abs(orders[power,jmerge])) then begin
;      j=where(index eq nmerge+1)              ;Discard order nmerge
;      index[j]=0
;    endif else begin
;      j=where(index eq jmerge+1)              ;Discard order jmerge
;      index[j]=0
;    endelse
;  endelse

  for i=0,nregions-2 do begin  ; Discard orders that stick out in the middle
    jj=where(index eq i+1, njj)
    if(njj gt 0) then begin
      y1=min(yy[*,i],iy1)
      y2=max(yy[*,i],iy2)
      if(yy[   0,i] ge 0 and yy[   0,i] lt ny and $ ; left and right values are within
         yy[nx-1,i] ge 0 and yy[nx-1,i] lt ny and $ ; range but min or max are outside
         (y1 le 0 or  y2 ge ny)) then index[jj]=0
      if(y1 lt  0 and (yy[   0,i] ge ny or $        ; min<0 and one of the ends or 
                       yy[nx-1,i] ge ny or $        ; max are out of the roof
                      y2 ge ny)) then index[jj]=0
      if(y2 ge ny and (yy[   0,i] le 0 or $         ; max is outside the top boundary
                       yy[nx-1,i] le 0 or $         ; and one of the ends is below zero 
                      y1 le 0)) then index[jj]=0
    endif
  endfor

  yy=fltarr(nx,nregions-1)
  ord_err=dblarr(nregions-1)
  orders=dblarr(power+1,nregions-1)
  ind=0
  cl_limits=intarr(2,nregions-1)   ; starting and ending columns of clusters
  for ireg=1,nregions-1 do begin
    j=where(index eq ireg, nj)
    if(nj gt 0) then begin ; Check if cluster is not empty
      if(n_elements(uniq(x[j])) gt power) then begin ; Check if cluster includes
        index[j]=ind+1                               ; enough different pixels
        coef=poly_fit(x[j],y[j],power,/double)
                                                     ; Check if this cluster is
                                                     ; not a vertical artifact
        if(abs((poly_fit(x(j),y[j],1,/double))[1]) lt 4.*ny/nx) then begin
          cl_limits[*,ind]=minmax(x[j])
          orders[*,ind]=coef
          yy[*,ind]=poly(xx,orders[*,ind])
          ord_err[ind]=sqrt(total((y[j]-poly(x[j],orders[*,ind]))^2)/nj)
          ind=ind+1
        endif
      endif else index[j]=0
    endif
  endfor
  nregions=ind+1

; Kill curvy lines with derivatives way of the average
  if(nregions gt 5) then begin
    repeat begin
      orders=orders[*,0:nregions-2]
      for i=power,1,-1 do begin
        dev=stddev(orders(i,*))
        jgood=where(abs(orders[i,*]-mean(orders[i,*])) lt 4.*dev, ngood)
        if(ngood gt 0) then begin
          orders=orders[*,jgood]
          cl_limits=cl_limits[*,jgood]
          ord_err=ord_err[jgood]
          yy=yy[*,jgood]
        endif else begin
          print,'This is terribly wrong'
          stop
        endelse
      endfor
      nregions=n_elements(orders[0,*])+1
    endrep until (ngood eq nregions-1)
  endif

  if(keyword_set(iplot)) then wdelete,1

  ymiddle=fltarr(nregions-1)
  xmiddle=0.5d0*nx
  for i=0,nregions-2 do ymiddle[i]=poly(xmiddle,orders[*,i])
  i=sort(ymiddle)
  yy=yy[*,i]
  orders=orders[*,i]
  cl_limits=cl_limits[*,i]
  ord_err=ord_err[i]
  if(keyword_set(iplot)) then begin
    display,im,/log,min=(median(im)+2.*min(im))/3.
    for i=0,nregions-2 do oplot,xx,yy[*,i],col=0
  endif
  return
end

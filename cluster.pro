Function Cluster,x,y,nx,ny,NREGIONS=nregions,THRES=threshold,PLOT=iplot
;
; Cluster takes two 1D integer arrays with X and Y coodinates of
; pixels somehow selected in a 2D (Nx,Ny) array and identifies
; clusters of pixels. X and Y should be Y-sorted, that is X should run faster
; while Y should increase monotonously. If that is not the case, it could be
; achieved with the following commends:
; i=sort(Y*2L^15+X) & Y=Y(i) & X=X(i)
; A cluster is defined so that every pixel in a cluster is adjacent to at
; list one other pixel from the same cluster (their X and Y coordinates do
; not differ by more then 1) and all the pixels that are adjacent to any
; pixel in a cluster belong ; to the same cluster.
; The function returns a 1D integer array of the same size as X and Y
; that contains cluster number for each pixel. Non-cluster members
; are marked with 0.
; Optional parameters:
;   nregions - (output) contain the number of identified clusters also
;              counting non-cluster pixels (if any) as a separate cluster
;   thres    - (input, default: 1) consider clusters with "threshold"
;              or smaller number of members to be non-clustered (index eq 0)
;   plot     - (input) if set color pixels as clusters are identified
;              (a bit slower), otherwise print % of completion
;
; History: 17-July-2000 N.Piskunov wrote the version optimized for
;                       clusters oriented preferentially along rows or
;                       columns.
;          21-July-2000 N.Piskunov modified to handle arbitrary oriented
;                       clusters in the optimal way. It is slower than the
;                       original version by about 10% for spectral orders
;                       which are nearly horizontal/vertical.
;
  if(n_params() lt 4) then begin
    print,'Syntax: ind=Cluster(x,y,nx,ny[,NREGIONS=n[,THRES=t[,/PLOT]]]'
    print,'where: x and y   are 1D arrays with column and row numbers of pixels'
    print,'                 to be analized for clustering,'
    print,'       nx and ny are the dimensions of the original image,'
    print,'       n         is the number of detected clusters+1'
    print,'       t         is the threshold for the size of a cluster'
    print,'       /PLOT     displays clusters as they are detected'
    print,'Cluster returns an array of the same size as x and y. The value'
    print,'of index associates pixel with the cluster number. Pixels not'
    print,'associated with clusters are marked with index 0.'
    return,0L
  endif
  
  if(n_elements(x) ne n_elements(y)) then begin
    print,'Cluster: X and Y arrays should have the same size'
    return,0L
  endif

  if(max(x) ge nx) then begin
    print,'Cluster: X cannot exceed NX-1'
    return,0L
  endif

  if(max(y) ge ny) then begin
    print,'Cluster: Y cannot exceed NY-1'
    return,0L
  endif

  if(keyword_set(iplot)) then plot,x,y,xs=3,ys=3,psym=3
  n=n_elements(x)                  ; Number of pixels in objects
  index=intarr(n)                  ; Created index array

                                     ; The initial X and Y are Y-sorted, which
                                     ; means that X runs faster while Y is
                                     ; increasing monotonously.
  ix_sort=sort( x*2L^15+y)           ; Pointers to X-sorted arrays
  iy_sort=sort((y*2L^15+x)(ix_sort)) ; Pointers to Y-sorted array from X-sorted
                                     ; arrays. This means that:
                                     ; x(ix_sort(iy_sort)) eq x and
                                     ; y(ix_sort(iy_sort)) eq y
  ; The double index cunstructed above allows us to find the nearest pixels
  ; directly without using "where" which is costly in large arrays.

  max_branch=1                     ; Number of found branches
  branches=[0L,0L]                 ; First/last pixel of each branch

  iper=0
  for i=0L,n-1L do begin           ; Loop through pixels
    if(not keyword_set(iplot) and i*10L/n gt iper) then begin
      iper=i*10L/n                 ; % complete
      print,iper*10,'% done'
    endif

    j=i                            ; Find pixels nearest to i in the same column
    j1=iy_sort(i)-1L
    j2=iy_sort(i)+1L
    if(j1 ge   0L) then j=[ix_sort(j1),j]
    if(j2 lt n-1L) then j=[j,ix_sort(j2)]

    if(i gt 0L)   then begin
      j=[i-1L,j]                   ; Find pixels nearest to i in the previous column
      j1=iy_sort(i-1L)-1L
      j2=iy_sort(i-1L)+1L
      if(j1 ge   0L) then j=[ix_sort(j1),j]
      if(j2 lt n-1L) then j=[j,ix_sort(j2)]
    endif
    if(i lt n-1L) then begin
      j=[j,i+1L]                   ; Find pixels nearest to i in the next column
      j1=iy_sort(i+1L)-1L
      j2=iy_sort(i+1L)+1L
      if(j1 ge   0L) then j=[ix_sort(j1),j]
      if(j2 lt n-1L) then j=[j,ix_sort(j2)]
    endif
    j=j(sort(j))

   ; Find immediate neighours. Where searches in 9 pixels or less.
    j=j(where(x(j) ge x(i)-1 and x(j) le x(i)+1 and $
              y(j) ge y(i)-1 and y(j) le y(i)+1, nj))
    nmax=max(j)

    if(nj gt 1) then begin         ; Check if this pixel has neighbours
      jj=where(index(j) gt 0, njj) ; Check for existing branches
      if(njj eq 0) then begin      ; None found (only unmarked pixels)
        free_branch=where(branches(0,*) eq -1, nfree)
        if(nfree gt 0) then begin
          free_branch=free_branch(0)
          branches(*,free_branch)=minmax(j)
          index(j)=free_branch
        endif else begin
          branches=[[branches],[minmax(j)]] ; First/last pixel of the current branch
          index(j)=max_branch        ; Start a new branch
          max_branch=max_branch+1    ; Increment branch number
        endelse
      endif else begin             ; This is part of the existing branch
        curr_index=index(j(jj))
        if(min(curr_index) eq max(curr_index)) then begin
          ind=curr_index(0)
          index(j)=ind
          branches(0,ind)=j(0)   <branches(0,ind)
          branches(1,ind)=j(nj-1)>branches(1,ind)
        endif else begin
          curr_branch=min(curr_index)
          curr_indices=curr_index(sort(curr_index))
          curr_indices=curr_indices(uniq(curr_indices))
; Paint all relevant branches with the same index
          for jjj=1,n_elements(curr_indices)-1 do begin
            ind=curr_indices(jjj)              ; Next non-zero index
            i1=branches(0,ind)                 ; First pixel of the ind branch
            i2=branches(1,ind)                 ; Last pixel of the ind branch
            ii=where(index(i1:i2) eq ind)+i1   ; Find all pixels in this branch
            index(ii)=curr_branch              ; and merge branches if needed
            branches(0,curr_branch)=i1<branches(0,curr_branch)
            branches(1,curr_branch)=i2>branches(1,curr_branch)
          endfor
; Remove the branches that were merged to a lower index (curr_branch)
          for jjj=1,n_elements(curr_indices)-1 do begin
            ind=curr_indices(jjj)
            branches(*,ind)=[-1,-1]
          endfor
; Index unmarked pixels and verify they do not change the pointers
          index(j)=curr_branch     ; Mark immediate neighbours
          branches(0,curr_branch)=j(0)   <branches(0,curr_branch)
          branches(1,curr_branch)=j(nj-1)>branches(1,curr_branch)
        endelse
      endelse
      if(keyword_set(iplot)) then oplot,x(j),y(j),psym=3,col=(index(i)+1) mod 128
    endif
  endfor
  if(not keyword_set(iplot)) then print,100,'% done'

  if(keyword_set(iplot)) then plot,x,y,xs=3,ys=3,psym=3,/nodata
  ind=0
  if(not keyword_set(threshold)) then threshold=1L
  for i=0,max_branch-1 do begin
    i1=branches(0,i)                     ; First pixel of the ind branch
    i2=branches(1,i)                     ; Last pixel of the ind branch
    if(i1 ge 0) then begin
      j=where(index(i1:i2) eq i, nj)+i1
      if(nj gt threshold) then begin
        ind=ind+1
        index(j)=ind
        if(keyword_set(iplot)) then oplot,x(j),y(j),psym=3,col=(ind mod 128)
      endif else if(nj gt 0) then begin
        index(j)=0
      endif
    endif
  endfor
  nregions=ind+1

  return,index
end

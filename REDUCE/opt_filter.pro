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

Function Opt_Filter,y,par,par1,DOUBLE=dbl,MAXIT=maxiter $
                   ,WEIGHT=weight,LAM2=lam2
;
; Optimal filtering of 1D and 2D arrays.
; Uses tridiag in 1D case and sprsin and linbcg in 2D case.
; Written by N.Piskunov 8-May-2000
;
common bandsolv,band_solv_name

  if(n_params() lt 2) then begin
    print,'Optimal filtering routine:'
    print,'Syntax: r=opt_filter(f,xwidth[,ywidth[,WEIGHT=weight[,/DOUBLE[,MAXIT=maxiter]]]])'
    print,'where:  f      - 1D or 2D array of type I,F or D'
    print,'        xwidth - filter width (for 2D array width in X direction (1st index)'
    print,'        ywidth - (for 2D array only) filter width in Y direction (2nd index)'
    print,'                 if ywidth is missing for 2D array, it set equal to xwidth'
    print,'        weight - an array of the same size(s) as f containing values between 0 and 1'
    print,'        DOUBLE - perform calculations in double precision'
    print,'        maxiter- maximum number of iteration for filtering of 2D array'
    print,'        weight - weight for the function (values between 0 and 1)'
    print,'   opt_filter solves the optimization problem for r:'
    print,'        total(weight*(f - r)^2) + width*total((r(i) - r(i-1))^2) = min'
    if(keyword_set(y)) then return,y else return,0
  endif

  if(not keyword_set(band_solv_name)) then begin
    help,calls=a
    a=strmid(a(0),strpos(a(0),'<')+1,strlen(a(0))-strpos(a(0),'<')-1)
    if(strpos(a,'/nfs') eq 0) then a=strmid(a,4)
    if(strpos(a,'/export') eq 0) then a=strmid(a,7)
    if(!VERSION.OS eq 'Win32') then delimiter='\' else delimiter='/'
    prefix=strmid(a,0,strpos(a,delimiter,/REVERSE_SEARCH))   ;slit_func directory
    if(strpos(a,delimiter) lt 0) then cd,CURRENT=prefix
    band_solv_name=prefix+delimiter+'bandsol.so.'   $
                                   +!version.os+'.' $
                                   +!version.arch+'.' $
                                   +strtrim(!version.MEMORY_BITS,2)
  endif

  if((size(y))(0) eq 1 or $
    ((size(y))(0) eq 2 and ((size(y))(1) eq n_elements(y) or $
                            (size(y))(2) eq n_elements(y)))) then begin
    if(par lt 0) then return,y
    n=n_elements(y)
    if(keyword_set(dbl)) then begin
      if(keyword_set(weight)) then wgt=double(weight) else wgt=replicate(1.d0,n)
      if(keyword_set(lam2)) then lambda2=double(lam2) else lambda2=-1.d0
      if(lambda2 gt 0) then begin
        aij=dblarr(n,5)
;
; 2nd lower subdiagonal
        aij[2:n-1,0]=lambda2
;
; Lower subdiagonal
        aij[1    ,1]=-par-2.d0*lambda2
        aij[2:n-2,1]=-par-4.d0*lambda2
        aij[n-1  ,1]=-par-2.d0*lambda2
;
; Main diagonal
        aij[0    ,2]=wgt[0    ]+par+lambda2
        aij[1    ,2]=wgt[1    ]+2.d0*par+5.d0*lambda2
        aij[2:n-3,2]=wgt[2:n-3]+2.d0*par+6.d0*lambda2
        aij[n-2  ,2]=wgt[n-2  ]+2.d0*par+5.d0*lambda2
        aij[n-1  ,2]=wgt[n-1  ]+par+lambda2
;
; Upper subdiagonal
        aij[0    ,3]=-par-2.d0*lambda2
        aij[1:n-3,3]=-par-4.d0*lambda2
        aij[n-2  ,3]=-par-2.d0*lambda2
;
; 2nd lower subdiagonal
        aij[0:n-3,4]=lambda2
;
; RHS
        b=double(reform(wgt*y))

        i=CALL_EXTERNAL(band_solv_name, 'bandsol', $
                        aij, b, long(n), 5L)
        f=b
      endif else begin
        a=replicate(-double(abs(par)),n)
        b=[wgt[0]+abs(par),wgt[1:n-2]+replicate(2.d0*abs(par),n-2),wgt[n-1]+abs(par)]
        f = trisol(a, b, a, double(reform(wgt*y)), DOUBLE=dbl)
      endelse
    endif else begin
      if(keyword_set(weight)) then wgt=weight else wgt=replicate(1.0,n)
      if(keyword_set(lam2)) then lambda2=lam2 else lambda2=-1.0
      if(lambda2 gt 0) then begin
        aij=dblarr(n,5)
;
; 2nd lower subdiagonal
        aij[2:n-1,0]=lambda2
;
; Lower subdiagonal
        aij[1    ,1]=-par-2.d0*lambda2
        aij[2:n-2,1]=-par-4.d0*lambda2
        aij[n-1  ,1]=-par-2.d0*lambda2
;
; Main diagonal
        aij[0    ,2]=wgt[0    ]+par+lambda2
        aij[1    ,2]=wgt[1    ]+2.d0*par+5.d0*lambda2
        aij[2:n-3,2]=wgt[2:n-3]+2.d0*par+6.d0*lambda2
        aij[n-2  ,2]=wgt[n-2  ]+2.d0*par+5.d0*lambda2
        aij[n-1  ,2]=wgt[n-1  ]+par+lambda2
;
; Upper subdiagonal
        aij[0    ,3]=-par-2.d0*lambda2
        aij[1:n-3,3]=-par-4.d0*lambda2
        aij[n-2  ,3]=-par-2.d0*lambda2
;
; 2nd lower subdiagonal
        aij[0:n-3,4]=lambda2
;
; RHS
        b=double(reform(wgt*y))
        i=CALL_EXTERNAL(band_solv_name, 'bandsol', $
                        aij, b, long(n), 5L)
        f=b
      endif else begin
        a=replicate(-double(abs(par)),n)
        c=a
        b=[wgt[0]+abs(par),wgt[1:n-2]+replicate(2.0*abs(par),n-2),wgt[n-1]+abs(par)]
        f = trisol(a, b, c, double(reform(wgt*y)), DOUBLE=dbl)
      endelse
    endelse
    return,f
  endif else if((size(y))(0) eq 2) then begin
    if(not keyword_set(par1)) then par1=par
    if(par eq 0 and par1 eq 0) then return,y
    n=n_elements(y)
    nc=(size(y))(1)
    nr=(size(y))(2)

    if(keyword_set(dbl)) then begin
      adiag=double(abs(par ))
      bdiag=double(abs(par1))
    endif else begin
      adiag=float(abs(par ))
      bdiag=float(abs(par1))
    endelse
                                                    ;Main diagonal first:
    aa=[1.+adiag+bdiag,                           $ ; upper-left corner
        replicate(1.+2.*adiag+   bdiag,nc-2),     $ ; first row
        1.+adiag+bdiag,                           $ ; upper-right corner
        replicate(1.+2.*adiag+2.*bdiag,n-2*nc),   $ ; all other points
        1.+adiag+bdiag,                           $ ; lower-left corner
        replicate(1.+2.*adiag+   bdiag,nc-2),     $ ; last row
        1.+adiag+bdiag,                           $ ; lower-right corner
        replicate(-adiag,n-1),                    $ ;Upper sub-diagonal for X
        replicate(-adiag,n-1),                    $ ;Lower sub-diagonal for X
        replicate(-bdiag,n-nc),                   $ ;Upper sub-diagonal for Y
        replicate(-bdiag,n-nc)]                     ;Lower sub-diagonal for Y

    col=lindgen(nr-2)*nc+nc                         ;Special cases:
    aaa=replicate(1.+adiag+2.*bdiag,nr-2)
    aa(col     )=aaa                                ; last columns
    aa(col+nc-1)=aaa                                ; first column
    col=n+lindgen(nr-1)*nc+nc-1
    aa(col    )=0.
    aa(col+n-1)=0.

    col=[lindgen(n),           $                    ;Main diagonal
         lindgen(n-1)+1,       $                    ;Upper sub-diagonal for X
         lindgen(n-1),         $                    ;Lower sub-diagonal for X
         lindgen(n-nc)+nc,     $                    ;Upper sub-diagonal for Y
         lindgen(n-nc)]                             ;Lower sub-diagonal for Y

    row=[lindgen(n),           $                    ;Main diagonal
         lindgen(n-1),         $                    ;Upper sub-diagonal for X
         lindgen(n-1)+1,       $                    ;Lower sub-diagonal for X
         lindgen(n-nc),        $                    ;Upper sub-diagonal for Y
         lindgen(n-nc)+nc]                          ;Lower sub-diagonal for Y

    aaa=sprsin(col,row,aa,n,DOUBLE=dbl,THRESH=-2.*(adiag>bdiag))
    col=bdiag
    row=adiag
    aa=reform(y,n)                ;Start with an initial guess at the solution.

    if(keyword_set(maxiter)) then maxit=maxiter else maxit=50
    aaa=linbcg(aaa,reform(y,n),aa,ITMAX=maxit,DOUBLE=dbl);Solve the linear system Ax=b.
    aaa=reform(aaa,nc,nr)                        ;Restore the shape of the result.
    return,aaa
  endif else begin
    return,y
  endelse
end

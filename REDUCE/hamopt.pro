pro HamOpt, im, back, ycen, yslitf, slitf, bincen, s, u,  $
            ccd_gain, ccd_dark, ccd_rdnoise, $
            x_left_lim=cole0, x_right_lim=cole1
;
;  Optimal extraction routine for the Hamilton reduction package.
;  Inputs:
;     im - The current image being reduced.
;     back - Background array for the current *order* being extracted.
;            Background is stored in (ncol,2) array where back(*,0) is the
;            value in the center of the order and bacl(*,1) is the slop
;            per row.
;     ycen - A vector giving the center of the order for each column.
;     yslitf - Vector giving the oversampled relative pixel scale for slitf.
;     slitf - An array giving the oversampled slit function at a number of
;               swaths in the image.
;     bincen - The column of the bin centers for the slit functions.
;  Outputs:
;     s - The optimally extracted spectrum.
;     u - The uncertainty in s.
;
; 31-Mar-1998  CMJ, JAV  Written
; 26-Jan-2000  NP, removed common ham_aux, replaced with data from
;              inst_setup structure available in ham.common
;  9-Oct-2000  NP, added processing of 2D background
;

  if n_params() lt 11 then begin
     print,'Syntax is: HamOpt,im,back,ycen,yslitf,slitf,bincen,s,u'
     retall
  endif

  ncol = (size(im))(1)				;get number of columns in image
  nrow = (size(im))(2)				;get number of rows in image
  nbin = n_elements(bincen)         ;get number of bins
  nysf = n_elements(yslitf)         ;size of oversampled slitf
  if(not keyword_set(cole0)) then cole0 = 0      ;first column to extract
  if(not keyword_set(cole1)) then cole1 = ncol-1 ;last column to extract
  ncole = cole1 - cole0 +1          ;number of columns to extract
  s = fltarr(ncole)			  	    ;init spectrum
  u = fltarr(ncole)				    ;init uncertainty vector

; Linearly interpolate slitf onto columns.
  bincol = -0.5+findgen(ncole)/(ncole-1.)*nbin + cole0 ;column number in bin centers
  sf = interpolate(slitf, findgen(nysf) $
                        , bincol, /grid)        ;slit function at each column
  osamp = (nysf-1.)/(yslitf(nysf-1)-yslitf(0))  ;calculate oversampling (pix)

; Now loop through columns and optimally extract spectra.
  ysfmin = min(yslitf, max=ysfmax)
  pbg = sqrt(ccd_dark + ccd_rdnoise^2.)
  for i = cole0, cole1 do begin
    ibeg = ceil(ycen(i)+ysfmin) > 0
    iend = floor(ycen(i)+ysfmax) < (nrow - 1)
    bgg = back(i,0) + back(i,1)*(findgen(iend - ibeg + 1) + ysfmin)
;    data = reform(im(i, ibeg:iend) - back(i))
    data = reform(im(i, ibeg:iend)) - bgg
    idata = findgen(iend - ibeg + 1) + ibeg - ycen(i)
    modl = interpolate(sf(*,i-cole0), (idata-ysfmin)*osamp)
    hamopt_c, data, modl, ccd_gain, pbg, int, sint, tot, stot
    s(i-cole0) = int
    u(i-cole0) = sint
  endfor
end

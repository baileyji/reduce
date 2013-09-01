pro sumbias, list1, list2, inst_setting, bias, head, DEBUG=debug $
           , XR=xr, YR=yr, ERR=err, EXTEN=exten
;Combine bias frames, determine read noise, reject bad pixels.
;Read noise calculation only valid if both lists yield similar noise.
;
;30-Aug-2000 NP Removed all traces of the instrument-specific behavour, modified
;               to use the new sumfits routine and included new error handling
;               via string ERR parameter.
;16-Nov-2005 NP Added extension number keyword to be passed to sumfits.
;               This way sumfits can handle fiels with multiple images.

  err = ''
  if n_params() lt 4 then begin
    print, 'syntax: sumbias, list1, list2, inst_setting, bias[, head[, /DEBUG[, ERR=err]]]'
    err = '%SUMBIAS: Insufficient number of parameters'
    retall
  endif

;Lists of images.
  list = [ list1, list2 ]
  n1 = n_elements(list1)
  n2 = n_elements(list2)
  n  = n1 + n2

;Separately images in two groups.

  sumfits, list1, inst_setting, bias1, head1, XR=xr, YR=yr, ERR=err, EXTEN=exten
  if(keyword_set(err)) then begin
    err = '%SUMBIAS: ERROR ' + err
    return
  endif
  bias1 = clipnflip(bias1 / float(n1), head1)

  sumfits, list2, inst_setting, bias2, head2, XR=xr, YR=yr, ERR=err, EXTEN=exten
  if(keyword_set(err)) then begin
    err = '%SUMBIAS: ERROR ' + err
    return
  endif
  bias2 = clipnflip(bias2 / float(n2), head2)

  if((size(bias1))(0) ne (size(bias2))(0) or(size(bias1))(1) ne (size(bias2))(1)) then begin
    err = '%SUMBIAS: Bias frames in two lists have different dimensions'
    return
  endif

;Make sure we know the gain.
  head = head2
  gain = sxpar(head, 'E_GAIN*')
  if(n_elements(gain) gt 1) then gain = 1

;Construct unnormalized sum.
  bias  = bias1 * float(n1) + bias2 * float(n2)

;Normalize the sum.
  bias  = bias  / float(n)

;Compute noise in difference image by fitting Gaussian to distribution.
  diff = 0.5 * (bias1 - bias2)              ;0.5 like the mean...
  if(min(diff) ne max(diff)) then begin
    crude = median(abs(diff))                 ;estimate of noise
    hmin = -5.0 * crude
    hmax = +5.0 * crude
    bsiz = 2.0 / float(n)
    h = histogram(diff, min=hmin, max=hmax, bin=bsiz)
;  nh = ceil((hmax - hmin) / bsiz)          ;# histogram points
    nh = n_elements(h)                        ;# histogram points
    xh = hmin + bsiz * (dindgen(nh) + 0.5)
    hfit = gaussfit(xh, h, par, nterm=3)
    noise = abs(par(2))                       ;noise in diff, bias

;Determine where wings of distribution become significantly non-Gaussian.
    contam = (h - hfit) / sqrt(hfit>1.)
    imid = where(abs(xh) lt 2*noise)
    consig = stdev(contam(imid))
    smcontam = gaussbroad(xh, contam, 0.1*noise)
    thresh = 3.0 * consig
    igood = where(smcontam lt thresh)
    gmin = min(xh(igood), max=gmax)

;Find and fix bad pixels.
    ibad = where(diff le gmin or diff ge gmax, nbad)

    diff=0

;stop
;    bias = (float(n1) / float(n)) * bias1 + (float(n2) / float(n)) * bias2
    if(nbad gt 0) then bias(ibad) = bias1(ibad) < bias2(ibad)

;Compute read noise.
    biasnoise = gain * noise
    bgnoise = biasnoise * sqrt(n)

;Print diagnostics.
    print, 'Change in bias between image sets= ' $
         + strtrim(string(gain*par(1), form='(f9.2)'), 2) + ' electrons'
    print, 'Measured background noise per image= ' $
         + strtrim(string(bgnoise, form='(f9.2)'), 2) + ' electrons/pixel'
    print, 'Background noise in combined image= ' $
         + strtrim(string(biasnoise, form='(f9.2)'), 2) + ' electrons/pixel'
    print, 'Fixing ' + strtrim(nbad, 2) + ' bad pixels'

;Plot noise distribution.
    if keyword_set(debug) then begin
      !p.multi = [0, 1, 2]
      plot, xh, h, xsty=3, ps=10,title='Noise distribution'
      colors
      oplot, xh, hfit, co=2
      oplot, gmin+[0,0], !y.crange, co=6
      oplot, gmax+[0,0], !y.crange, co=6

;Plot contamination estimation.
      plot, xh, contam, xsty=3, ps=10,title='Contamination estimation'
      oplot, !x.crange, thresh+[0,0], co=2
      oplot, xh, smcontam, co=3
      oplot, gmin+[0,0], !y.crange, co=6
      oplot, gmax+[0,0], !y.crange, co=6
      !p.multi = 0
      print, form='(a,$)', 'Push space to continue...'
      junk = get_kbrd(1)
      print, ''
    endif
  endif else begin
    diff=0
    biasnoise=1.
    nbad=0L
  endelse

;Toss blank header lines.
  iwhr = where(strtrim(head, 2) ne '')
  head = head(iwhr)

  obslist = strtrim(list[0], 2)
  for i=1, n_elements(list)-1 do $
    obslist = obslist + ' ' + strtrim(list[i], 2)
  sxdelpar, head, 'TAPELIST'
  sxaddpar, head, 'BZERO', 0.0
  sxaddpar, head, 'BSCALE', 1.0
  sxaddpar, head, 'OBSLIST', obslist
  sxaddpar, head, 'NIMAGES', n $
          , ' number of images summed'
  sxaddpar, head, 'NPIXFIX', nbad $
          , ' pixels corrected for cosmic rays'
  sxaddpar, head, 'BGNOISE', biasnoise $
          , ' noise in combined image, electrons'
  return

end

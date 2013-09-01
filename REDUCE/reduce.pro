;This is a generic script to run reduction of echelle spectra using REDUCE
;package. This package was written by Jeff Valenti (STScI) and
;Nikolai Piskunov (Uppsala University). Read "New algorithms for reducing
;cross-dispersed echelle spectra", N.Piskunov, J.Valenti 2002 A&A 385, 1095.
;It is based on HAMSPEC package written by JV for Hamilton spectrograph,
;on 2 mathematical algorithms proposed and implemented by NP (cluster
;analysis and spectral order decomposition), and on 2+1D wavelength
;calibration tool which is a joint effort by JV and NP.
;
;WARNING: REDUCE was never ment to be a complete, user friendly package
;         but rather a demonstration platform for a few new algorithms.
;
;Current status: 2004-09-30
; REDUCE was successfully used for reducing the data from
; the following instruments (see modeinfo.pro):
;    - ESO UVES (both arms and all CCDs, all slicers, FLAMES)
;    - ESO CES (slit, all slicers)
;    - McDonald 2D
;    - Zelentchuk NES (slit, slicer, polarimeter)
;    - SAAO Giraffe
;    - Elodie
;    - FEROS
;    - SARG (spectroscopy and spectropolarimetry)
;    - S. Korean BOES
;    - Keck HiRes
; 
;REDUCE consists of a number of rather independent subroutines implementing
;different steps of reduction. All information is passed via mandatory and
;optional parameters (no common blocks) but to keep things under control
;recommend certain convention in file naming and data processing sequence.
;This is maintained by the control script and here is an example of such
;script. Note that it does not include wavelength calibration which is a
;separate widget-based program.
;
;This sample script is tuned for reducing UVES data ("uves_middle" stands
;for the "blue" CCD in the red arm of UVES)
;
;All the raw frames after reading are going through a standard procedure:
;
; a) header is modified by a call to modeinfo. New keywords start with "E_".
; b) frame is cropped and rotated into standard REDUCE orientation with
;    clipnflip (orders run horizontally etc. See modeinfo for details)
;
;This procedure is also performed over master bias and master flat but
;in the first case it is hidden inside sumbias while sumfits only contains
;step a)
;
;=======================New UVES data format================================
;In 2004 ESO has changed the format of the data obtained with the red arm
;of UVES. Now the "middle" and the "red" CCDs are stored in a single FITS
;file with two extensions: Ext #1 - red; Ext #2 - middle.
;This would have no consequences for REDUCE except for sumfits which is also
;called from sumbias. sumfits reads FITS files gradually and thus requires
;to know which extension it should read. The cheap solution is to split all
;the files form red arm into two and treat them the old way. This can be done
;with a following simple loop:
;
;  file_list=file_search('*.fits.gz')      ;Get the list of files for red arm
;  for i=0,n_elements(file_list)-1 do begin
;    fits_read,file_list(i),im,h,exten=1   ;Read ext. #1
;    writefits,'RED/'+file_list(i),im,h    ;Store it in subdirectory RED
;    fits_read,file_list(i),im,h,exten=2   ;Read ext. #2
;    writefits,'MIDDLE/'+file_list(i),im,h ;Store it in subdirectory MIDDLE
;  endfor
;
;In parallel NP included modifications to sumfits and sumbias which can handle
;the new format directly but no thorough testing was performed so ... you
;may find small problems. Also, the modeinfo does not yet distinguish between
;the two formats and does not supply the correct extension number. It must
;be put manually


  inst_mode = 'uves_red' ;Spectrometer identification to modeinfo
  exten=1
  prefix  = 'Reduced_Red/' ;Prefix for intermediate products (e.g. master flat)
  path    = './'     ;Top directory for the data
  
;Automatic file classification. Alternatively one can put bias, flats, ThAr and
;science frames in separate subdirectories
  file_list=file_search(path+'*.fits.gz') ;Get file list
                                  ;REDUCE can take gzipped FITS directly
  file_type=file_list
  instrum=file_list
  for i=0,n_elements(file_list)-1 do begin;Find object description
    h=headfits(file_list[i])
    file_type[i]=hierarch(h, 'HIERARCH ESO DPR TYPE')
  endfor

;Compile a list of BIAS frames
  biaslist=file_list[where(strpos(file_type,'BIAS') ge 0)]
  bias_file = prefix+inst_mode+'.bias.fits'         ;Master bias filename

;Compile a list of FF frames
  flatlist=file_list[where(strpos(file_type,'LAMP,FLAT') ge 0)]
  flat_file=prefix+inst_mode+'.flat.fits'           ;Master flat filename
  norm_flat_file=prefix+inst_mode+'.flat.norm.fits' ;Normalized master flat filename

;Order definition file
  dord = file_list[where(strpos(file_type,'LAMP,ORDERDEF') ge 0)]

;Compile a list of ThAr frames
  tharlist=file_list[where(strpos(file_type,'LAMP,WAVE') ge 0)]

;Assume that all remaining frames are science data
  spec=file_list[where(file_type eq 'OBJECT,POINT')]

  mask_file='~/IDL/REDUCE/mask_uves_red.fits.gz';Mask of cosmetic defects on the CCD
;Mask structure: it is a FITS file with the same sizes as the raw CCD frames,
;byte-type pixels (8-bit) and values of 1B for good pixels and 0B for bad
;onces. Mask helps handling situations when extended defects are alsmost
;parallel to spectral orders. If you don't have a mask and do not want to
;make one just keep mask_file undefined and comment out the corresponding
;read/clipnflip statements below.

  opower=4      ; Power of the order fitting polynomial

;goto,skip  ;simple trick to re-run parts of the reduction without repeating
           ; the whole thing. Just put label skip: where you want to restart

;Compile and save master bias
;Master bias by splitting the list of bias frames in two groups (e.g. before
;and after the the night). In each group median bias is produced using sumfits
;routine. The two groups are then compared to each other to identify possible
;trends. The net result is the mean (master) bias with most of outlayers
;removed 
  nbias=n_elements(biaslist)
  sumbias,biaslist[0:nbias/2-1],biaslist[nbias/2:nbias-1],inst_mode $
         ,bias,bhead,XR=50+[0,2047],ERR=err,EXTEN=exten
  if(keyword_set(err)) then stop
  writefits,bias_file,bias,bhead  ;Write out master bias

;Compile and save master flat
;Master flat is done for one group is the internal differences in counts
;between individual frames may indicate exposure differences. Note that
;sumfits actually SUMS signals in all frames, therefore, we sut divide by the
;number of flats before subtracting master bias.
;skip:
  nflats=n_elements(flatlist)
  sumfits,flatlist,inst_mode,flat,fhead,XR=50+[0,2047],ERR=err,EXTEN=exten
  if(keyword_set(err)) then stop
  flat=clipnflip(flat, fhead)
  flat=flat/nflats
  sz1=size(bias)
  sz2=size(flat)
  if(sz1(0) ne sz2(0) or sz1(1) ne sz2(1) or sz1(2) ne sz2(2)) then begin
    print,' BIAS frames should have the same dimensions as the FF frames'
    stop
  endif

  bias=readfits(bias_file,bhead)
  flat=(flat-bias)>1.             ;Subtract master bias
  writefits,flat_file,flat,fhead  ;Write out master flat
;  stop

  colors

;==========================================================================
;skip:
;Trace order locations and find default extraction windows.
  dord=dord[0]
  fits_read,dord,ordim,hdr,exten=exten ; Read order definition frame NEW UVES format
; ordim=readfits(dord, hdr)            ; Read order definition frame OLD UVES format
  hdr=modeinfo(hdr, inst_mode, XR=50+[0,2047])
  ordim=clipnflip(ordim, hdr)
  mask=readfits(mask_file,mhead)   ;Get mask if any
  mhead=modeinfo(mhead, inst_mode)
  mask=clipnflip(mask, mhead)

;Locate spectral orders
;Parameter(s) to be tuned: filter - unsharp filter parameter in
;                                   cross-dispersion direction
  hamdord,ordim,orders,or_err=ord_err,or_range=or_range,MASK=mask, $
          or_column=col_range,filter=200.,power=opower,/PLOT

;Determine extraction windows above and below order maximum.
  getxwd,ordim,orders,def_xwd,def_sxwd,colrange=col_range ;get extraction width

;Save image format description for all orders
  save,orders,or_range,ord_err,col_range,def_xwd,def_sxwd $
      ,file=prefix+'.ord.sav'

  ordim = 0

;==========================================================================
;Construct normalized flat field.
;skip:
  erase
  device,decomposed=0

  restore,prefix+'.ord.sav'
;  orders[0,*]=orders[0,*]-35

  flat=readfits(flat_file,fhead)   ;Get master flat
  mask=readfits(mask_file,mhead)   ;Get mask if any
  mhead=modeinfo(mhead, inst_mode)
  mask=clipnflip(mask, mhead)
;
;Each sp. order in master flat is devided in equal vertical swaths of 
;selected width. Each such swath is decomposed into slit illumination
;function and "spectrum" (see REDUCE paper) with additional smoothing
;applied. Smoothing in the slit direction is needed for proper handling
;of low signal parts while smoothing in the dispersion direction defines
;division flat field features between the normalized flat and the blaze
;function ("spectrum"). The normalized flat is made by dividing each column
;in each swath with the corresponding value of a blaze function.
;Ideally one would like to see the spectrum to hold only the blaze function
;while the normalized flat should have pixel-to-pixel variations.
;In practice, one has to adjust this empirically for each instrument setting
;as on encounters such things as fringes with relatively low spatial
;frequencies but which cannot be reproduced by the 1D blaze function.
;Filter parameter in dispersion direction affects the division of signal
;between the 1D blaze and 2D normalized flat in spatial frequencies: higher
;value restricts blaze to lower frequency components.

;Normalized flat is set to 1 where S/N is less then 100
;Parameter(s) to be tuned: swath_width - swath width in columns
;                          SF_SMOOTH - smoothing accross dispersion
;                          SP_SMOOTH - smoothing in dispersion direction
;                          OSAMPLE   - slit function is reconstructed on
;                                      subpixel grid with stepsize OSAMPLE
;                                      times smaller than CCD pixels. Larger
;                                      OSAMPLE require more computing time and
;                                      larger SF_SMOOTH. Plots of residuals
;                                      allow to verify if the OSAMPLE is good
;                                      enough.
  hamflat, flat, fhead, orders, blzcoef, colrange=col_range, FXWD=def_xwd $
         , SF_SMOOTH=20., SP_SMOOTH=20., OSAMPLE=12, swath_width=400 $
         , MASK=mask, /PLOT
  writefits,norm_flat_file,flat,fhead ;Write out normalized flat

;Add blaze functions to the save file
  save,orders,or_range,ord_err,col_range,def_xwd,def_sxwd,blzcoef $
      ,file=prefix+'.ord.sav'
  mask=0
  flat=0
  stop

;==========================================================================
;Extract thorium spectra.
;skip:
  restore,prefix+'.ord.sav'

;  tharlist=spec
  for n=0,n_elements(tharlist)-1 do begin
    fits_read,tharlist[n],im,head,exten=exten ; Read ThAr NEW UVES format
;   im = readfits(tharlist[n], head)          ; Read ThAr OLD UVES format
    head = modeinfo(head, inst_mode, XR=50+[0,2047]) ; Modify header with instrument setup
    im=clipnflip(im, head)                     ; Clip & flip

;Extraction of ThAr. Non-optimal extraction is used, so no parameters to tune.
    hamspec,im,head,orders,def_xwd,def_sxwd,or_range(0),thar $
           ,sig=sunc,colrange=col_range,/thar
;    display,im,/log
;    for i=0,n_elements(orders(0,*))-1 do $
;      oplot,poly(dindgen(4096),orders(*,i)),col=0
;    im=0

    sxaddpar, head, 'OBASE', or_range[0], ' base order number'
;    fileout=strmid(tharlist[n],strlen(path))
    fileout=tharlist[n]
    fileout=strmid(fileout,0,strpos(fileout,'.',/REVERSE_SEARCH)) $
           +'.thar.ech'
    fileout=prefix+fileout
    wdech,fileout,head,thar,/overwrite    ;save spectrum to disk
    print,'Extracted ThAr was written to: "'+fileout+'"'
  endfor
  stop

;==========================================================================
;Extract science spectra.
skip:
  restore,prefix+'.ord.sav'
  nord=n_elements(orders[0,*])

  bias=readfits(bias_file,bhead)
  mask=readfits(mask_file,mhead)
  mhead=modeinfo(mhead, inst_mode)
  mask=clipnflip(mask, mhead)
  mean_bias=median(bias)
  flat=readfits(norm_flat_file,fhead)
  npix=n_elements(flat[*,0])

  sp       =dblarr(npix,nord)
  cont     =sp+1.d0

;Image can be smaller than BIAS?
;Scatter light is affected by ThAr
;Optimal extraction has not enough S/N
; Use the whole order to make Valenti-style slit function
; Make low_signal_extraction routine as an alternative to slit_func

  for n=0,n_elements(spec)-1 do begin
    fits_read,spec[n],im,head,exten=exten  ; Read target spectrum NEW UVES format
;   im = readfits(spec[n], head)           ; Read target spectrum OLD UVES format
    head = modeinfo(head, inst_mode, XR=50+[0,2047]) ; Modify header with instrument setup
    im=clipnflip(im, head)                 ; Clip & flip
    im=im-bias                             ; Bias subtraction

;Extract frame information from the header
    readn = sxpar(head, 'E_READN')
    dark  = sxpar(head, 'E_BACKG')
    gain  = sxpar(head, 'E_GAIN')


; Subtract scattered light. The approximation is returned in 2D array
; bg for each inter-order troff.
;
;Parameter(s) to be tuned: swath_width - swath width in columns
;                          LAMBDA_SF   is the same as SF_SMOOTH (see hamflat)
;                          LAMBDA_SP   is the same as SP_SMOOTH (see hamflat)
;                          OSAMPLE   - pixel oversampling (see hamflat)
;
;    orders[0,*]=orders[0,*]-35
    mkscatter, im, orders, bg, ybg, colrange=col_range, LAMBDA_SF=100. $
             , SWATH_WIDTH=600, OSAMPLE=10, MASK=mask, LAMBDA_SP=100. $
             , GAIN=gain, READN=readn, /SUBTRACT
;Flat fielding
    im=im/flat

; Optimally extract science spectrum
;Parameter values should be similar to those used in the calls to mkscatter
;and hamspec
;
    getxwd, im, orders, xwd, sxwd, colrange=col_range ;get extraction width
    xwd=replicate(30.,2,nord)

    hamspec, im, head, orders, xwd, sxwd, or_range(0), sp $
           , sig=sig, colrange=col_range, SF_SMOOTH=10., OSAMPLE=12 $
           , swath_width=400, MASK=mask, /PLOT, FILENAME=spec[n]
    im=0
    sxaddpar, head, 'OBASE', or_range[0], ' base order number';, before='COMMENT'

    erase
    !p.multi=[4,1,4]
    for iord=0,nord-1 do begin
      c0=col_range(0,iord)
      c1=col_range(1,iord)
      x=indgen(c1-c0+1)+c0
;      if iord/3 gt 0 and iord/3 mod 4 eq 0 then j=get_kbrd(1)

      cc=blzcoef[*,iord]+mean_bias
      s=sp[*,iord]/cc
      c=top(s,2.d5,eps=0.002)
      cc=cc*c & s=s/c
      c=top(s,2.d5,eps=0.002)
      cc=cc*c & s=s/c
      cont[x,iord]=(sp[x,iord]/s)>1.
;      cont[x,iord]=blzcoef[x,iord]>1.
;      cc=median(sp[x,iord]/cont[x,iord])
;      cont[x,iord]=cont[x,iord]*cc

      plot,x,sp[*,iord]/cont[*,iord],xs=1,ys=3,xr=[0,npix-1] $
          ,charsize=2,title='Order:'+strtrim(iord,2)
    endfor

    fileout=strmid(spec[n],0,strpos(spec[n],'.',/REVERSE_SEARCH))+'.sp.ech'
    fileout=prefix+fileout
    wdech,fileout,head,sp,cont=cont,sig=sig $
         ,/overwrite ;save spectrum to disk
    print,'Extracted spectrum was written to: "'+fileout+'"'

    !p.multi=0
  endfor

;Done
end

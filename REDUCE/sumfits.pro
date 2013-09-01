Pro SumFits, List, Instrument, Summ, Head, DEBUG=debug, ERR=err $
           , THRES=thres, WIN=hwin, MASK=mask, XR=xr, YR=yr $
           , EXTEN=extension
;Subroutine to correct cosmic rays blemishes, while adding otherwise
;  similar images.
; List - (input string array) list of the filenames
; Summ - (output array(# columns,# rows)) coadded, cosmic ray corrected image.
; Head - (optional output string vector(# cards x 80)) FITS header associated
;        with coadded image.
; Optional parameters:
; debug - if set, plot scatter diagrams for cosmic ray corrections
; err   - string, if set, it contains an error message
;
;08-Aug-91 JAV  Fixed bug in signed to unsigned I*2 logic.
;25-Jun-92 Fen  Space padding logic for file name fixed.
;05-Sep-92 JAV  Fixed sign error in application of BZero.
;29-Mar-93 JAV  Renamed sumfits.pro to avoid conflict w/ intrinsic coadd().
;       Use FITS standard BZERO and BSCALE transformation:
;          Image = (BSCALE * ScaledImage) + BZERO
;       Detect older, nonstandard tranformations and abort.
;01-Apr-94 JAV  Equate FRAMENO (hires) and OBSNUM (Lick). Also equate ELAPTIME
;        (hires) and EXPOSURE (Lick). Use addcard to add cards to the
;        FITS header.
;10-Apr-94 JAV  Make index of exposure card a scalar for use in conditional.
;05-Dec-94 CMJ  Equate EXPTIME (McDonald) with ELAPTIME (hires) and
;                EXPOSURE (Lick).
;27-Mar-97 CMJ  Equate ITIME (CSHELL) with EXPOSURE (Lick).
;27-Mar-97 CMJ,JAV  Added logic to force the header to a multiple of 2880.
;                     Done by checking length of HSList.
;30-Aug-2000 NP Removed all traces of the instrument-specific behavour. Read all
;               files in the list line-by-line and use pipes for the compressed files.
;               Included new error handling via string ERR parameter.
;16-Nov-2005 NP Replaced all readfits with fits_read and pass extension number in exten
;               keyword. This way sumfits can handle fiels with multiple images.
;10-Oct-2011 NP Incorporating separate treatment of amplifiers when the CCD is read
;               using more than one port.
  err = ''
  if N_Params() lt 3 then begin
    Message,/Info,'Syntax: SumFits, List, Instrument, Im [,Head [, /DEBUG ' $
                 +'[, ERR=err[, THRES=thres[, WIN=win]]]]]'
    err = '%SUMFITS: Insufficient number of parameters'
    return
  endif

;Default values for optional parameters
  if(not keyword_set(thresh)) then thresh = 3.5       ;bad pixel threshold (sigmas)
  if(not keyword_set(hwin)) then hwin = 50            ;half window for column avg
  if(not keyword_set(extension)) then exten=0 else exten=extension

;Verify sensibility of passed parameters.
  Lst = List[sort(List)]
  Lst = Lst[uniq(Lst)]
  nFil = N_Elements(Lst)             ;number of files in list

;===========================================================================
  if nFil lt 2 then begin            ;true: only one image
    fits_read,Lst[0],Summ,head,exten=exten,/PDU
    Head = modeinfo(Head, Instrument, TIME=Exposure, GAIN=gain $
                        , READN=rdnoise, ORIENT=orient, XR=xr, YR=yr)
    if n_elements(Summ) eq 1 then begin
      err = '%SUMFITS: Failed reading file '+Lst[0]
    endif
    return
;===========================================================================
  endif else if nFil eq 2 then begin ;true: only two images
    fits_read,Lst[0],  im, h, exten=exten,/PDU
    h = modeinfo(h, Instrument, TIME=Exp1)
    fits_read,Lst[1], Summ, head, exten=exten,/PDU
    Summ=Summ+im
    im=0
    h = modeinfo(head, Instrument, TIME=Exp2, GAIN=gain $
                     , READN=rdnoise, ORIENT=orient, XR=xr, YR=yr)
    TotalExpo = Exp1 + Exp2
    if n_elements(Summ) eq 1 then begin
      err = '%SUMFITS: Failed reading file '+Lst[0]+Lst[1]
    endif
;Add info to header.
    sxaddpar, head, 'BZERO', 0.0
    sxaddpar, head, 'BSCALE', 1.0
    sxaddpar, head, 'EXPTIME', TotalExpo
    sxaddpar, head, 'DARKTIME', TotalExpo
;Because we do not devide the signal by the number of files the
;read-out noise goes up by the square root of the number of files
    sxaddpar, head, 'RDNOISE', rdnoise * sqrt(nFil) $
            , ' noise in combined image, electrons'
    sxaddpar, head, 'NIMAGES', nFil $
            , ' number of images summed'
    sxaddpar, head, 'HISTORY', 'Images coadded by sumfits.pro on ' + systime()
    head = modeinfo(Head, Instrument, XR=xr, YR=yr)
    return
  endif

;===========================================================================
;Initialize header information lists (one entry per file).
  BZList = FltArr(nFil)             ;BZERO list (0.0=default)
  BSList = BZList + 1.0             ;BSCALE list (1.0=default)
  ExpoList = BZList                 ;exposure time list
  TotalExpo = 0.d0                  ;init total exposure time
  LUnit = LonArr(nFil)              ;Array of units

;Loop through files in list, grabbing and parsing FITS headers.
  FNLen = max(StrLen(List))                     ;length of longest filename
  print,'  FILE' + String(Replicate(32b,FNLen-3)) + 'OBS COLS ROWS  OBJECT  EXPOSURE'
  LUnit_pool=128

  FName = Lst[0]
  fits_read, FName, i, Head, EXTEN=exten $      ;Read main header
           , /HEADER_ONLY, /PDU
  Head = modeinfo(Head, Instrument, TIME=Exposure, GAIN=gain $
                      , READN=rdnoise, ORIENT=orient, XR=xr, YR=yr)

  nCol_old = sxpar(Head, 'NAXIS1')              ;# columns
  nRow_old = sxpar(Head, 'NAXIS2')              ;# rows
  
  n_ampl = sxpar(Head, 'E_AMPL', count=count)   ; Check if we deal with multiple amplifiers
  if(count eq 0) then n_ampl = 1
  if(n_ampl eq 1) then begin
    xlo = sxpar(Head, 'E_XLO', count=count)     ; section(s) of the detector to process
    if(count eq 1) then xlo=[xlo]               ; convert scalar to a 1-element array
    xhi = sxpar(Head, 'E_XHI', count=count)
    if(count eq 1) then xhi=[xhi]
    ylo = sxpar(Head, 'E_YLO', count=count)
    if(count eq 1) then ylo=[ylo]
    yhi = sxpar(Head, 'E_YHI', count=count)
    if(count eq 1) then yhi=[yhi]
  endif else begin
    xlo = sxpar(Head, 'E_XLO*', count=count)    ; section(s) of the detector to process
    if(count eq 1) then xlo=[xlo]               ; convert scalar to a 1-element array
    xhi = sxpar(Head, 'E_XHI*', count=count)
    if(count eq 1) then xhi=[xhi]
    ylo = sxpar(Head, 'E_YLO*', count=count)
    if(count eq 1) then ylo=[ylo]
    yhi = sxpar(Head, 'E_YHI*', count=count)
    if(count eq 1) then yhi=[yhi]
  endelse
  
  nFix = 0                                      ;init fixed pixel counter
  Summ = FltArr(nCol_old,nRow_old)              ;init R*4 output image array

  linear = sxpar(Head, 'E_LINEAR', count=count) ; Check if non-linearity correction
  if(count eq 0) then linear = 1 eq 1           ; is needed and prepare for it.
  if(not linear) then begin
    pref = sxpar(Head, 'E_PREFMO', count=count)
  endif

  for amplifier=1,n_ampl do begin               ;Outer loop through amplifiers (note: 1,2 ...)
    for iFil=0,nFil-1 do begin                  ;loop though files
      FName = Lst[iFil]                         ;take the next filename

      fits_read, FName, i, Head, EXTEN=exten $  ;Read main header
               , /HEADER_ONLY, /PDU

      Head = modeinfo(Head, Instrument, TIME=Exposure, GAIN=gain $
                          , READN=rdnoise, ORIENT=orient, XR=xr, YR=yr)

      if(n_elements(Head) eq 1) then begin
        err = '%SUMFITS: ' + Head + ' in file '+FName
        if iFil gt 0 then for jFil=0,iFil-1 do close, LUnit[jFil]
        return
      endif

      nCol = sxpar(Head, 'NAXIS1')                ;# columns
      nRow = sxpar(Head, 'NAXIS2')                ;# rows

      if(nRow ne nRow_old or nCol ne nCol_old) then begin
        err = '%SUMFITS: file ' + FName + ' has different dimensions'
        if iFil gt 0 then for jFil=0,iFil-1 do close, LUnit[jFil]
        return
      endif
      BitPix = sxpar(Head, 'BITPIX')              ;# rows
      Object = sxpar(Head, 'OBJECT')              ;get object ID
      BSList[iFil] = sxpar(Head, 'BSCALE')        ;pixel scale
      if BSList[iFil] eq 0 then BSList[iFil] = 1
      BZList[iFil] = sxpar(Head, 'BZERO')         ;pixel offset

      LUnit_pool=LUnit_pool-1
      if(LUnit_pool eq 101) then $                ;skip LUnits 100 because header
        LUnit_pool=LUnit_pool-2                   ;reading above needs one unit

      if(LUnit_pool eq 100 or $
         LUnit_pool eq 106) then $                ;skip LUnits 101, 106 because header
        LUnit_pool=LUnit_pool-1                   ;reading above needs one unit

      close,LUnit_pool                            ;Make sure this unit is available

      LUnit[iFil] = fxposit(FName,exten,/READONLY $ ;Open file and position at the start
                           ,LUnit=LUnit_pool)       ;of the correct extension
      if(LUnit[iFil] lt 0) then begin
        err = '%SUMFITS: file '+ FName +' not found'
        if iFil gt 0 then for jFil=0,iFil-1 do close, LUnit[jFil]
        return
      endif

      mrd_hread, LUnit[iFil], h, status        ;Skip header

      TotalExpo = TotalExpo + Exposure            ;accumulate exposure times
      Blanks = ''                                 ;null insert string
      if StrLen(FName) lt FNLen then begin        ;true: need to pad with spaces
        Blanks = String(Replicate(32b,FNLen - StrLen(FName)))
      endif
      print,format='(2X,2A,I4,I5,I5,2A,I6)', $
        FName,Blanks,iFil,nCol,nRow,'  ',StrTrim(Object), Exposure ;summarize header info
    endfor

    xleft   = xlo[amplifier-1]
    xright  = xhi[amplifier-1]
    ybottom = ylo[amplifier-1]
    ytop    = yhi[amplifier-1]

    if(n_ampl gt 1) then begin
      gain_amp = gain[amplifier-1]
      rdnoise_amp = rdnoise[amplifier-1]
    endif else begin
      gain_amp = gain
      rdnoise_amp = rdnoise
    endelse
  
    case BitPix of
       8:   IDL_type = 1                         ; Byte
      16:   IDL_type = 2                         ; UInteger*2
      32:   IDL_type = 3                         ; UInteger*4
     -32:   IDL_type = 4                         ; Real*4
     -64:   IDL_type = 5                         ; Real*8
    else:   begin
              err = '%SUMFITS: Illegal value of BITPIX (= ' +  $
                               strtrim(BitPix,2) + ') in FITS header'
              for iFil=0,nFil-1 do close, LUnit[iFil]
              return
            end
    endcase

    if orient eq 0 or orient eq 2 or $    ;No need to rotate the image, proceed row-by-row
       orient eq 5 or orient eq 7 then begin

      nCol_a = xright - xleft + 1
      nRow_a = ytop - ybottom + 1

      MBuff  = FltArr(nCol_a,nFil) ;init multi-file buffer
      prob   = FltArr(nCol_a,nFil) ;init probability function
      mFit   = FltArr(nCol_a,nFil) ;init fit to mbuff
      Buff   = FltArr(nCol_a)      ;init final image data buffer
      Block  = make_array(nCol, type = IDL_type)

      ON_IOERROR,ioerr

      if(ybottom gt 0) then begin                 ; skip ybottom rows
        for iFil=0, nFil-1 do begin               ; loop though files
          readu, LUnit[iFil], Block
        endfor
      endif

      for iRow = ybottom, ytop do begin          ;loop through rows
        if iRow mod 100 eq 0 and iRow gt 0 then begin ;report status
          Message,/Info,StrTrim(iRow,2) + ' rows processed - ' + $
          StrTrim(String(nFix),2) + ' pixels fixed so far.'
        endif
        for iFil=0, nFil-1 do begin              ;loop though files
          readu, LUnit[iFil], Block
          Blockf = float(Block) * BSList[iFil] + BZList[iFil]
          if(not linear) then begin           ; Linearity needs fixing
;print,'1) Calling nonlinear correction'
            Blockf = call_function('nonlinear_' + pref, Blockf, Head $
                                , AMPLIFIER=amplifier, GAIN=gain_amp)
          endif
          if(keyword_set(mask)) then begin
            mBuff[*,iFil] = (Blockf * mask[*,iRow])[xleft:xright]
          endif else begin
            mBuff[*,iFil] = Blockf[xleft:xright]
          endelse
        endfor

;Construct a probability function based on mbuff data.
        for iCol=hwin, nCol_a-hwin-1 do begin
;          filwt = total(mBuff[iCol-hwin:iCol+hwin,*], 1)
          filwt = median(mBuff[iCol-hwin:iCol+hwin,*], dimension=1)

          tot_filwt = total(filwt)                        ;norm for probability
          if(tot_filwt gt 0.0) then filwt = filwt / tot_filwt
          prob[icol,*] = filwt
        endfor

        for iCol=0,hwin-1 do begin
          prob[iCol,*] = 2.0 * prob[hwin,*] - prob[2*hwin-iCol,*]
        endfor
;        iCol=indgen(hwin)
;        prob[iCol,*] = 2.0 * prob[hwin,*] - prob[2*hwin-iCol,*]

        for iCol=nCol_a-hwin, nCol_a-1 do begin
          prob[iCol,*] = 2.0 * prob[nCol_a-hwin-1,*] $
                          - prob[2*(nCol_a-hwin-1)-iCol,*]
        endfor
;        iCol=nCol_a-hwin+indgen(hwin)
;        prob[iCol,*] = 2.0 * prob[nCol_a-hwin-1,*] $
;                     - prob[2*(nCol_a-hwin-1)-iCol,*]

;Loop through columns, fitting data and constructing mfit.
;Reject cosmic rays in amplitude fit by ignoring highest and lowest point.
        for iCol=0, nCol_a-1 do begin
;          rat = mBuff[iCol,*] / prob[iCol,*]
;          iSort = sort(rat)
;          amp = total(rat[iSort[1:nFil-2]]) / (nFil-2)
;          mFit[iCol,*] = amp * prob[iCol,*]
          pr = reform(prob[iCol, *])
          mFit[iCol,*] = 0.0
          iProb = where(pr gt 0.0, nProb)
          if(nProb gt 0) then begin
            rat = reform(mBuff[iCol, iProb]) / pr[iProb]
            amp = (total(rat) - min(rat) - max(rat)) / (nFil-2)
            mFit[iCol,iProb] = amp * pr[iProb]
          endif
        endfor

;Construct noise model.
        predsig = sqrt(rdnoise_amp^2 + abs(mFit/gain_amp))
  
;Identify outliers.
        iBad = where(mBuff-mFit gt thresh*predsig, nBad)

;Debug plot.
        if keyword_set(debug) and iRow mod 10 eq 0 then begin
          plot, mBuff-mFit, xsty=3, ysty=3, ps=3 $
            , ytit='Data - Model  (ADU)' $
            , tit= 'Row = ' + strtrim(iRow, 2) $
                 + ',  Threshold = ' + string(thresh, form='(f9.4)')
          colors
          oplot, thresh*predsig, co=4
          oplot, -thresh*predsig, co=4
          if nBad gt 0 then begin
            oplot, iBad, mBuff[ibad]-mFit[ibad], ps=2, syms=1.0, co=2
          endif
          print, form='(a,$)', 'Push space to continue...'
          junk = get_kbrd(1)
          print, ''
        endif

;Fix bad pixels, if any.
        if nBad gt 0 then begin
          mBuff[iBad] = mFit[iBad]
          nFix = nFix + nBad
        endif

;Construct the summed FITS.
        Buff = total(mBuff, 2); / nfil
;plot,Buff,xs=3,ys=3,title='Row '+strtrim(irow+ybottom,2)
        Summ[xleft:xright,iRow] = Buff   ;insert final data into image
      endfor                             ; end of the iRow loop 
    endif else if orient eq 1 or orient eq 3 or $  ;Orders are vertical, read 2*hwin rows first
                  orient eq 4 or orient eq 6 then begin

      nCol_a = xright - xleft + 1
      mRow   = 2 * hwin + 1             ;# of rows in the FIFO buffer
      mBuff  = fltarr(nCol_a, mRow, nFil) ;init multi-file FIFO buffer
      mmBuff = fltarr(nCol_a, hwin, nFil) ;init multi-file sub-buffer
      prob   = fltarr(nCol_a,nFil)        ;init probability function
      mFit   = fltarr(nCol_a,nFil)        ;init fit to mbuff
      Buff   = fltarr(nCol_a)             ;init final image data buffer
      Block  = make_array(nCol, type = IDL_type)

      if(ytop - ybottom + 1 lt mRow) then begin
        err = '%SUMFITS: the number of rows should be larger than 2 * win = '+strtrim(mRow,2)
        for iFil=0,nFil-1 do begin                ;loop thru open files
          close, LUnit[iFil]                      ;free logical unit
        endfor
        return
      endif

;      ON_IOERROR,ioerr

      if(ybottom gt 0) then begin                 ;Skip rows that do not belong to this amplifier
        for jRow = 0, ybottom-1 do begin
          for iFil=0, nFil-1 do begin
            readu, LUnit[iFil], Block
          endfor
        endfor
      endif

      for jRow = 0, 2*hwin-1 do begin             ;fill out the initial window
        iRow = ybottom + jRow
        for iFil=0, nFil-1 do begin               ;loop though files
          readu, LUnit[iFil], Block
          Blockf = float(Block) * BSList[iFil] + BZList[iFil]
          if(not linear) then begin            ; Linearity needs fixing
;print,'2) Calling nonlinear correction'
            Blockf = call_function('nonlinear_' + pref, Blockf, Head $
                                , AMPLIFIER=amplifier, GAIN=gain_amp)
          endif
          if(keyword_set(mask)) then begin
            mBuff[*,jRow,iFil] = (Blockf * mask[*,iRow])[xleft:xright]
          endif else begin
            mBuff[*,jRow,iFil] = Blockf[xleft:xright]
          endelse
        endfor
      endfor
      mmBuff = mBuff[*,0:hwin-1,*]                        ;save that for probability extrapolation
;Analyse this!
      jRow = 2*hwin-1                                     ;last read row in the FIFO buffer
      cRow = hwin - 1                                     ;current row pointer in the FIFO buffer
;time1=0.0
;time2=0.0
      for iRow = ybottom+hwin, ytop-hwin do begin         ;loop through rows

;time0=systime(1)

        if iRow mod 100 eq 0 and iRow gt 0 then begin     ;report status
          Message,/Info,StrTrim(iRow,2) + ' rows processed - ' $
          + StrTrim(String(nFix),2) + ' pixels fixed so far.'; $
;          + StrTrim(time1)+StrTrim(time2)
        endif
        jRow = (jRow + 1) mod mRow                        ;increment FIFO new row pointer
        cRow = (cRow + 1) mod mRow                        ;increment FIFO current row pointer
        for iFil=0, nFil-1 do begin                       ;loop though files
          readu, LUnit[iFil], Block
          Blockf = float(Block) * BSList[iFil] + BZList[iFil]
          if(not linear) then begin                       ; Linearity needs fixing
;print,'3) Calling nonlinear correction'
            Blockf = call_function('nonlinear_' + pref, Blockf, Head $
                                 , AMPLIFIER=amplifier, GAIN=gain_amp)
          endif
          if(keyword_set(mask)) then begin
            mBuff[*,jRow,iFil] = (Blockf * mask[*,iRow])[xleft:xright]
          endif else begin
            mBuff[*,jRow,iFil] = Blockf[xleft:xright]
          endelse
        endfor

;Construct a probability function based on mBuff data.
        filwt = reform(total(mBuff, 2))                    ;boxcar average for iRow
;        filwt = median(mBuff, dimension=2)
        tot_filwt = total(filwt, 2)                        ;norm for probability
        inorm = where(tot_filwt gt 0, nnorm)
        if(nnorm gt 0) then filwt[inorm,*] = filwt[inorm,*] $
                                           / (tot_filwt[inorm]#replicate(1,nFil))
;        for iCol = 0, nCol_a-1 do filwt(iCol,*) = filwt(iCol,*) / total(filwt(iCol,*))
        prob = filwt                                       ;probability estimate for iRow

;Preparation for handling rows below hwin and above ytop-hwin
        if iRow eq ybottom + hwin then prob1 = prob ;we will use it later for iRow's
                                                    ;between 0 and hwin-1

;time1=time1+systime(1)-time0
;time0=systime(1)

        for iCol = 0, nCol_a-1 do begin
          pr = reform(prob[iCol, *])
          mFit[iCol,*] = 0.0
          iProb = where(pr gt 0.0, nProb)
          if(nProb gt 0) then begin
            rat = reform(mBuff[iCol, cRow, iProb]) / pr[iProb]
;            iSort = sort(rat)
;            amp = total(rat[iSort[1:nFil-2]]) / (nFil-2)
            amp = (total(rat) - min(rat) - max(rat)) / (nFil-2)
            mFit[iCol,iProb] = amp * pr[iProb]
          endif
        endfor

;time2=time2+systime(1)-time0
;stop

;Construct noise model.
        predsig = sqrt(rdnoise_amp^2 + (mFit/gain_amp))

;Identify outliers.
        Buff = reform(mBuff[*,cRow,*])
        iBad = where(Buff - mFit gt thresh*predsig, nBad)

;Debug plot.
        if keyword_set(debug) and iRow mod 10 eq 0 then begin
          plot, Buff-mFit, xsty=3, ysty=3, ps=3 $
              , ytit='Data - Model  (ADU)' $
              , tit= 'Row = ' + strtrim(iRow, 2) $
                   + ',  Threshold = ' + string(thresh, form='(f9.4)')
          colors
          oplot, thresh*predsig, co=4
          oplot, -thresh*predsig, co=4
          if nBad gt 0 then begin
            oplot, iBad, Buff[ibad]-mFit[ibad], ps=2, syms=1.0, co=2
          endif
          print, form='(a,$)', 'Push space to continue...'
          junk = get_kbrd(1)
          print, ''
        endif

;Fix bad pixels, if any.
        if nBad gt 0 then begin
          Buff[iBad] = mFit[iBad]
          nFix = nFix + nBad
        endif

;Construct the summed flat.
        Buff = total(Buff, 2); / nfil
;plot,Buff,xs=3,ys=3,title='Row '+strtrim(irow,2)
        Summ[xleft:xright,iRow] = Buff                 ;insert final data into image

;1st special case: rows less than hwin from the 0th row
        if iRow gt ybottom+hwin and iRow le ybottom+2*hwin then begin ;extrapolate probability
                                                                      ;from hwin < iRow < 2*hwin+1
          lRow = 2*hwin - iRow + ybottom                              ;in the FIFO buffer
          for iCol = 0, nCol_a-1 do begin
            prob2 = reform(2.0 * prob1[iCol, *] - prob[iCol, *])
;            rat = reform(mmBuff[iCol, lRow, *]) / prob2
;            iSort = sort(rat)
;            amp = total(rat[iSort[1:nFil-2]]) / (nFil-2)
;            mFit[iCol,*] = amp * prob2
            mFit[iCol,*] = 0.0
            iProb = where(prob2 gt 0.0, nProb)
            if(nProb gt 0) then begin
              rat = reform(mmBuff[iCol, 2*hwin-cRow, iProb]) / prob2[iProb]
              amp = (total(rat) - min(rat) - max(rat)) / (nFil-2)
              mFit[iCol,iProb] = amp * prob2[iProb]
            endif
          endfor
;Construct noise model.
          predsig = sqrt(rdnoise_amp^2 + (mFit/gain_amp))

;Identify outliers.
          Buff = reform(mmBuff[*,lRow,*])
          iBad = where(Buff - mFit gt thresh*predsig, nBad)

;Fix bad pixels, if any.
          if nBad gt 0 then begin
            Buff[iBad] = mFit[iBad]
            nFix = nFix + nBad
          endif

;Construct the summed flat.
          Buff = total(Buff, 2); / nfil
          Summ[xleft:xright,ybottom + hwin - (iRow - ybottom - hwin)] = Buff ;insert final data into image
        endif
;Save probabilities for ytop-2*hwin < iRow < ytop-hwin
        if iRow ge ytop - 2*hwin and iRow lt ytop - hwin then begin
          mmBuff[*,iRow - (ytop - 2*hwin), *] = prob
        endif
      endfor

;2nd special case: rows greater than ytop-hwin
      prob1 = prob                             ;save probability from nRow - hwin
      for iRow = ytop - hwin + 1, ytop do begin
        prob = reform(mmBuff[*,ytop - iRow,*])
        cRow = (cRow + 1) mod mRow             ;the corresponding row location in the FIFO buffer
        for iCol = 0, nCol_a-1 do begin
          prob2 = reform(2.0 * prob1[iCol, *] - prob[iCol, *])
;          rat = reform(mBuff[iCol, cRow, *]) / prob2
;          iSort = sort(rat)
;          amp = total(rat[iSort[1:nFil-2]]) / (nFil-2)
;          mFit[iCol,*] = amp * prob2
          mFit[iCol,*] = 0.0
          iProb = where(prob2 gt 0.0, nProb)
          if(nProb gt 0) then begin
            rat = reform(mBuff[iCol, cRow, iProb]) / prob2[iProb]
            amp = (total(rat) - min(rat) - max(rat)) / (nFil-2)
            mFit[iCol,iProb] = amp * prob2[iProb]
          endif
        endfor

;Construct noise model.
        predsig = sqrt(rdnoise_amp^2 + (mFit/gain_amp))

;Identify outliers.
        Buff = reform(mBuff[*,cRow,*])
        iBad = where(Buff - mFit gt thresh*predsig, nBad)

;Fix bad pixels, if any.
        if nBad gt 0 then begin
          Buff[iBad] = mFit[iBad]
          nFix = nFix + nBad
        endif

;Construct the summed flat.
        Buff = total(Buff, 2); / nfil
        Summ[xleft:xright,iRow] = Buff  ;insert final data into image
;print,2*(ytop-hwin)-irow,irow
;plot,Summ[xleft:xright,2*(ytop-hwin)-irow],xs=1,yr=[2.5d4,3.d4]&oplot,Summ[xleft:xright,iRow],col=c24(2)
;stop
      endfor
    endif
    For iFil=0,nFil-1 Do Begin        ;loop thru files
      close, LUnit[iFil]       ;free logical unit
    EndFor
  endfor

  ON_IOERROR,NULL

  Message, /Info $
    , 'Total cosmic ray hits identified and removed: ' $
    + StrTrim(String(nFix), 2)

;Add info to header.
  sxaddpar, head, 'BZERO', 0.0
  sxaddpar, head, 'BSCALE', 1.0
  sxaddpar, head, 'EXPTIME', TotalExpo
  sxaddpar, head, 'DARKTIME', TotalExpo
;Because we do not devide the signal by the number of files the
;read-out noise goes up by the square root of the number of files
;  sxaddpar, head, 'RDNOISE', rdnoise / sqrt(nFil) $
;          , ' noise in combined image, electrons'
  if(n_elements(rdnoise) gt 1) then begin
    for amplifier=1,n_elements(rdnoise) do begin
      sxaddpar, head, 'RDNOISE'+suffix(amplifier,1) $
              , rdnoise[amplifier-1] * sqrt(nFil) $
              , ' noise in combined image, electrons'
    endfor
  endif else begin
    sxaddpar, head, 'RDNOISE', rdnoise * sqrt(nFil) $
            , ' noise in combined image, electrons'
  endelse
  sxaddpar, head, 'NIMAGES', nFil $
          , ' number of images summed'
  sxaddpar, head, 'NPIXFIX', nFix $
          , ' pixels corrected for cosmic rays'
  sxaddpar, head, 'HISTORY', 'Images coadded by sumfits.pro on ' + systime()

  if(not linear) then begin ; Non-linearity was fixed. Mark this in the header
    i = where(strpos(head, 'E_LINEAR') ge 0, ni)
    if(ni gt 0) then head = [head[0:i-1], head[i+1:*]]
    sxaddpar, head, 'E_LINEAR', 'T'  $
            , ' Image corrected of non-linearity'

    ii = where(strpos(head, 'E_GAIN*') ge 0, ni)
    if(ni gt 0) then begin
      for i=0,ni-1 do begin
        k=ii[i]
        head = [head[0:k-1], head[k+1:*]]
      endfor
    endif          
    sxaddpar, head, 'E_GAIN', 1  $
            , ' Image was converted to e-'
  endif

  Head = modeinfo(head, Instrument, XR=xr, YR=yr)

  return

ioerr:
  err = '%SUMFITS: Input error reading file ' + Lst[iFil]
  For iFil=0,nFil-1 Do Begin        ;loop thru files
    close, LUnit[iFil]       ;free logical unit
  EndFor
  ON_IOERROR,NULL
  return

End


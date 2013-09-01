function modeinfo, header, specmode, xr=xr, yr=yr, orient=orient $
  , gain=gain, readn=readn, backg=backg, time=time $
  , ra2000=ra2000, de2000=de2000, jd=jd $
  , obslon=obslon, obslat=obslat, obsalt=obsalt $
  , helcorr=helcorr
;+
;Modify FITS header to include a uniform description of spectrograph
; characteristics relevant to spectral extraction. Characteristics are
; either passed by keyword or extracted from the FITS header using
; instrument specific parsing code.
;
;Inputs:
; [specmode] (string) Standard name of spectrograph mode. If not specified,
;   routine will attempt to deduce the spectrograph mode from information
;   in the FITS header.
; [orient=] (integer) values of orient (E_ORIENT) to put in header.
; [xr=] (vector(2)) value of trim range [E_XLO,E_XHI] to put in header.
; [yr=] (vector(2)) values of trim range [E_YLO,E_YHI] to put in header.
; [gain=] (scalar) value of gain (E_GAIN) to put in header.
; [readn=] (scalar) value of read noise (E_READN) to put in header.
; [backg=] (scalar) value of background (E_BACKG) to put in header.
; [time=] (scalar) value of exposure time (E_TIME) to put in header.
; [ra2000=] (scalar) value of right ascension of object (E_RA2000) in hrs to put in header.
; [de2000=] (scalar) value of declination of object (E_DE2000) in deg to put in header.
; [jd=] (scalar) value of Julian date - 2400000 at middle of exposure (E_JD) to put in header.
; [obslon=] (scalar) value of longitude of observatory (E_OBSLON) in deg to put in header.
; [obslat=] (scalar) value of latitude of observatory (E_OBSLAT) in deg to put in header.
; [obsalt=] (scalar) value of altitude of observatory (E_OBSALT) in m to put in header.
; [helcorr=] (scalar) value of heliocentric correction (E_HELCOR) in km/s to put in header.
;
;Input/Output:
; header (string array) FITS header associated with image. Input header will
;   be modified by modifying or adding the following header cards:
;     E_ORIENT (integer)  used to reorient: image=rotate(image, E_ORIENT)
;     E_XLO    (integer)  used to trim: image=image(E_XLO:E_XHI,E_YLO:E_YHI)
;     E_YLO    (integer)    [see previous]
;     E_XHI    (integer)    [see previous]
;     E_YHI    (integer)    [see previous]
;     E_GAIN   (float)    conversion factor (electrons/ADU)
;     E_READN  (float)    RMS read noise (electrons)
;     E_BACKG  (float)    mean background counts (electrons/pixel)
;     E_TIME   (float)    exposure time (seconds)
;     E_RA2000 (float)    right ascension of object for epoch 2000.0 (hrs)
;     E_DE2000 (float)    declination of object for epoch 2000.0 (deg)
;     E_JD     (float)    Julian date - 2400000 at middle of exposure
;     E_OBSLON (float)    longitude of observatory (deg)
;     E_OBSLAT (float)    latitude of observatory (deg)
;     E_OBSALT (float)    altitude of observatory (deg)
;     E_HELCOR (float)    heliocentric correction (km/s)
;     E_HJD    (float)    heliocentric Julian date -2400000 at middle of exposure
;
;Output:
; Return value of function is either the modified FITS header image (string
;   vector, E_* cards added) or an error message (scalar string).
;
;Checklists for Programmers:
;  Adding a new spectrograph mode:
;    1) Add new prefix and name in file "known_modes.txt".
;    2) Create additional procedure file "modeinfo_prefix.pro"
;    3) Set instrument parameters in this file, as needed.
;
;  Adding new header parameter for all spectrograph modes:
;    1) Add keyword in procedure definition at top of routine.
;    2) Add description of change to History section in header.
;    3) Update syntax print statement in syntax check section.
;    4) Add description and default value in defaults section.
;    5) Add keyword override logic in keyword handling section.
;    6) Add element to inst_setup structure in build section.
;    7) Add code in other package routines to interpret new element.

;Notes on E_ORIENT:
;
;Reorientation occurs immediately after image data is read from disk.
;Images must be transformed so that echelle dispersion is mainly along
;the X-axis, and cross-dispersion is mainly along the Y-axis. After
;reorientation, wavelength should increase with both X-index (from left
;to right) and Y-index (bottom to top, assuming !order=0). This implies
;that red orders will be at the top of an image displayed with !order=0.
;The reorientation flag value is passed as the second argument to the
;IDL function rotate(). The codes required to correct various initial
;orientations are given below:
;
; Use             [original frame orientation]
;orient
; 0:    red orders at top,    orders run left to right
; 1:    red orders at right,  orders run top to bottom
; 2:    red orders at bottom, orders run right to left
; 3:    red orders at left,   orders run bottom to top
; 4:    red orders at right,  orders run bottom to top
; 5:    red orders at top,    orders run right to left
; 6:    red orders at left,   orders run top to bottom
; 7:    red orders at bottom, orders run left to right
;
;Notes on E_XLO, E_YLO, E_XHI, E_YHI:
;
;Coordinates (X,Y) of the botton left (BL) and top right (TR) corners of the
;subimage that is to be reduced, i.e. the portion on the image that is not
;overscan, physically masked, subject to amplifier ringing, etc. Note that
;indices for both axes begin at 0, not 1. Also, the coordinates must be
;specified in the coordinate system BEFORE APPLYING THE RE-ORIENTATION
;DESCRIBED ABOVE.
;
;Notes on E_RA2000, E_DE2000, E_JD, E_OBSLON, E_OBSLAT, E_OBSALT, E_HJD
;and E_HELCOR.
;
;Coordinates of object and observatory as well as the point in time of the
;observation are required for calculating the heliocentric correction (helcorr.pro).
;The correction is already calculated here (see OHP section) but all the
;keywords are added to show the values used for the calculation. This is only
;done for files of scientific frames (stellar spectra). The correction is
;applied (i.e the keyword BARYCORR modified or added) in wavecal.pro only
;for the files already having the E_HELCOR keyword in their header section.
;
;History:
; 1999-Nov-27 Valenti      Wrote.
; 2000-Aug-24 Valenti      Adapted from inst_setup.pro.
; 2001-Mar-16 Piskunov     Changed the sequence of Clip and Flip and added
;                          the relation between ORIENT and the original frame
;                          layout.
; 2003-Nov-08 Mittermayer  Added E_RA2000, E_DE2000, E_JD, E_OBSLON, E_OBSLAT,
;                          E_OBSALT, E_HJD, E_HELCOR for heliocentric correction.
; 2004-Feb-19 Kochukhov    Added setup for SARG, corrected format for E_JD and E_HJD
; 2004-May    Kudryavtsev  Fixed wrong calculation of JDs (for NES) in some
;                          specific cases
; 2005-Apr-26 Kochukhov    Fixed wrong extraction of Observatory coordinates and JD
;                          for UVES observations
; 2005-Dec-19 Piskunov     Added HET HRS spectrograph description for the red CCD 
; 2006-Dec-11 Piskunov     Added Gemini-South bHROS spectrometer (single CCD in the extention)
; 2008-Jan-16 Piskunov     Added Keck HiRes spectrograph description for the blue, middle and red CCDs 
; 2008-May-07 Kochukhov    Known modes and procedure prefixes are read from file
; 2011-Oct-07 Piskunov     Added the number of amplifiers and arrays of
;                          amplifier-dependent parameters when detectors are
;                          read through more than one port
; 2013-Jan-30 Piskunov     Added prototype description for the NIR arm of the ESO X-Shooter
;-

  common known_modes,prefixes,valid_modes

  if n_params() lt 1 then begin
    print, 'syntax: newhead = modeinfo(header [,specmode ,orient= ,gain=]'
    print, '          [,readn= ,backg= ,time=])'
    print, '          [,ra2000=, de2000=, obslon=, obslat=, obsalt=])'
    print, '          [,jd=, helcorr=])'
    return, 'Syntax error'
  endif

;Check that FITS header has proper structure.
  sz = size(header)
  if sz(0) ne 1 or sz(2) ne 7 then begin
    message, /info, 'Header must be a string vector'
    return, 'Unable to put spectrograph mode information into header'
  endif

  newhead = header

;Make sure specmode is defined.
  if not keyword_set(specmode) then begin
    return, 'Unable to automatically detect spectrograph mode'
  endif

;Read list of valid modes and corresponding procedure prefixes from file 'modes.txt'
  if not keyword_set(valid_modes) then begin
    help,calls=a
    a=strmid(a(0),strpos(a(0),'<')+1,strlen(a(0))-strpos(a(0),'<')-1)
    if(strpos(a,'/nfs') eq 0) then a=strmid(a,4)
    if(strpos(a,'/export') eq 0) then a=strmid(a,7)
    if(!VERSION.OS eq 'Win32') then delimiter='\' else delimiter='/'
    path=strmid(a,0,strpos(a,delimiter,/REVERSE_SEARCH))

    str = ' ' 
    openr, munit, path+delimiter+'known_modes.txt', /get_lun
    while not eof(munit) do begin
      readf, munit, str

;Not am empty line, not a comment ; #
      if strlen(strtrim(str,2)) gt 1 and strmid(str,0,1) ne ';' and strmid(str,0,1) ne '#' then begin

;Procedure prefix name starts from 1st position, valid mode name starts after space
        if strmid(str,0,1) ne ' ' then prefix = strtrim(str,2) else begin
          prefixes = keyword_set(prefixes) ? [prefixes,prefix] : prefix
          valid_modes = keyword_set(valid_modes) ? [valid_modes,strtrim(str,2)] : strtrim(str,2)
        endelse
      endif
    endwhile
    close, munit
  endif
  
;Check that spectrograph mode is valid.
  mode = strlowcase(strtrim(specmode, 2))       ;canonical form
  imode = where(valid_modes eq mode, nmode)     ;check if known mode
  if nmode ne 1 then begin
    message, /info, 'Unknown spectrograph mode: ' + mode
    maxwid = 65                     ;width of display
    outstr = 'Valid spectrograph modes: '
    for i=0, n_elements(valid_modes)-1 do begin
      if strlen(outstr + valid_modes(i)) gt maxwid-2 then begin
        message, /info, outstr
        outstr = '... '
      endif
      outstr = outstr + valid_modes(i) + ' '
    endfor
    message, /info, outstr
    return, 'Unable to put spectrograph mode information into header'
  endif

;Call routine defined by in prefixes array
  call_procedure,'modeinfo_' + prefixes[imode[0]],mode,newhead,object,orient,xlo,xhi,ylo,yhi,gain $
                ,readn,backg,time,obslat,obslon,obsalt,ra2000,de2000,jd,NONLINEAR=nonlinear

;Put spectrograph mode information into output header.
  if n_elements(obslon) gt 0 then begin             ; Add heliocentric correction
    helcorr,obslon,obslat,obsalt,ra2000,de2000,jd $ ; calculating heliocentric
           ,corr,hjd                                ; correction and HJD
    sxaddpar, newhead, 'E_OBSLON'  , obslon     $
            , ' Longitude of observatory in degrees'
    sxaddpar, newhead, 'E_OBSLAT'  , obslat     $
            , ' Latitude of observatory in degrees'
    sxaddpar, newhead, 'E_OBSALT'  , obsalt     $
            , ' Altitude of observatory in meters'
    sxaddpar, newhead, 'E_RA2000'  , ra2000     $
            , ' Right ascension of object for 2000 in hours'
    sxaddpar, newhead, 'E_DE2000'  , de2000     $
            , ' Declination of object for 2000 in degrees'
    sxaddpar, newhead, 'E_JD'  , jd     $
            , ' Julian date - 2400000 at middle of exposure', form = '(F13.6)'
    sxaddpar, newhead, 'E_HJD'  , hjd     $
            , ' Heliocentric Julian date - 2400000 at middle of exposure', form = '(F13.6)'
    sxaddpar, newhead, 'E_HELCOR'  , corr     $
            , ' Baricentric correction in km/s'
  endif
  sxaddpar, newhead, 'E_HVERS', 1.1       $
          , ' Header keyword version number'
  sxaddpar, newhead, 'E_SPMODE', specmode $
          , ' Instrument mode'
  sxaddpar, newhead, 'E_PREFMO', prefixes[imode[0]] $
          , ' Instrument mode'
  sxaddpar, newhead, 'E_ORIENT', orient   $
          , ' Reorientation flag for extraction'
  n1 = n_elements(xlo)
  n2 = n_elements(ylo)
  n3 = n_elements(xhi)
  n4 = n_elements(yhi)
  n5 = n_elements(gain)
  n6 = n_elements(readn)
  n7 = n_elements(backg)
  if(n1 ne n2 or n1 ne n3 or n1 ne n4 or $
     n1 ne n5 or n1 ne n6 or n1 ne n7) then begin
    message,'Inconsistent dimensions of CCD parameters'
  endif
  
  if(keyword_set(nonlinear)) then begin
    linear = sxpar(Newhead, 'E_LINEAR', count=count)
    if(count eq 0) then sxaddpar, newhead, 'E_LINEAR' , 'F'  $
                                , ' Image needs linearity correction'
  endif

  sxaddpar, newhead, 'E_AMPL' , n1  $
              , ' Number of amplifiers used for readout'
  if(n1 gt 1) then begin ; Multiple amplifiers
    for amp=0,n1-1 do begin
      sf = suffix(amp+1, 1)
      sxaddpar, newhead, 'E_XLO'+sf   , xlo[amp]  $
              , ' Lowest X index to extract (0 base)'
      sxaddpar, newhead, 'E_YLO'+sf   , ylo[amp]  $
              , ' Lowest Y index to extract (0 base)'
      sxaddpar, newhead, 'E_XHI'+sf   , xhi[amp]  $
              , ' Highest X index to extract (0 base)'
      sxaddpar, newhead, 'E_YHI'+sf   , yhi[amp]  $
              , ' Highest Y index to extract (0 base)'
      sxaddpar, newhead, 'E_GAIN'+sf  , gain[amp] $
              , ' Gain (electrons/ADU)'
      sxaddpar, newhead, 'E_READN'+sf , readn[amp]$
              , ' Read noise (electrons)'
      sxaddpar, newhead, 'E_BACKG'+sf , backg[amp]$
              , ' Total background dark and sky (electrons)'
    endfor
  endif else begin
    sxaddpar, newhead, 'E_XLO'   , xlo      $
            , ' Lowest X index to extract (0 base)'
    sxaddpar, newhead, 'E_YLO'   , ylo      $
            , ' Lowest Y index to extract (0 base)'
    sxaddpar, newhead, 'E_XHI'   , xhi      $
            , ' Highest X index to extract (0 base)'
    sxaddpar, newhead, 'E_YHI'   , yhi      $
            , ' Highest Y index to extract (0 base)'
    sxaddpar, newhead, 'E_GAIN'  , gain     $
            , ' Gain (electrons/ADU)'
    sxaddpar, newhead, 'E_READN' , readn    $
            , ' Read noise (electrons)'
    sxaddpar, newhead, 'E_BACKG' , backg    $
            , ' Total background dark and sky (electrons)'
  endelse
  sxaddpar, newhead, 'E_TIME'  , time     $
          , ' Exposure time (seconds)'

  return, newhead

end

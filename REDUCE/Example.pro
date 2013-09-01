!error = 0				;clear any existing errors
@precompile				;precompile routines for nicer logs

;------------------------- Make changes below line -------------------------

;Specify filenames here. Do not include the (mandatory) ".fits" or ".dsk"
;  extensions. Comment out definitions of narr,wide,thar, and/or spec to
;  suppress processing of these images.
  set = 'df-1'				;setting name
  prefix = 'df'				;prefix for narr,wide,thar,spec
    narr = '3'				;  default order locations
    wide = '4'				;  wide flat
    thar = '5'				;  thorium-argon
    spec = ['3','6']			;  target spectra
  suffix = ''				;suffix for narr,wide,that,spec

;Change global variables (trace,image,dord,id,xwid,bg,bin,spike,roff,coff) here.
  hamset,id=8

;------------------------- Make changes above line -------------------------

;REMEMBER: This is a batch file, not a procedure, so begin/end blocks
;	     are not permitted....

;Print filename specifications for log file.
  print,'% set = "'+ set + '"'
  print,'% prefix = "' + prefix + '"'
  if keyword_set(narr) $
    then print,'%   narr = "' + narr + '"' $
    else print,'% No default order location image specified.'
  if keyword_set(wide) $
    then print,'%   wide = "' + wide + '"' $
    else print,'% No wide flat image specified.'
  if keyword_set(thar) $
    then print,'%   thar = "' + thar +'"' $
    else print,'% No Thorium-Argon Lamp image specified.'
  if keyword_set(spec) $
    then print,'%   spec =' , '"' + spec + '"' $
    else print,'% No target spectra specified.'
  print,'% suffix = "' + suffix + '"'

;Prepare to reduce spectra. Errors here are fatal and kill the batch job (or
;  exit IDL, if this file was invoked interactively with the "@" command).
  if !error ne 0 then exit			;check for errors
  if keyword_set(narr) $
    then hamdord,set,prefix+narr+suffix		;find default orders
  if !error ne 0 then exit			;check for errors
  if keyword_set(wide) $
    then hamflat,set,prefix+wide+suffix		;normalize wide flat
  if !error ne 0 then exit			;check for errors

;Reduce individual spectra, including the Th-Ar lamp spectrum.
  if keyword_set(thar) $
    then hamspec,set,prefix+thar+suffix,/thar	;reduce Th-Ar lamp spectrum
  nspec = n_elements(spec)			;number of spectra to reduce
  if nspec gt 0 $
    then for i=0,nspec-1 $
      do hamspec,set,prefix+spec(i)+suffix

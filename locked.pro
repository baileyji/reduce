function locked,file,DURATION=duration
; Function "locked" makes a temporary file with specific to
; allow concurrent processing of many files by
; several parallel IDL programs
  if(n_params() lt 1) then begin
    print,'Usage: if(not locked(<file_name>[,DURATION=lock_lifetime])) then ...'
    print,'   Function "locked" creates a temporary binary file which has'
    print,'   one record:'
    print,'"LOCK",x,"LOCK",y,"LOCK"'
    print,'   where "LOCK" is a 4-character string, x is a 32-bit'
    print,'   unsigned integer giving the creation time of the lock'
    print,'   in seconds and y is a 32-bit unsigned integer giving the'
    print,'   lock lifetime in seconds.'
    print,'   The lock lifetime is given by an optional parameter. The'
    print,'   default value is 3600 seconds. If the duration time is'
    print,'   set explicitely it supersedes the value in the lock file.'
    print,'A call to lock returns "false" if the lock file does not'
    print,'exist or already expired. In this case a new lock file is'
    print,'created. If the file exists but it has structure different'
    print,'from a lock file or if the lock did not expire the function'
    print,'returns "true".'
    return,0B
  endif
  
  if(not keyword_set(duration)) then begin
    yy=3600UL ; default duration time
    yy_set=0B
  endif else begin
    yy=ulong(duration)>1UL
    yy_set=1B
  endelse

  if(file_test(file)) then begin ; File exists
    openr,un,file,/GET_LUN,/SWAP_IF_BIG_ENDIAN
    s1=bytarr(4) & s2=s2 & s3=s1
    x=0UL & y=0UL
    readu,un,s1,x,s2,y,s3
    free_lun,un
    if(string(s1) eq 'LOCK' and string(s2) eq 'LOCK' and $
       string(s3) eq 'LOCK') then begin  ; This is a lock file
      curr_time=ulong(systime(/SECONDS))
      if(yy_set) then y=yy
      if(x+yy gt curr_time) then begin   ; We have an active lock
        return,1B
      endif else begin                   ; We have an expired lock
        openw,un,file,/GET_LUN,/SWAP_IF_BIG_ENDIAN
        writeu,un,'LOCK',ulong(systime(/SECONDS)),'LOCK',yy,'LOCK'
        free_lun,un
        return,0B
      endelse
    endif else return,1B                 ; This is not a lock file
  endif else begin                       ; Lock file does not exist, make one
    openw,un,file,/GET_LUN,/SWAP_IF_BIG_ENDIAN
    writeu,un,'LOCK',ulong(systime(/SECONDS)),'LOCK',yy,'LOCK'
    free_lun,un
    return,0B
  endelse
end
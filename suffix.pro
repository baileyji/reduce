Function suffix,v,suf_len,SYMBOL=suf_char
  if(n_params() lt 1) then begin
    print,' Suffix forms a suffix of a given length out of a number or an'
    print,' array of strings if the input is an array of numbers.'
    print,' Syntax: s=suffix(value,len[,SYMBOL=''_''])'
    print,'         where s is a string and len is a positive integer.'
    print,' If the length of value converted to string is larger then len'
    print,' it is trimmed on the left otherwise it is padded with zeroes on'
    print,' the left. If len is absent, the string containing value is returned.'
    print,' If SYMBOL is set it is converted to a string and the first'
    print,' symbol is used for padding.'
    return,''
  endif
  
  s=strtrim(v,2)
  if(n_params() eq 1) then return,s
  if(suf_len le 0) then return,s

  if(keyword_set(suf_char)) then symb=strmid(string(suf_char),0,1) else symb='0'
  n=n_elements(s)
  if(n eq 1) then begin ; v is a scalar
    if(strlen(s) eq suf_len) then return,s
    if(strlen(s) gt suf_len) then begin
      return,strmid(s,strlen(s)-suf_len,suf_len)
    endif else if(strlen(s) lt suf_len) then begin
      s=string(replicate((byte(symb))[0],suf_len-strlen(s)))+s
      return,s
    endif
  endif else begin                  ; v is an array
    for i=0l,n-1L do begin
      if(strlen(s[i]) gt suf_len) then begin
        s[i]=strmid(s[i],strlen(s[i])-suf_len,suf_len)
      endif else if(strlen(s[i]) lt suf_len) then begin
        s[i]=string(replicate((byte(symb))[0],suf_len-strlen(s[i])))+s[i]
      endif
    endfor
    return,s
  endelse
end

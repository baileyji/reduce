function line_label_position,wl,wlrange,FLEX=flexibility,index=l $
                            ,nindex=nline,lbl_shift=shift

  if(not keyword_set(flexibility)) then flex=30. else flex=1./(flexibility>0.001)
  w1=double(min(wlrange))
  w2=double(max(wlrange))
  l=where(wl ge w1 and wl le w2, nline)
  if(nline eq 0) then return,w1-1        ;No lines found in the range  
  if(nline eq 1) then return,wl[l]       ;For a single line return the actual
                                         ;location

  dx=(w2-w1)/nline
  shift=0.05*dx
  k=dindgen(nline)
  pos=((w1+dx*(k+0.5))*nline/flex+wl(l))/(nline/flex+1.)

  return,pos
end

pro dsp_defaults,font,nfont,default_device,box_clr
  case !version.os of
   'linux': begin
;             font ='-adobe-courier-medium-o-normal--14-100-100-100-m-90-iso8859-1'
;             nfont='-adobe-courier-bold-r-normal--14-100-100-100-m-90-iso8859-1'
              DEFAULT_DEVICE='X'
              box_clr=255
            end
   'Win32': begin
              font='arialb*bold*12'
              nfont='courier*10'
              DEFAULT_DEVICE='WIN'
              box_clr=255
            end
   'hp-ux': begin
              font='-adobe-times-medium-r-normal--14-100-100-100-p-74-iso8859-1'
              nfont='-adobe-courier-medium-r-normal--0-120-75-75-m-0-iso8859-1'
              DEFAULT_DEVICE='X'
              box_clr=255
            end
   'darwin': begin
              font='-adobe-times-medium-r-normal--14-100-100-100-p-74-iso8859-1'
              nfont='-adobe-courier-medium-r-normal--0-120-75-75-m-0-iso8859-1'
              DEFAULT_DEVICE='X'
              box_clr=255
            end
     'AIX': begin
              font='ncenr18'
              nfont='courier*bold*10'
              DEFAULT_DEVICE='X'
              box_clr=0
            end
  'ultrix': begin
              font='newcenturyschlbk_roman18'
              nfont='courier'
              DEFAULT_DEVICE='X'
              box_clr=0
            end
   'sunos': begin
              font='-adobe-courier-bold-r-normal--14-140-75-75-m-90-iso8859-1'
              nfont='courier'
              DEFAULT_DEVICE='X'
              box_clr=0
            end
     else : begin
              font='-*-courier-bold-r-*-*-12-*-*-*-*-*-*-*'
              nfont='courier'
              DEFAULT_DEVICE=!D.NAME
              box_clr=255
            end
    endcase
  return
end

pro draw,cs,csscale,csscale_thar,line,nline,title,print=print, $
    plot=default_device,reset=reset,thar_spec=thar_spec, $
    thar_flag=thar_flag,thar_list=thar_list,lambda_scale=wl

  i=where(cs gt 0, ni)
  csmax=max(cs(i)) & csmin=min(cs(i))
  if(not keyword_set(csscale) or keyword_set(reset)) then begin
    csscale=[0.2*(csmax-csmin)/(median(cs(i))-csmin), $
             csmin*0.2*(csmax-csmin)/(median(cs(i))-csmin)-csmax]
  endif
  if(keyword_set(print)) then begin
    set_plot,'ps'
    psfile=strmid(title,21,strlen(title)-21)
    if(strpos(psfile,'.') ge 0) then psfile=strmid(psfile,0,strpos(psfile,'.'))
    psfile=psfile+'.ps'
    device,/landscape,file=psfile
    if(thar_flag eq 1) then begin
      !p.position=[0.05,0.05,0.98,0.25]
      if(keyword_set(wl)) then xr=minmax(wl) else xr=minmax(thar_spec(*,0))
      plot,thar_spec(*,0),thar_spec(*,1)*csscale_thar,yr=[0,max(thar_spec(*,1))] $
          ,col=200,xs=1,xr=xr,yticklen=0.005,ytickname=replicate(' ',10),yticks=2
      if(keyword_set(thar_list)) then begin
        ii=where(thar_list(*,0) ge xr(0) and thar_list(*,0) le xr(1), nii)
        if(nii gt 0) then begin
          ylim=max(thar_spec(*,1))*4.
          xyouts,thar_list(ii,0),(thar_list(ii,1)*csscale_thar)<ylim,$
                 thar_list(ii,0),orient=90,col=200
        endif
      endif
      !p.position=[0.05,0.25,0.98,0.99]
    endif else !p.position=0
    if(keyword_set(wl)) then begin
      plot,wl,cs,xs=1,ys=3,yticklen=0.005,yr=csmin+[0,(csmax-csmin)*2],$
           title=title,noerase=thar_flag
      oplot,wl,cs*csscale(0)-csscale(1)
      if(nline gt 0) then begin
        wlmin=min(wl)
        wlmax=max(wl)
        for i=0,nline-1 do begin
          x=line(i).wll
          if(x ge wlmin and x lt wlmax) then begin
            xx=min(abs(x-wl),ix)
            oplot,[x,x],cs(ix)+[0,0.1*(!y.crange(1)-!y.crange(0))]
            if(line(i).flag eq 0) then $
              xyouts,x,cs(ix)+0.1*(!y.crange(1)-!y.crange(0)), $
                     string(line(i).wll,'(F10.4)'),orientation=90 $
            else $
              xyouts,x,cs(ix)+0.1*(!y.crange(1)-!y.crange(0)), $
                     string(line(i).wll,'(F10.4,''*'')'),orientation=90
          endif
          x=line(i).wlc
          xx=min(abs(x-wl),ix)
          if(x ge wlmin and x lt wlmax) then $
            oplot,[x,x],[cs(ix)+0.1*(!y.crange(1)-!y.crange(0)),!y.crange(1)]
        endfor
      endif
    endif else begin
      plot,cs,xs=1,ys=3,yticklen=0.005,yr=csmin+[0,(csmax-csmin)*2],$
           title=title,noerase=thar_flag
      oplot,cs*csscale(0)-csscale(1)
      if(nline gt 0) then begin
        for i=0,nline-1 do begin
          x=round(line(i).posm)
          if(x ge 0 and x lt n_elements(cs)) then begin
            oplot,[x,x],cs(x)+[0,0.1*(!y.crange(1)-!y.crange(0))]
            if(line(i).flag eq 0) then $
              xyouts,x,cs(x)+0.1*(!y.crange(1)-!y.crange(0)), $
                     string(line(i).wll,'(F10.4)'),orientation=90 $
            else $
              xyouts,x,cs(x)+0.1*(!y.crange(1)-!y.crange(0)), $
                     string(line(i).wll,'(F10.4,''*'')'),orientation=90
          endif
          x=round(line(i).posc)
          if(x ge 0 and x lt n_elements(cs)) then $
            oplot,[x,x],[cs(x)+0.1*(!y.crange(1)-!y.crange(0)),!y.crange(1)]
        endfor
      endif
    endelse
    device,/close
    print,'Figure plotted to the file "'+psfile+'"'
    set_plot,default_device
  endif else begin
    if(thar_flag eq 1) then begin
      !p.position=[0.05,0.05,0.98,0.25]
      if(keyword_set(wl)) then xr=minmax(wl) else xr=minmax(thar_spec(*,0))
      plot,thar_spec(*,0),thar_spec(*,1)*csscale_thar,yr=[0,max(thar_spec(*,1))] $
          ,col=4,xs=1,xr=xr,yticklen=0.005,ytickname=replicate(' ',10),yticks=2
      if(keyword_set(thar_list)) then begin
        ii=where(thar_list(*,0) ge xr(0) and thar_list(*,0) le xr(1), nii)
        if(nii gt 0) then begin
          ylim=max(thar_spec(*,1))*4.
          xyouts,thar_list(ii,0),(thar_list(ii,1)*csscale_thar)<ylim $
                ,thar_list(ii,0),orient=90,col=4
        endif
      endif
      !p.position=[0.05,0.25,0.98,0.99]
    endif else !p.position=0
    if(keyword_set(wl)) then begin
      plot,wl,cs,xs=1,ys=3,col=0,yticklen=0.005,yr=csmin+[0,(csmax-csmin)*2] $
          ,noerase=thar_flag
      oplot,wl,cs*csscale(0)-csscale(1),col=6
      if(nline gt 0) then begin
        wlmin=min(wl)
        wlmax=max(wl)
        for i=0,nline-1 do begin
          x=line(i).wll
          if(x ge wlmin and x lt wlmax) then begin
            xx=min(abs(x-wl),ix)
            oplot,[x,x],cs(ix)+[0,0.1*(!y.crange(1)-!y.crange(0))],col=3
            if(line(i).flag eq 0) then $
              xyouts,x,cs(ix)+0.1*(!y.crange(1)-!y.crange(0)), $
                     string(line(i).wll,'(F10.4)'),orientation=90,col=3 $
            else $
              xyouts,x,cs(ix)+0.1*(!y.crange(1)-!y.crange(0)), $
                    string(line(i).wll,'(F10.4,''*'')'),orientation=90,col=5
          endif
          x=line(i).wlc
          xx=min(abs(x-wl),ix)
          if(x ge wlmin and x lt wlmax) then $
            oplot,[x,x],[cs(ix)+0.1*(!y.crange(1)-!y.crange(0)),!y.crange(1)],col=1
        endfor
      endif
    endif else begin
      plot,cs,xs=1,ys=3,col=0,yticklen=0.005,yr=csmin+[0,(csmax-csmin)*2] $
          ,noerase=thar_flag
      oplot,cs*csscale(0)-csscale(1),col=6
      if(nline gt 0) then begin
        for i=0,nline-1 do begin
          x=round(line(i).posm)
          if(x ge 0 and x lt n_elements(cs)) then begin
            oplot,[x,x],cs(x)+[0,0.1*(!y.crange(1)-!y.crange(0))],col=3
            if(line(i).flag eq 0) then $
              xyouts,x,cs(x)+0.1*(!y.crange(1)-!y.crange(0)), $
                     string(line(i).wll,'(F10.4)'),orientation=90,col=3 $
            else $
             xyouts,x,cs(x)+0.1*(!y.crange(1)-!y.crange(0)), $
                    string(line(i).wll,'(F10.4,''*'')'),orientation=90,col=5
          endif
          x=round(line(i).posc)
          if(x ge 0 and x lt n_elements(cs)) then $
            oplot,[x,x],[cs(x)+0.1*(!y.crange(1)-!y.crange(0)),!y.crange(1)],col=1
        endfor
      endif
    endelse
  endelse
  return
end

pro mark,xfirst,xlast,cs,csscale,pos,wll,col=col,lambda_scale=wl
  if(keyword_set(wl)) then x=wl(indgen(xlast-xfirst+1)+xfirst) $
  else                     x=indgen(xlast-xfirst+1)+xfirst
  oplot,x,cs,col=col,thick=3
  oplot,x,cs*csscale(0)-csscale(1),col=col,thick=3
;  oplot,x,cs,col=1,thick=1
  if(pos ge 0 and pos lt n_elements(cs)) then begin
    if(keyword_set(wl)) then xpos=wl(pos) else xpos=pos
    oplot,[xpos,xpos],cs(pos-xfirst)+[0,0.1*(!y.crange(1)-!y.crange(0))],col=3
    xyouts,xpos,cs(pos-xfirst)+0.1*(!y.crange(1)-!y.crange(0)), $
           string(wll,'(F10.4)'),orientation=90,col=3
  endif
  oplot,x,cs,col=0
  oplot,x,cs*csscale(0)-csscale(1),col=0
  return
end

pro dsp_fit,cs,line,nline,iline,disp_poly,top_event,GROUP_LEADER=wGroup ; Show computed wavelength
;
;  Setup a new widget
;
  line_base=widget_base(title='RMS of the last fit',/COLUMN,/FRAME,GROUP_LEADER=wGroup)
  base=widget_base(line_base,/ROW)
  b1=widget_button(base,value='EXIT ',uvalue='EXIT')
  base1=widget_base(base,/NONEXCLUSIVE)
  b2=widget_button(base1,value='Flag line ',uvalue='FLAG')
  device,get_screen_size=scr_sz
  scr_sz[0]=scr_sz[0]<1280
  XPlotdisplay=widget_draw(line_base,XSIZE=scr_sz(0)-100,YSIZE=scr_sz(1)-200, $
                           /FRAME,/BUTTON_EVENTS,colors=255)
  !p.position=[0.12,0.12,0.88,0.96]

  widget_control, line_base,/REALIZE ; Realize the widget
  if(line(iline).flag) then widget_control,b2,SET_BUTTON=1
  y=line(0:nline-1).wll-line(0:nline-1).wlc
  x=line(0:nline-1).posm
  flag=line(0:nline-1).flag
  erase
  xrr=[0,n_elements(cs)-1]
  plot,poly(dindgen(n_elements(cs)),disp_poly), $ ; Plot the calibration curve
       xs=3,ys=11,col=0,xr=xrr,ytitle='Wavelength'
  plot,x,y,xs=3,ys=-2,col=0,                    $ ; Plot RMS of the calibration curve
       /nodata,/noerase,xr=xrr
  axis,yaxis=1,/noerase, ytitle='Residuals',ys=11,col=0

  oplot,!x.crange,[0,0],col=0
  ; Plot individual lines
  for i=0,nline-1 do begin
    if(flag(i) eq 0) then oplot,[x(i),x(i)],[y(i),y(i)],col=3,psym=4 $
    else                  oplot,[x(i),x(i)],[y(i),y(i)],col=3,psym=1
  endfor
  ; Mark currently selected line with color
  if(line(iline).flag eq 0) then $
     oplot,[x(iline),x(iline)],[y(iline),y(iline)],col=2,psym=4 $
  else $
     oplot,[x(iline),x(iline)],[y(iline),y(iline)],col=2,psym=1

  while(1) do begin
    event=widget_event(line_base,BAD_ID=bad)
    if(bad ne 0) then begin
      posm=-1.d0
      return
    endif
    if(event.id eq XPlotdisplay) then begin
      if(event.press eq 1) then begin   ; Left mouse button pressed
        ; Unmark currently selected line
        if(line(iline).flag eq 0) then  $
          oplot,[x(iline),x(iline)],[y(iline),y(iline)],col=3,psym=4 $
        else $
          oplot,[x(iline),x(iline)],[y(iline),y(iline)],col=3,psym=1
        ; Step back in iline
        iline=iline-1 & if(iline lt 0) then iline=nline-1
        ; Mark newly selected line with red color and correct flag button
        if(line(iline).flag eq 0) then $
          oplot,[x(iline),x(iline)],[y(iline),y(iline)],col=2,psym=4 $
        else $
          oplot,[x(iline),x(iline)],[y(iline),y(iline)],col=2,psym=1
        if(line(iline).flag eq 1) then widget_control,b2,SET_BUTTON=1 $
        else                           widget_control,b2,SET_BUTTON=0
      endif else if(event.press gt 1) then begin   ; Right mouse button pressed
        ; Unmark currently selected line
        if(line(iline).flag eq 0) then $
          oplot,[x(iline),x(iline)],[y(iline),y(iline)],col=3,psym=4 $
        else $
          oplot,[x(iline),x(iline)],[y(iline),y(iline)],col=3,psym=1
        ; Step back in iline
        iline=iline+1 & if(iline ge nline) then iline=0
        ; Mark newly selected line with red color and correct flag button
        if(line(iline).flag eq 0) then begin
          oplot,[x(iline),x(iline)],[y(iline),y(iline)],col=2,psym=4
          widget_control,b2,SET_BUTTON=0
        endif else begin
          oplot,[x(iline),x(iline)],[y(iline),y(iline)],col=2,psym=1
          widget_control,b2,SET_BUTTON=1
        endelse
      endif
    endif else begin
      widget_control,event.id,GET_UVALUE=userid
      case userid of
       'FLAG':begin
                line(iline).flag=(line(iline).flag+1) mod 2
                if(line(iline).flag eq 0) then begin
                  oplot,[x(iline),x(iline)],[y(iline),y(iline)],col=!p.background,psym=1
                  oplot,[x(iline),x(iline)],[y(iline),y(iline)],col=2,psym=4
                  widget_control,b2,SET_BUTTON=0
                endif else begin
                  oplot,[x(iline),x(iline)],[y(iline),y(iline)],col=!p.background,psym=4
                  oplot,[x(iline),x(iline)],[y(iline),y(iline)],col=2,psym=1
                  widget_control,b2,SET_BUTTON=1
                endelse
              end
       'EXIT':begin
                widget_control,line_base,/DESTROY
                return
              end
      endcase
    endelse
  endwhile
  widget_control,line_base,/DESTROY
  return
end

pro dsp_mark,cs,line,nline,iline,wllist,top_event, $ ; Zoom onto selected line
             box_px,box_py,zoom_ON,box_clr,nfont, $  ; Do the line fitting
             default_approx,csscale,csscale_thar,disp_poly, $
             COMPUTED=COMPUTED,thar_flag=thar_flag, $
             thar_list=thar_list,thar_spec=thar_spec, $
             lambda_scale=wl,GROUP_LEADER=wGroup

  posm=-1.d0
  xy=convert_coord(top_event.x,top_event.y,/DEVICE,/TO_DATA)
  if((zoom_ON ne 1) and $                         ; No mouse button pressed
     (top_event.type ne 0)) then return $
  else if(top_event.type eq 0) then begin         ; Now button was pressed
    zoom_ON=1
    device,set_graphics=6                         ; Set XOR
    xx=xy(0) & yy=xy(1)                           ; Select one corner
    box_px(0:4)=xx
    box_py(0:4)=yy
    oplot,box_px,box_py, $                        ; Plot new box
          thick=1,/noclip,lines=0,col=box_clr
  endif else if((zoom_ON eq 1) and $              ; Mouse moved while button
                (top_event.type eq 2)) then begin ; is still pressed
    oplot,box_px,box_py, $                        ; Erase old box
    thick=1,/noclip,lines=0,col=box_clr
    box_px(1:2)=xy(0)                             ; Update box coordinates
    box_py(2:3)=xy(1)
    oplot,box_px,box_py, $                        ; Plot new box
          thick=1,/noclip,lines=0,col=box_clr
  endif else if((zoom_ON eq 1) and $              ; Button was released
                (top_event.type eq 1)) then begin
    zoom_ON=0
    oplot,box_px,box_py, $                        ; Erase old box
          thick=1,/noclip,lines=0,col=box_clr
    device,set_graphics=3                         ; Return to the normal mode
    if(keyword_set(wl)) then begin
      a=min(abs(wl-(box_px(0)<box_px(1))),xfirst)
      a=min(abs(wl-(box_px(0)>box_px(1))),xlast)
    endif else begin
      xfirst=max([0,round(min([box_px(0),box_px(1)]))])
      xlast =min([round(max([box_px(0),box_px(1)])),n_elements(cs)-1])
    endelse
    if(xlast-xfirst lt 5) then return
    ii=where(line.posm ge xfirst and line.posm le xlast, nii)
    if(nii gt 0) then begin
      ilin=ii(0) & iline=ilin & approx=line(ilin).approx
      wll=line(ilin).wll & height=line(ilin).height & width=line(ilin).width
    endif else begin
      ilin=-1 & approx=default_approx
    endelse
;
;  Setup a new widget
;
    line_base=widget_base(title='Selected line',/COLUMN,/FRAME,GROUP_LEADER=wGroup)
    base=widget_base(line_base,/ROW)
    buttons=widget_base(base,/COLUMN)
    b1=widget_button(buttons,value='Add ',uvalue='ADD')
    b2=widget_button(buttons,value='Ignore ',uvalue='IGNORE')
    b3  =widget_base(base,title='Wavelength:',/COLUMN,/FRAME)
    b3a =widget_base(b3,/ROW)
    b3_a=widget_label(b3a,value='Lab.WL:   ',font=nfont)
    isel=0
    if(ilin ge 0) then begin
      if(line(iline).wll le 0 and line(iline).wlc gt 0) then begin
        isel=(where(abs(line(iline).wlc-wllist) eq $
                min(abs(line(iline).wlc-wllist))))(0)
        line(iline).wll=wllist(isel)
      endif
      wll=line(iline).wll
      b3_1=widget_text(b3a,value=string(line(ilin).wll,'(F10.4)'), $
                       uvalue='WLLAB',XSIZE=12,/EDITABLE,/KBRD_FOCUS_EVENTS,font=nfont)
    endif else begin
      b3_1=widget_text(b3a,value='    0.0000',   $
                       uvalue='WLLAB',XSIZE=12,/EDITABLE,/KBRD_FOCUS_EVENTS,font=nfont)
      wll=-1.d0
    endelse
    b3b =widget_base(b3,/ROW)
    b3_b=widget_label(b3b,value='Line list ',font=nfont)
    if(iline ge 0) then ii=sort(abs(line(iline).wll-wllist)) $
    else                ii=indgen(n_elements(wllist)<40)
    nii=n_elements(ii)<40 & ii=ii(0:nii-1)
    wltmplist=wllist(ii)
    isel=(where(sort(wltmplist) eq 0))(0)
    wltmplist=wltmplist(sort(wltmplist))
    b3_2=widget_droplist(b3b,value=string(wltmplist,'(F10.4)'),/DYNAMIC_RESIZE, $
         uvalue='WLLIST',font=nfont)
    b4  =widget_base(base,title='Approximation:',COLUMN=2,/FRAME,/EXCLUSIVE)
    b4_1=widget_button(b4,value='Gauss     ',uvalue='APPGAU')
    b4_2=widget_button(b4,value='Parabola  ',uvalue='APPPAR')
    b4_3=widget_button(b4,value='Box+Gauss ',uvalue='APPCEN')
    b4_4=widget_button(b4,value='Lorentz   ',uvalue='APPLOR')
    XPlotdisplay =widget_draw(line_base,XSIZE=400,YSIZE=200,/FRAME,colors=255)
    !p.position=[0.12,0.12,0.96,0.96]

    widget_control, line_base,/REALIZE ; Realize the widget

    widget_control,b3_2,set_droplist_select=isel

    erase

    case approx of ; Compute line center position
                   ; according to the default approximation
     'G':begin
           x=xfirst+indgen(xlast-xfirst+1)
           xapp=xfirst+dindgen((xlast-xfirst)*10+1)/10
           yapp=disp_gaussfit(x,cs(xfirst:xlast),a,xapp)
           height=abs(a(0))
           posm=a(1)
           if(keyword_set(computed)) then begin
             wlc=poly(posm,disp_poly)
             wll=min(abs(wlc-wllist),isel1)
;             isel1=(where(abs(wlc-wllist) eq min(abs(wlc-wllist))))(0)
             wll=wllist(isel1)
             ii=sort(abs(wll-wllist))
             nii=n_elements(ii)<40 & ii=ii(0:nii-1)
             wltmplist=wllist(ii)
             isel=(where(sort(wltmplist) eq 0))(0)
             wltmplist=wltmplist(sort(wltmplist))
               widget_control,b3_2,set_value=string(wltmplist,'(F10.4)')
             widget_control,b3_2,set_droplist_select=isel
             widget_control,b3_1,set_value=string(wll,'(F10.4)')
           endif else wlc=-1.
           width=abs(a(2))
           height=abs(a(0))
           yr=[min(cs(xfirst:xlast)),max(cs(xfirst:xlast))>max(yapp)]
           plot,x,cs(xfirst:xlast),col=0,/NODATA,ys=3,yr=yr
           oplot,x,cs(xfirst:xlast),col=1
           oplot,xapp,yapp,col=3
           oplot,[posm,posm],yr,col=5
           widget_control,b4_1,/SET_BUTTON
         end
     'P':begin
           x=xfirst+indgen(xlast-xfirst+1)
           xapp=xfirst+dindgen((xlast-xfirst)*10+1)/10
           a=poly_fit(x-xfirst,cs(xfirst:xlast),2)
           yapp=poly(xapp-xfirst,a)
           posm=-a(1)/(2*a(2))+xfirst
           if(keyword_set(computed)) then begin
             wlc=poly(posm,disp_poly)
             wll=min(abs(wlc-wllist),isel1)
;             isel1=(where(abs(wlc-wllist) eq min(abs(wlc-wllist))))(0)
             wll=wllist(isel1)
             ii=sort(abs(wll-wllist))
             nii=n_elements(ii)<40 & ii=ii(0:nii-1)
             wltmplist=wllist(ii)
             isel=(where(sort(wltmplist) eq 0))(0)
             wltmplist=wltmplist(sort(wltmplist))
             widget_control,b3_2,set_value=string(wltmplist,'(F10.4)')
             widget_control,b3_2,set_droplist_select=isel
             widget_control,b3_1,set_value=string(wll,'(F10.4)')
           endif else wlc=-1.
           height=poly(posm-xfirst,a)-a(0)
           width=(xlast-xfirst)*0.5
           yr=[min(cs(xfirst:xlast)),max(cs(xfirst:xlast))>max(yapp)]
           plot,x,cs(xfirst:xlast),col=0,/NODATA,ys=3,yr=yr
           oplot,x,cs(xfirst:xlast),col=1
           oplot,xapp,yapp,col=3
           oplot,[posm,posm],yr,col=5
           widget_control,b4_2,/SET_BUTTON
         end
     'C':begin
           x=xfirst+indgen(xlast-xfirst+1)
           xapp=xfirst+dindgen((xlast-xfirst)*10+1)/10
           yapp=disp_gaussboxfit(x,cs(xfirst:xlast),a,xapp)
           height=abs(a(0))
           posm=a(1)
           if(keyword_set(computed)) then begin
             wlc=poly(posm,disp_poly)
             wll=min(abs(wlc-wllist),isel1)
;             isel1=(where(abs(wlc-wllist) eq min(abs(wlc-wllist))))(0)
             wll=wllist(isel1)
             ii=sort(abs(wll-wllist))
             nii=n_elements(ii)<40 & ii=ii(0:nii-1)
             wltmplist=wllist(ii)
             isel=(where(sort(wltmplist) eq 0))(0)
             wltmplist=wltmplist(sort(wltmplist))
             widget_control,b3_2,set_value=string(wltmplist,'(F10.4)')
             widget_control,b3_2,set_droplist_select=isel
             widget_control,b3_1,set_value=string(wll,'(F10.4)')
           endif else wlc=-1.
           width=abs(a(2))
           yr=[min(cs(xfirst:xlast)),max(cs(xfirst:xlast))>max(yapp)]
           plot,x,cs(xfirst:xlast),col=0,/NODATA,ys=3,yr=yr
           oplot,x,cs(xfirst:xlast),col=1
           oplot,xapp,yapp,col=3
           oplot,[posm,posm],yr,col=5
           widget_control,b4_3,/SET_BUTTON
         end
     'L':begin
           x=xfirst+indgen(xlast-xfirst+1)
           xapp=xfirst+dindgen((xlast-xfirst)*10+1)/10
           yapp=disp_lorentzfit(x,cs(xfirst:xlast),a,xapp)
           height=a(0)/(a(2)*a(2))
           posm=a(1)
           if(keyword_set(computed)) then begin
             wlc=poly(posm,disp_poly)
             wll=min(abs(wlc-wllist),isel1)
;             isel1=(where(abs(wlc-wllist) eq min(abs(wlc-wllist))))(0)
             wll=wllist(isel1)
             ii=sort(abs(wll-wllist))
             nii=n_elements(ii)<40 & ii=ii(0:nii-1)
             wltmplist=wllist(ii)
             isel=(where(sort(wltmplist) eq 0))(0)
             wltmplist=wltmplist(sort(wltmplist))
             widget_control,b3_2,set_value=string(wltmplist,'(F10.4)')
             widget_control,b3_2,set_droplist_select=isel
             widget_control,b3_1,set_value=string(wll,'(F10.4)')
           endif else wlc=-1.
           width=abs(a(2))
           yr=[min(cs(xfirst:xlast)),max(cs(xfirst:xlast))>max(yapp)]
           plot,x,cs(xfirst:xlast),col=0,/NODATA,ys=3,yr=yr
           oplot,x,cs(xfirst:xlast),col=1
           oplot,xapp,yapp,col=3
           oplot,[posm,posm],yr,col=5
           widget_control,b4_4,/SET_BUTTON
         end
    endcase

    while(1) do begin
      event=widget_event(line_base,BAD_ID=bad)
      if(bad ne 0) then begin
        posm=-1.d0
        return
      endif
      widget_control,event.id,GET_UVALUE=userid
      case userid of
        'ADD':begin
                if(ilin ge 0) then iline=ilin else iline=nline
                widget_control,line_base,/DESTROY
                line(iline).wll=wll
                line(iline).wlc=wlc
                line(iline).posm=posm
                line(iline).xfirst=xfirst
                line(iline).xlast =xlast
                line(iline).approx=approx
                line(iline).width=width
                line(iline).height=height
                line(iline).flag=0
                if(ilin lt 0) then begin
                  ii=sort(line.posm)
                  line=line(ii)
                  iline=(where(ii eq iline))(0)
                  nline=nline+1
                  line=[line,{csline,-1.d0,-1.d0,-1.d0,-1.d0,0,0,'G',0.d0,0,0.d0,0}]
                endif else begin
                  ii=sort(line(0:nline-1).posm)
                  line=[line(ii),{csline,-1.d0,-1.d0,-1.d0,-1.d0,0,0,'G',0.d0,0,0.d0,0}]
                  iline=(where(ii eq iline))(0)
                endelse
                ncs=n_elements(cs)
                !p.position=[80.d0/ncs,0.12,1-40.d0/ncs,0.96]
                !p.region  =[0.00,0.000,1.00,1.00]
                draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
                       ,thar_flag=thar_flag,thar_list=thar_list $
                       ,lambda_scale=wl
                mark,line(iline).xfirst,line(iline).xlast,     $
                     cs(line(iline).xfirst:line(iline).xlast), $
                     csscale,lambda_scale=wl,                  $
                     round(line(iline).posm),line(iline).wll,col=1
;                widget_control,line_base,/DESTROY
                return
              end
     'IGNORE':begin
                widget_control,line_base,/DESTROY
                ncs=n_elements(cs)
                !p.position=[80.d0/ncs,0.12,1-40.d0/ncs,0.96]
                !p.region  =[0.00,0.000,1.00,1.00]
                draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
                       ,thar_flag=thar_flag,thar_list=thar_list $
                       ,lambda_scale=wl
                if(nline gt 0) then $
                  mark,line(iline).xfirst,line(iline).xlast,     $
                       cs(line(iline).xfirst:line(iline).xlast), $
                       csscale,lambda_scale=wl,                  $
                       round(line(iline).posm),line(iline).wll,col=1
;                widget_control,line_base,/DESTROY
                return
              end
      'WLLAB':begin
                widget_control,b3_1,get_value=www
                wll=double(www(0))
                ii=sort(abs(wll-wllist))
                nii=n_elements(ii)<40 & ii=ii(0:nii-1)
                wltmplist=wllist(ii)
                isel=(where(sort(wltmplist) eq 0))(0)
                wltmplist=wltmplist(sort(wltmplist))
                widget_control,b3_2,set_value=string(wltmplist,'(F10.4)')
                widget_control,b3_2,set_droplist_select=isel
              end
     'WLLIST':begin
                wll=wltmplist(event.index)
                ii=sort(abs(wll-wllist))
                nii=n_elements(ii)<40 & ii=ii(0:nii-1)
                wltmplist=wllist(ii)
                isel=(where(sort(wltmplist) eq 0))(0)
                wltmplist=wltmplist(sort(wltmplist))
                widget_control,b3_2,set_value=string(wltmplist,'(F10.4)')
                widget_control,b3_2,set_droplist_select=isel
                widget_control,b3_1,set_value=string(wll,'(F10.4)')
              end
     'APPGAU':if(event.select eq 1) then begin
                x=xfirst+indgen(xlast-xfirst+1)
                xapp=xfirst+dindgen((xlast-xfirst)*10+1)/10
                yapp=disp_gaussfit(x,cs(xfirst:xlast),a,xapp)
                height=abs(a(0))
                posm=a(1)
                if(keyword_set(computed)) then begin
                  wlc=poly(posm,disp_poly)
                  wll=min(abs(wlc-wllist),isel1)
;                  isel1=(where(abs(wlc-wllist) eq min(abs(wlc-wllist))))(0)
                  wll=wllist(isel1)
                  ii=sort(abs(wll-wllist))
                  nii=n_elements(ii)<40 & ii=ii(0:nii-1)
                  wltmplist=wllist(ii)
                  isel=(where(sort(wltmplist) eq 0))(0)
                  wltmplist=wltmplist(sort(wltmplist))
                  widget_control,b3_2,set_value=string(wltmplist,'(F10.4)')
                  widget_control,b3_2,set_droplist_select=isel
                  widget_control,b3_1,set_value=string(wll,'(F10.4)')
                endif
                width=abs(a(2))
                yr=[min(cs(xfirst:xlast)),max(cs(xfirst:xlast))>max(yapp)]
                plot,x,cs(xfirst:xlast),col=0,/NODATA,ys=3,yr=yr
                oplot,x,cs(xfirst:xlast),col=1
                oplot,xapp,yapp,col=3
                oplot,[posm,posm],yr,col=5
                approx='G'
              end
     'APPLOR':if(event.select eq 1) then begin
                x=xfirst+indgen(xlast-xfirst+1)
                xapp=xfirst+dindgen((xlast-xfirst)*10+1)/10
                yapp=disp_lorentzfit(x,cs(xfirst:xlast),a,xapp)
                height=a(0)/(a(2)*a(2))
                posm=a(1)
                if(keyword_set(computed)) then begin
                  wlc=poly(posm,disp_poly)
                  wll=min(abs(wlc-wllist),isel1)
;                  isel1=(where(abs(wlc-wllist) eq min(abs(wlc-wllist))))(0)
                  wll=wllist(isel1)
                  ii=sort(abs(wll-wllist))
                  nii=n_elements(ii)<40 & ii=ii(0:nii-1)
                  wltmplist=wllist(ii)
                  isel=(where(sort(wltmplist) eq 0))(0)
                  wltmplist=wltmplist(sort(wltmplist))
                  widget_control,b3_2,set_value=string(wltmplist,'(F10.4)')
                  widget_control,b3_2,set_droplist_select=isel
                  widget_control,b3_1,set_value=string(wll,'(F10.4)')
                endif
                width=abs(a(2))
                yr=[min(cs(xfirst:xlast)),max(cs(xfirst:xlast))>max(yapp)]
                plot,x,cs(xfirst:xlast),col=0,/NODATA,ys=3,yr=yr
                oplot,x,cs(xfirst:xlast),col=1
                oplot,xapp,yapp,col=3
                oplot,[posm,posm],yr,col=5
                approx='L'
              end
     'APPPAR':if(event.select eq 1) then begin
                x=xfirst+indgen(xlast-xfirst+1)
                xapp=xfirst+dindgen((xlast-xfirst)*10+1)/10
                a=poly_fit(x-xfirst,cs(xfirst:xlast),2)
                height=poly(posm-xfirst,a)-a(0)
                width=(xlast-xfirst)*0.5
                yapp=poly(xapp-xfirst,a)
                posm=-a(1)/(2*a(2))+xfirst
                if(keyword_set(computed)) then begin
                  wlc=poly(posm,disp_poly)
                  wll=min(abs(wlc-wllist),isel1)
;                  isel1=(where(abs(wlc-wllist) eq min(abs(wlc-wllist))))(0)
                  wll=wllist(isel1)
                  ii=sort(abs(wll-wllist))
                  nii=n_elements(ii)<40 & ii=ii(0:nii-1)
                  wltmplist=wllist(ii)
                  isel=(where(sort(wltmplist) eq 0))(0)
                  wltmplist=wltmplist(sort(wltmplist))
                  widget_control,b3_2,set_value=string(wltmplist,'(F10.4)')
                  widget_control,b3_2,set_droplist_select=isel
                  widget_control,b3_1,set_value=string(wll,'(F10.4)')
                endif
                yr=[min(cs(xfirst:xlast)),max(cs(xfirst:xlast))>max(yapp)]
                plot,x,cs(xfirst:xlast),col=0,/NODATA,ys=3,yr=yr
                oplot,x,cs(xfirst:xlast),col=1
                oplot,xapp,yapp,col=3
                oplot,[posm,posm],yr,col=5
                approx='P'
              end
     'APPCEN':if(event.select eq 1) then begin
                x=xfirst+indgen(xlast-xfirst+1)
                xapp=xfirst+dindgen((xlast-xfirst)*10+1)/10
                yapp=disp_gaussboxfit(x,cs(xfirst:xlast),a,xapp)
                height=abs(a(0))
                posm=a(1)
                if(keyword_set(computed)) then begin
                  wlc=poly(posm,disp_poly)
                  wll=min(abs(wlc-wllist),isel1)
;                  isel1=(where(abs(wlc-wllist) eq min(abs(wlc-wllist))))(0)
                  wll=wllist(isel1)
                  ii=sort(abs(wll-wllist))
                  nii=n_elements(ii)<40 & ii=ii(0:nii-1)
                  wltmplist=wllist(ii)
                  isel=(where(sort(wltmplist) eq 0))(0)
                  wltmplist=wltmplist(sort(wltmplist))
                  widget_control,b3_2,set_value=string(wltmplist,'(F10.4)')
                  widget_control,b3_2,set_droplist_select=isel
                  widget_control,b3_1,set_value=string(wll,'(F10.4)')
                endif
                width=abs(a(2))
                yr=[min(cs(xfirst:xlast)),max(cs(xfirst:xlast))>max(yapp)]
                plot,x,cs(xfirst:xlast),col=0,/NODATA,ys=3,yr=yr
                oplot,x,cs(xfirst:xlast),col=1
                oplot,xapp,yapp,col=3
                oplot,[posm,posm],yr,col=5
                approx='C'
              end
      endcase
    endwhile
    widget_control,line_base,/DESTROY
    draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
        ,thar_flag=thar_flag,thar_list=thar_list,lambda_scale=wl
  endif
  return
end

function dsp_widget,cs,wllist,font,nfont,default_device,box_clr, $
                    line_id=line,log=log,poly_power=polynom_power, $
                    thar_spec=thar_spec,thar_list=thar_list, $
                    widget_title=widget_title,GROUP_LEADER=wGroup

  is_computed=0  ; Flag showing the the dispersion curve has been computed
  if(not keyword_set(line)) then begin
    line={csline, wlc:-1.d0,  $ ; computed wavelength
                  wll:-1.d0,  $ ; laboratory wavelength
                  posc:-1.d0, $ ; computed position
                  posm:-1.d0, $ ; measured position
                  xfirst:0,   $ ; starting pixel
                  xlast:0,    $ ; ending pixel
                  approx:'G', $ ; approximation type
                  width:0.d0, $ ; width in pixels
                  flag:0,     $ ; flag status
                  height:0.d0,$ ; line strength
                  order:0}      ; spectral order
    line=[line]
    nline=0
    iline=-1
  endif else begin
    nline=n_elements(line)
    iline=0
    line=[line,{csline,-1.d0,-1.d0,-1.d0,-1.d0,0,0,'G',0.d0,0,0.d0,0}]
  endelse

  device,get_screen_size=scr_sz
  scr_sz[1]=scr_sz[1]-30
  widget_control,DEFAULT_FONT=font
  if(keyword_set(widget_title)) then widtit=widget_title else widtit='Fit the dispersion curve'
  disp_curve=widget_base(TITLE=widtit,/COLUMN $
                         ,GROUP_LEADER=wGroup,XSIZE=scr_sz[0]-5)
  buttons=widget_base(disp_curve,/ROW)

  b6  =widget_base(buttons,/COLUMN)
  b6_1=widget_button(b6,value='EXIT',uvalue='EXIT')
  b6_2=widget_button(b6,value='Compute  ',uvalue='COMP')

  b1  =widget_base(buttons,/COLUMN)
  b1a =widget_base(b1,/NONEXCLUSIVE)
  b1_1=widget_button(b1a,value='Flag line  ',uvalue='FLAG')
  b1_2=widget_button(b1,value='Remove line ',uvalue='REM')
  b1b =widget_base(b1,/NONEXCLUSIVE)
  b1_3=widget_button(b1b,value='Show ThAr   ',uvalue='SHOW')

  b2  =widget_base(buttons,/COLUMN)
  b2_1=widget_button(b2,value='Next line ',uvalue='NEXT')
  b2_2=widget_button(b2,value='Prev.line ',uvalue='PREV')
  b2_3=widget_button(b2,value=' -> Wave  ',uvalue='PIXWAVE')

  b3  =widget_base(buttons,title='Approximation:',COLUMN=2,/FRAME,/EXCLUSIVE)
  b3_1=widget_button(b3,value='Gauss     ',uvalue='APPGAU')
  b3_2=widget_button(b3,value='Box+Gauss ',uvalue='APPCEN')
  b3_3=widget_button(b3,value='Parabola  ',uvalue='APPPAR')
  b3_4=widget_button(b3,value='Lorentz   ',uvalue='APPLOR')

  b4  =widget_base(buttons,COLUMN=3,/FRAME)
  b4a =widget_base(b4,/ROW)
  b4_1a=widget_label(b4a,value='WL lab ',font=nfont)
  b4_1=widget_text(b4a,value='    0.0000', uvalue='WLLAB',/EDITABLE $
                      ,/KBRD_FOCUS_EVENTS,xsize=10,font=nfont)
  b4b =widget_base(b4,/ROW)
  b4_2a=widget_label(b4b,value='WL comp',font=nfont)
  b4_2=widget_text(b4b, value='    0.0000',uvalue='WLCOM',xsize=10, font=nfont)
  b4c =widget_base(b4,/ROW)
  b4_3a=widget_label(b4c,value='Meas.pos.',font=nfont)
  b4_3=widget_text(b4c,value='    0.0000', uvalue='MESPOS',xsize=10,font=nfont)
  b4d =widget_base(b4,/ROW)
  b4_4a=widget_label(b4d,value='Comp.pos.',font=nfont)
  b4_4=widget_text(b4d,value='    0.0000', uvalue='COMPOS',xsize=10,font=nfont)
  b4e =widget_base(b4,/ROW)
  b4_5a=widget_label(b4e,value='Width ',font=nfont)
  b4_5=widget_text(b4e,value=' 0.0000', uvalue='WIDTH',xsize=8,font=nfont)
  b4f =widget_base(b4,/ROW)
  b4_6a=widget_label(b4f,value='Line #',font=nfont)
  b4_6=widget_text(b4f, value='  0',uvalue='ILINE',xsize=8, font=nfont)

  b5  =widget_base(buttons,/COLUMN)
  b5_1=widget_label(b5,value='Polynom ')
  b5_2=widget_label(b5,value='power   ')
  b5_3=widget_droplist(b5,value=['1','2','3','4','5','6','7'],uvalue='POWER')

  bmsg=widget_base(disp_curve,/ROW)
  b7_0=widget_button(bmsg,value='Auto ID',uvalue='AUTOID')
  b7_1=widget_button(bmsg,value='Redetermine centers',uvalue='REDO')
  b7_2=widget_button(bmsg,value='Print',uvalue='PRINT')
  b7_3=widget_label(bmsg,value='     Selected lines: ',/align_right)
  b7_4=widget_text(bmsg,value='  0',font=nfont,xsize=4)
  b7_5=widget_label(bmsg,value='     Flagged lines: ',/align_right)
  b7_6=widget_text(bmsg,value='  0',font=nfont,xsize=4)
  b7_7=widget_label(bmsg,value='     Mean deviation: ',/align_right)
  b7_8=widget_text(bmsg,value='-100.000000 A',font=nfont,xsize=14)

  b8  =widget_base(disp_curve,/ROW) ; Wavelength shift slider
  sldr=widget_slider(b8,MINIMUM=-410,MAXIMUM=410,XSIZE=scr_sz(0)-50,$
       value=0,uvalue='SHIFT',/SUPPRESS_VALUE)
  pos_shift0=0.d0

  scly=widget_slider(buttons,MINIMUM=1,MAXIMUM=100,YSIZE=100,$
       value=100,uvalue='SCALE',/SUPPRESS_VALUE,/VERTICAL)

  sclthar=widget_slider(buttons,MINIMUM=1,MAXIMUM=100,YSIZE=100,$
       value=10,uvalue='SCALEThAr',/SUPPRESS_VALUE,/VERTICAL)

  device,decomposed=0
  tvlct,r,g,b,/GET
  BLACK=0 & RED=1 & GREEN=2 & BLUE=3 & WHITE=n_elements(r)-1
  r(BLACK)=  0B & r(RED)=255B & r(GREEN)=  0B & r(BLUE)=  0B & r(WHITE)=255B
  g(BLACK)=  0B & g(RED)=  0B & g(GREEN)=255B & g(BLUE)=  0B & g(WHITE)=255B
  b(BLACK)=  0B & b(RED)=  0B & b(GREEN)=  0B & b(BLUE)=255B & b(WHITE)=255B
  tvlct,r,g,b

  !p.background=!d.n_colors-1
  ncs=n_elements(cs)
  xfullsize=round(ncs*1.3)
  xviewport=scr_sz(0)-50
  xfullsize=round(ncs*1.3)>xviewport
  XPlotdisplay =WIDGET_DRAW(disp_curve,/SCROLL,RETAIN=2,          $
                            XSIZE=xfullsize,                      $
                            YSIZE=scr_sz(1)-300,                  $
                            X_SCROLL_SIZE=xviewport,              $
                            Y_SCROLL_SIZE=scr_sz(1)-300,          $
                            /MOTION_EVENTS,/BUTTON_EVENTS,/FRAME, $
                            colors=255,                           $
                            uvalue='MARK')

  widget_control, disp_curve,/REALIZE ; Realize the widget
  widget_control,XPlotdisplay,get_value=temp & plot_window=temp
  device,set_graphics=3
  erase
  !p.position=[80.d0/ncs,0.12,1-40.d0/ncs,0.96]
  !p.region  =[0.00,0.000,1.00,1.00]

  if(not keyword_set(thar_spec)) then begin
    widget_control,b1_3,SENSITIVE=0
    widget_control,sclthar,SENSITIVE=0
    Show_ThAr_flag=0
    widget_control,sclthar,SENSITIVE=0
  endif else begin
    Show_ThAr_flag=1
    widget_control,b1_3,SET_BUTTON=1
    widget_control,sclthar,SENSITIVE=1
  endelse
  widget_control,b2_3,SENSITIVE=0
  wl_scale=0

  csscale_thar=10.
  draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
      ,thar_flag=Show_ThAr_flag,thar_list=thar_list $
      ,lambda_scale=wl_scale
  csscale0=csscale
  csmooth=middle(cs,1.d3) ; bottom envelope of cs used to identify the edges
                          ; of sp. lines

  if(iline ge 0) then mark,line(iline).xfirst,line(iline).xlast,     $
                      cs(line(iline).xfirst:line(iline).xlast),      $
                      csscale,lambda_scale=wl_scale,                 $
                      round(line(iline).posm),line(iline).wll,col=RED


  widget_control,b4_1,SET_BUTTON=0
  widget_control,b1_1,SENSITIVE=0
  if(not keyword_set(polynom_power)) then polynom_power=2
  polynom_power=(polynom_power>1)<7
  widget_control,b5_3,SET_DROPLIST_SELECT=polynom_power-1
  widget_control,b7_0,SENSITIVE=0
  zoom_ON=0
  box_px=dindgen(5) & box_py=box_px
  default_approx='G'
  flag=0
  error_message=''
;
;  Loop through the events
;
  disp_poly=double([0,1])
  userid=''
  while(1) do begin
    if(iline ge 0) then begin
      widget_control,b1_1,set_button=line(iline).flag
      widget_control,b4_1,set_value=string(line(iline).wll,'(F10.4)')
      widget_control,b4_2,set_value=string(line(iline).wlc,'(F10.4)')
      widget_control,b4_3,set_value=string(line(iline).posm,'(F10.4)')
      widget_control,b4_4,set_value=string(line(iline).posc,'(F10.4)')
      case line(iline).approx of
        'G':begin
              widget_control,b3_1,SET_BUTTON=1
              widget_control,b3_2,SET_BUTTON=0
              widget_control,b3_3,SET_BUTTON=0
              widget_control,b3_4,SET_BUTTON=0
              default_approx='G'
            end
        'C':begin
              widget_control,b3_1,SET_BUTTON=0
              widget_control,b3_2,SET_BUTTON=1
              widget_control,b3_3,SET_BUTTON=0
              widget_control,b3_4,SET_BUTTON=0
              default_approx='C'
            end
        'P':begin
              widget_control,b3_1,SET_BUTTON=0
              widget_control,b3_2,SET_BUTTON=0
              widget_control,b3_3,SET_BUTTON=1
              widget_control,b3_4,SET_BUTTON=0
              default_approx='P'
            end
        'L':begin
              widget_control,b3_1,SET_BUTTON=0
              widget_control,b3_2,SET_BUTTON=0
              widget_control,b3_3,SET_BUTTON=0
              widget_control,b3_4,SET_BUTTON=1
              default_approx='L'
            end
      endcase
      widget_control,b7_4,set_value=string(nline,'(I2)')
      nflaggged=0
      if(nline gt 0) then i=where(line(0:nline-1).flag eq 1, nflagged)
      widget_control,b7_6,set_value=string(nflagged,'(I3)')
    endif
    if(nline gt 0) then widget_control,b1_1,SENSITIVE=1
    event=widget_event(disp_curve,BAD_ID=bad)
    if(bad ne 0) then begin
;      widget_control,disp_curve,/DESTROY
      loadct,0
      !p.position=0
      !p.background=0
      return,disp_poly
    endif
    widget_control,event.id,GET_UVALUE=userid
    case userid of
      'EXIT':begin
               answ=widget_message('Do you really want to exit?',TITLE='Disp', $
                                   /QUESTION,/CANCEL)
               if(answ eq 'Yes') then begin
                 widget_control,disp_curve,/DESTROY
;                 for i=0,nline-1 do print,form='(4f10.4,2i6,4x,A1,i3)',line(i)
;                 print,iline
                 loadct,0
                 !p.position=0
                 !p.background=0
                 return,disp_poly
               endif else if(answ eq 'Cancel') then begin
                 widget_control,disp_curve,/DESTROY
                 loadct,0
                 !p.position=0
                 !p.background=0
                 return,-1
               endif else widget_control,b5_1,SET_BUTTON=0
             end
      'MARK':begin
               dsp_mark,cs,line,nline,iline,wllist,event, $
                        box_px,box_py,zoom_ON,box_clr,nfont,default_approx, $
                        csscale,csscale_thar,disp_poly,COMPUTED=is_computed, $
                        thar_flag=Show_ThAr_flag,thar_list=thar_list, $
                        thar_spec=thar_spec,lambda_scale=wl_scale,GROUP_LEADER=disp_curve
               if(nline gt 0) then begin
                 widget_control,b4_5,set_value=string(line(iline).width,'(f7.4)')
                 widget_control,b4_6,set_value=string(iline,'(I7)')
               endif
             end
       'REM':if(nline gt 0) then begin
               if(iline eq nline) then begin
                 nline=nline-1
                 iline=iline-1
                 line=line(0:nline)
               endif else if(iline gt 0) then begin
                 line=[line(0:iline-1),line(iline+1:nline)]
                 nline=nline-1
                 iline=iline-1
               endif else begin
                 line=[line(1:nline)]
                 nline=nline-1
               endelse
               draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
                   ,thar_flag=Show_ThAr_flag,thar_list=thar_list $
                   ,lambda_scale=wl_scale
               if(nline gt 0) then begin
                 if(line(iline).xfirst ge 0 and                   $
                    line(iline).xlast  lt n_elements(cs)) then    $
                   mark,line(iline).xfirst,line(iline).xlast,     $
                        cs(line(iline).xfirst:line(iline).xlast), $
                        csscale,lambda_scale=wl_scale,            $
                        round(line(iline).posm),line(iline).wll,col=RED
                 if(line(iline).flag eq 1) then widget_control,b1_1,SET_BUTTON=1 $
                 else                           widget_control,b1_1,SET_BUTTON=0
               endif
               if(n_elements(wl_scale) eq 1) then $
                 xy=convert_coord(line(iline).posm,0,/DATA,/TO_DEVICE) $
               else $
                 xy=convert_coord(line(iline).wll,0,/DATA,/TO_DEVICE)
               widget_control,XPlotdisplay,get_draw_view=xview
               widget_control,XPlotdisplay, $
                              set_draw_view=[(0>(xy(0)-xviewport*0.5))< $
                                                (xfullsize-xviewport-20),0]
               widget_control,b4_5,set_value=string(line(iline).width,'(f7.4)')
               widget_control,b4_6,set_value=string(iline,'(I7)')
             endif
      'NEXT':if(nline gt 0 and iline lt nline-1) then begin
               if(line(iline).xfirst ge 0 and                   $
                  line(iline).xlast  lt n_elements(cs)) then    $
                 mark,line(iline).xfirst,line(iline).xlast,     $
                      cs(line(iline).xfirst:line(iline).xlast), $
                      csscale,lambda_scale=wl_scale,            $
                      round(line(iline).posm),line(iline).wll,col=!p.background
               iline=iline+1
               if(line(iline).xfirst ge 0 and                   $
                  line(iline).xlast  lt n_elements(cs)) then    $
                 mark,line(iline).xfirst,line(iline).xlast,     $
                      cs(line(iline).xfirst:line(iline).xlast), $
                      csscale,lambda_scale=wl_scale,            $
                      round(line(iline).posm),line(iline).wll,col=RED
               if(line(iline).flag eq 1) then widget_control,b1_1,SET_BUTTON=1 $
               else                           widget_control,b1_1,SET_BUTTON=0
               if(n_elements(wl_scale) eq 1) then $
                 xy=convert_coord(line(iline).posm,0,/DATA,/TO_DEVICE) $
               else $
                 xy=convert_coord(line(iline).wll,0,/DATA,/TO_DEVICE)
               widget_control,XPlotdisplay,get_draw_view=xview
               widget_control,XPlotdisplay, $
                              set_draw_view=[(0>(xy(0)-xviewport*0.5))< $
                                                (xfullsize-xviewport-20),0]
               widget_control,b4_5,set_value=string(line(iline).width,'(f7.4)')
               widget_control,b4_6,set_value=string(iline,'(I7)')
             endif
      'PREV':if(iline gt 0)  then begin
               if(line(iline).xfirst gt 0 and $
                  line(iline).xlast lt ncs-1) then $
                 mark,line(iline).xfirst,line(iline).xlast,     $
                      cs(line(iline).xfirst:line(iline).xlast), $
                      csscale,lambda_scale=wl_scale,            $
                      round(line(iline).posm),line(iline).wll,  $
                      col=!p.background
               iline=iline-1
               if(line(iline).xfirst gt 0 and $
                  line(iline).xlast lt ncs-1) then $
                 mark,line(iline).xfirst,line(iline).xlast,     $
                      cs(line(iline).xfirst:line(iline).xlast), $
                      csscale,lambda_scale=wl_scale,            $
                      round(line(iline).posm),line(iline).wll,col=RED
               if(line(iline).flag eq 1) then widget_control,b1_1,SET_BUTTON=1 $
               else                           widget_control,b1_1,SET_BUTTON=0
               if(n_elements(wl_scale) eq 1) then $
                 xy=convert_coord(line(iline).posm,0,/DATA,/TO_DEVICE) $
               else $
                 xy=convert_coord(line(iline).wll,0,/DATA,/TO_DEVICE)
               widget_control,XPlotdisplay,get_draw_view=xview
               widget_control,XPlotdisplay, $
                              set_draw_view=[(0>(xy(0)-xviewport*0.5))< $
                                                (xfullsize-xviewport-20),0]
               widget_control,b4_5,set_value=string(line(iline).width,'(f7.4)')
               widget_control,b4_6,set_value=string(iline,'(I7)')
             endif
      'FLAG':if(iline ge 0) then begin
               line(iline).flag=(line(iline).flag+1) mod 2
               if(line(iline).flag eq 1) then widget_control,b1_1,SET_BUTTON=1 $
               else                           widget_control,b1_1,SET_BUTTON=0
               draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
                   ,thar_flag=Show_ThAr_flag,thar_list=thar_list $
                   ,lambda_scale=wl_scale
               if(nline gt 0) then begin
                 if(line(iline).xfirst ge 0 and                   $
                    line(iline).xlast  lt n_elements(cs)) then    $
                   mark,line(iline).xfirst,line(iline).xlast,     $
                        cs(line(iline).xfirst:line(iline).xlast), $
                        csscale,lambda_scale=wl_scale,            $
                        round(line(iline).posm),line(iline).wll,col=RED
               endif
             endif
      'SHOW':begin
               if(keyword_set(thar_list)) then begin
                 Show_ThAr_flag=(Show_ThAr_flag+1) mod 2
               endif
               if(Show_ThAr_flag eq 1) then begin
                 widget_control,b1_3,SET_BUTTON=1
                 widget_control,sclthar,SENSITIVE=1
               endif else begin
                 widget_control,b1_3,SET_BUTTON=0
                 widget_control,sclthar,SENSITIVE=0
               endelse
               draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
                   ,thar_flag=Show_ThAr_flag,thar_list=thar_list $
                   ,lambda_scale=wl_scale
               if(nline gt 0) then begin
                 if(line(iline).xfirst ge 0 and                   $
                    line(iline).xlast  lt n_elements(cs)) then    $
                   mark,line(iline).xfirst,line(iline).xlast,     $
                        cs(line(iline).xfirst:line(iline).xlast), $
                        csscale,lambda_scale=wl_scale,            $
                        round(line(iline).posm),line(iline).wll,col=RED
               endif
             end
   'PIXWAVE':begin
               if(n_elements(wl_scale) eq 1 and is_computed eq 1) then begin
                 wl_scale=poly(dindgen(n_elements(cs)),disp_poly)
                 widget_control,b2_3,set_value=' -> Pixel '
               endif else begin
                 wl_scale=0
                 widget_control,b2_3,set_value=' -> Wave  '
               endelse
               draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
                   ,thar_flag=Show_ThAr_flag,thar_list=thar_list $
                   ,lambda_scale=wl_scale
               if(nline gt 0) then begin
                 if(line(iline).xfirst ge 0 and                   $
                    line(iline).xlast  lt n_elements(cs)) then    $
                   mark,line(iline).xfirst,line(iline).xlast,     $
                        cs(line(iline).xfirst:line(iline).xlast), $
                        csscale,lambda_scale=wl_scale,            $
                        round(line(iline).posm),line(iline).wll,col=RED
               endif
             end
     'POWER':begin
               polynom_power=event.index+1
               is_computed=0
               widget_control,b7_0,SENSITIVE=0
             end
      'COMP':if(n_elements(where(line(0:nline-1).flag eq 0 and $
                                 line(0:nline-1).wll gt 0)) gt polynom_power) then begin
               ii=where(line(0:nline-1).flag eq 0 and $
                        line(0:nline-1).wll gt 0, nii)
               disp_poly=poly_fit(line(ii).posm,line(ii).wll,polynom_power,/DOUBLE)
               disp_poly=reform(disp_poly)
               line.wlc=poly(reform(line.posm),disp_poly)
               rms=sqrt(total((line(ii).wll-line(ii).wlc)^2)/(nii*(nii-1)))
               log=[log(0)]
               log=[log,string(polynom_power,'(''Polynom power='',I1)')]
               log=[log,string(disp_poly,'(3E21.12)')]
               for i=0,nline-1 do begin
                 posc=fz_roots([disp_poly(0)-line(i).wll,disp_poly(1:polynom_power)],/DOUBLE)
                 line(i).posc=posc((where(min(abs(posc-line(i).posm)) eq $
                                              abs(posc-line(i).posm)))(0))
                 log=[log,string(i,line(i).posm,line(i).wll,line(i).wlc,line(i).flag, $
                 '(''k='',I2,'', Meas.Pos.='',F8.2,'', Lab.WL='',F10.4,'', Comp.WL='',F10.4,'', Flag='',I1)')]
               endfor
               log=[log,string(rms,'(''Standard deviation='',f11.8)')]
               widget_control,b7_8,set_value=string(rms,'(f11.6,'' A'')')
               if(n_elements(wl_scale) gt 1) then wl_scale=poly(dindgen(n_elements(cs)),disp_poly)
               draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
                   ,thar_flag=Show_ThAr_flag,thar_list=thar_list $
                   ,lambda_scale=wl_scale
               mark,line(iline).xfirst,line(iline).xlast,     $
                    cs(line(iline).xfirst:line(iline).xlast), $
                    csscale,lambda_scale=wl_scale,            $
                    round(line(iline).posm),line(iline).wll,col=RED
               if(is_computed eq 0) then begin
                 is_computed=1
                 widget_control,b7_0,SENSITIVE=1
                 widget_control,b2_3,SENSITIVE=1
               endif
               dsp_fit,cs,line,nline,iline,disp_poly,event
               draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
                   ,thar_flag=Show_ThAr_flag,thar_list=thar_list $
                   ,lambda_scale=wl_scale
               if(n_elements(wl_scale) eq 1) then $
                 xy=convert_coord(line(iline).posm,0,/DATA,/TO_DEVICE) $
               else $
                 xy=convert_coord(line(iline).wll,0,/DATA,/TO_DEVICE)
               widget_control,XPlotdisplay,get_draw_view=xview
               widget_control,XPlotdisplay $
                             ,set_draw_view=[(0>(xy(0)-xviewport*0.5))< $
                                                (xfullsize-xviewport-20),0]
               if(line(iline).xfirst ge 0 and                   $
                  line(iline).xlast  lt n_elements(cs)) then    $
                 mark,line(iline).xfirst,line(iline).xlast,     $
                      cs(line(iline).xfirst:line(iline).xlast), $
                      csscale,lambda_scale=wl_scale,            $
                      round(line(iline).posm),line(iline).wll,col=RED
             endif else bbb=widget_message('Not enough lines to fit the polynomial', $
                                           TITLE='Disp',/ERROR)
    'APPGAU':if(event.select eq 1 and iline ge 0) then begin
               xxx=line(iline).xfirst+dindgen(line(iline).xlast-line(iline).xfirst+1)
               yapp=disp_gaussfit(xxx,cs(line(iline).xfirst:line(iline).xlast),a,xxx)
               line(iline).height=abs(a(0))
               line(iline).posm=a(1)
               line(iline).width=abs(a(2))
               line(iline).approx='G'
               widget_control,b4_5,set_value=string(line(iline).width,'(f7.4)')
             endif else if(event.select eq 1) then default_approx='G'
    'APPCEN':if(event.select eq 1 and iline ge 0) then begin
               xxx=line(iline).xfirst+dindgen(line(iline).xlast-line(iline).xfirst+1)
               yapp=disp_gaussboxfit(xxx,cs(line(iline).xfirst:line(iline).xlast),a,xxx)
               line(iline).height=abs(a(0))
               line(iline).posm=a(1)
               line(iline).width=abs(a(2))
               line(iline).approx='C'
               widget_control,b4_5,set_value=string(line(iline).width,'(f7.4)')
             endif else if(event.select eq 1) then default_approx='C'
    'APPPAR':if(event.select eq 1 and iline ge 0) then begin
               a=poly_fit(line(iline).xfirst+dindgen(line(iline).xlast-line(iline).xfirst+1),$
                          cs(line(iline).xfirst:line(iline).xlast),2)
               yapp=poly(line(iline).xfirst+dindgen(line(iline).xlast-line(iline).xfirst+1),a)
               line(iline).height=poly(line(iline).posm-line(iline).xfirst,a)-a(0)
               line(iline).width=(line(iline).xlast-line(iline).xfirst)*0.5
               line(iline).posm=-a(1)/(2*a(2))
               line(iline).approx='P'
               widget_control,b4_5,set_value=' 0.0000'
             endif else if(event.select eq 1) then default_approx='P'
    'APPLOR':if(event.select eq 1 and iline ge 0) then begin
               xxx=line(iline).xfirst+dindgen(line(iline).xlast-line(iline).xfirst+1)
               yapp=disp_lorentzfit(xxx,cs(line(iline).xfirst:line(iline).xlast),a,xxx)
               line(iline).height=a(0)/(a(2)*a(2))
               line(iline).posm=a(1)
               line(iline).width=abs(a(2))
               line(iline).approx='L'
               widget_control,b4_5,set_value=string(line(iline).width,'(f7.4)')
             endif else if(event.select eq 1) then default_approx='L'
     'WLLAB':if(iline ge 0) then begin
               widget_control,b4_1,get_value=www
               line(iline).wll=double(www(0))
               draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
                   ,thar_flag=Show_ThAr_flag,thar_list=thar_list $
                   ,lambda_scale=wl_scale
             end
    'AUTOID':if(is_computed) then begin
               w1=poly(0.d0,disp_poly) & w2=poly(double(n_elements(cs)),disp_poly)
               ww1=w1<w2 & ww2=w1>w2
               ii=where(wllist gt ww1-1 and wllist lt ww2+1, nlist)
               ncs=n_elements(cs)
               jmax=where(cs(1:ncs-2) gt cs(0:ncs-3) and $   ; Locate main local maxima
                          cs(1:ncs-2) gt cs(2:ncs-1), njmax)
               jjj=(reverse(sort(cs(jmax+1))))(0:(njmax<nlist)-1)
               jmax=jmax(jjj)+1
               njmax=n_elements(jmax)
               jmin=where(cs(1:ncs-2) lt cs(0:ncs-3) and $   ; Locate all local minima
                          cs(1:ncs-2) lt cs(2:ncs-1), njmin)+1
               for j=0,njmax-1 do begin
                 if(jmax(j) gt 10 and jmax(j) lt ncs-10) then begin
                   xx1=max(where(jmax(j) gt jmin))
                   xx2=min(where(jmax(j) lt jmin))
                   xxx=dindgen(jmin(xx2)-jmin(xx1)+1)+jmin(xx1)
                                       ; Check if there is an overlap with earlier identification
                   if(nline gt 0) then ll=where(jmax(j) ge line(0:nline-1).xfirst and $
                                                jmax(j) le line(0:nline-1).xlast, nll) $
                   else                nll=0
                   if(n_elements(xxx) ge 8 and nll eq 0) then begin
                     yapp=disp_gaussfit(xxx,cs(xxx),a,xxx)
                     posm=a(1)
                     wlc=poly(posm,disp_poly)
                     if(posm gt 0 and posm lt ncs-1) then begin
                       if(min(abs(line.posm-posm)) gt 0.5 and $
                          min(abs(wllist(ii)-wlc),i) lt 0.02) then begin
                         p=fz_roots([disp_poly(0)-wllist(ii(i)), $
                                     disp_poly(1:polynom_power)],/DOUBLE)
                         posc=p((where(p ge 0 and p lt ncs-1))(0))
                         line(nline).xfirst=min(xxx)
                         line(nline).xlast=max(xxx)
                         line(nline).posm=posm
                         line(nline).posc=posc
                         line(nline).wll =wllist(ii(i))
                         line(nline).wlc =wlc
                         line(nline).height=abs(a(0))
                         line(nline).width=abs(a(2))
                         line(nline).approx='G'
                         line(nline).flag=1
                         line=[line,{csline,-1.d0,-1.d0,-1.d0,-1.d0,0,0,'G',0.d0,0,0.d0,0}]
                         nline=nline+1
                       endif
                     endif
                   endif
                 endif
               endfor
               ii=[sort(line(0:nline-1).posm),nline]
               line=line(ii)
               iline=(where(ii eq iline))(0)
               draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
                   ,thar_flag=Show_ThAr_flag,thar_list=thar_list $
                   ,lambda_scale=wl_scale
               mark,line(iline).xfirst,line(iline).xlast,     $
                    cs(line(iline).xfirst:line(iline).xlast), $
                    csscale,lambda_scale=wl_scale,            $
                    round(line(iline).posm),line(iline).wll,col=RED
               widget_control,b4_5,set_value=string(line(iline).width,'(f7.4)')
               widget_control,b4_6,set_value=string(iline,'(I7)')
             endif
     'PRINT':draw,cs,csscale,csscale_thar,line,nline,log(0),thar_spec=thar_spec $
                 ,thar_flag=Show_ThAr_flag,thar_list=thar_list     $
                 ,/PRINT,PLOT=default_device,lambda_scale=wl_scale
      'REDO':if(nline gt 0) then begin
               offset=0.d0
               for i=0,nline-1 do begin
;                 xf=(line(i).xfirst-3)>0
;                 xl=(line(i).xlast+3)<(n_elements(cs)-1)
;                 dummy=max(cs(xf:xl),ix0)
;                 ix0=ix0+xf
;                 xfl=indgen(ix0-xf+1)+xf
;                 ii=where(cs(xfl+1) gt cs(xfl) and cs(xfl-1) gt cs(xfl), nii)
;                 if(nii gt 0) then line(i).xfirst=xfl(max(ii))
;                 xfl=indgen(xl-ix0-1)+ix0+1
;                 ii=where(cs(xfl+1) gt cs(xfl) and cs(xfl-1) gt cs(xfl), nii)
;                 if(nii gt 0) then line(i).xlast=xfl(min(ii))
                 line(i).xfirst=line(i).xfirst>0
                 line(i).xlast =line(i).xlast<(n_elements(cs)-1)
                 case line(i).approx of
                   'G':begin
                         xxx=line(i).xfirst+dindgen(line(i).xlast- $
                                line(i).xfirst+1)
                         if(line(i).xfirst ge 0 and $
                           line(i).xlast lt n_elements(cs)) then begin
                           yapp=disp_gaussfit(xxx, $
                                cs(line(i).xfirst:line(i).xlast),a,xxx)
                           offset=offset+(line(i).posm-a(1))
                           line(i).height=abs(a(0))
                           line(i).posm=a(1)
                           line(i).width=abs(a(2))
                         endif
                       end
                   'C':begin
                         xxx=line(i).xfirst+dindgen(line(i).xlast- $
                                line(i).xfirst+1)
                         if(line(i).xfirst ge 0 and $
                           line(i).xlast lt n_elements(cs)) then begin
                           yapp=disp_gaussboxfit(xxx, $
                                cs(line(i).xfirst:line(i).xlast),a,xxx)
                           offset=offset+(line(i).posm-a(1))
                           line(i).height=abs(a(0))
                           line(i).posm=a(1)
                           line(i).width=abs(a(2))
                         endif
                       end
                   'P':begin
                         a=poly_fit(dindgen(line(i).xlast-line(i).xfirst+1), $
                                    cs(line(i).xfirst:line(i).xlast),2)
                         offset=offset+(line(i).posm+a(1)/(2*a(2))-line(i).xfirst)
                         line(i).posm=-a(1)/(2*a(2))+line(i).xfirst
                         line(i).height=poly(line(i).posm-line(i).xfirst,a)-a(0)
                         line(i).width=(line(i).xlast-line(i).xfirst)*0.5
                       end
                   'L':begin
                         xxx=line(i).xfirst+dindgen(line(i).xlast- $
                                line(i).xfirst+1)
                         if(line(i).xfirst ge 0 and $
                            line(i).xlast lt n_elements(cs)) then begin
                           yapp=disp_lorentzfit(xxx,cs(line(i).xfirst:line(i).xlast),a,xxx)
                           offset=offset+(line(i).posm-a(1))
                           line(i).height=a(0)/(a(2)*a(2))
                           line(i).posm=a(1)
                           line(i).width=abs(a(2))
                         endif
                       end
                 endcase
               endfor
               offset=offset/nline
;               line.xfirst=line.xfirst-offset
;               line.xlast =line.xlast -offset
               draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
                   ,thar_flag=Show_ThAr_flag,thar_list=thar_list $
                   ,lambda_scale=wl_scale
               if(line(iline).xfirst ge 0 and                   $
                  line(iline).xlast  lt n_elements(cs)) then    $
                 mark,line(iline).xfirst,line(iline).xlast,     $
                      cs(line(iline).xfirst:line(iline).xlast), $
                      csscale,lambda_scale=wl_scale,            $
                      round(line(iline).posm),line(iline).wll,col=RED
               widget_control,b4_5,set_value=string(line(iline).width,'(f7.4)')
               widget_control,b4_6,set_value=string(iline,'(I7)')
             endif
     'SHIFT':begin
               WIDGET_CONTROL,sldr,get_value=shft
               pos_shift=double(shft(0))/410.D0*0.2d0*(n_elements(cs))
               if(nline gt 0) then begin
                 line(0:nline-1).posm  =line(0:nline-1).posm  +pos_shift-pos_shift0
                 line(0:nline-1).xfirst=line(0:nline-1).xfirst+pos_shift-pos_shift0
                 line(0:nline-1).xlast =line(0:nline-1).xlast +pos_shift-pos_shift0
                 pos_shift0=pos_shift
               endif
               draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
                   ,thar_flag=Show_ThAr_flag,thar_list=thar_list $
                   ,lambda_scale=wl_scale
               if(nline gt 0) then begin
                 if(line(iline).xfirst ge 0 and                   $
                    line(iline).xlast  lt n_elements(cs)) then    $
                   mark,line(iline).xfirst,line(iline).xlast,     $
                        cs(line(iline).xfirst:line(iline).xlast), $
                        csscale,lambda_scale=wl_scale,            $
                        round(line(iline).posm),line(iline).wll,col=RED
                 widget_control,b4_5,set_value=string(line(iline).width,'(f7.4)')
                 widget_control,b4_6,set_value=string(iline,'(I7)')
               endif
             end
     'SCALE':begin
               WIDGET_CONTROL,scly,get_value=newscl
               csscale(0)=csscale0(0)*(newscl-1)/99.+(100-newscl)/99.
               csscale(1)=csscale0(1)*(newscl-1)/99.
               draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
                   ,thar_flag=Show_ThAr_flag,thar_list=thar_list $
                   ,lambda_scale=wl_scale
               if(nline gt 0) then begin
                 if(line(iline).xfirst ge 0 and                   $
                    line(iline).xlast  lt n_elements(cs)) then    $
                   mark,line(iline).xfirst,line(iline).xlast,     $
                        cs(line(iline).xfirst:line(iline).xlast), $
                        csscale,lambda_scale=wl_scale,            $
                        round(line(iline).posm),line(iline).wll,col=RED
               endif
             end
 'SCALEThAr':begin
               WIDGET_CONTROL,sclthar,get_value=newscl
               csscale_thar=99*(newscl-1)/99.+1.
               draw,cs,csscale,csscale_thar,line,nline,thar_spec=thar_spec $
                   ,thar_flag=Show_ThAr_flag,thar_list=thar_list $
                   ,lambda_scale=wl_scale
               if(nline gt 0) then begin
                 if(line(iline).xfirst ge 0 and                   $
                    line(iline).xlast  lt n_elements(cs)) then    $
                   mark,line(iline).xfirst,line(iline).xlast,     $
                        cs(line(iline).xfirst:line(iline).xlast), $
                        csscale,lambda_scale=wl_scale,            $
                        round(line(iline).posm),line(iline).wll,col=RED
               endif
             end
        else:widget_control,disp_curve,/CLEAR_EVENTS
    endcase
  endwhile
  widget_control,disp_curve,/DESTROY
  loadct,0
  !p.position=0
  !p.background=0
  return,disp_poly
end

; Identify comparison spectrum lines and
; fit the dispersion curve with a polynomial

function disp,cs,                        $ ; Comparison spectrum
              line_id=lines,             $ ; Line identification structure
              log=log,                   $ ; Fitting log file
              wl_range=wlrange,          $ ; Wavelength region dblarr(2)
              poly_power=poly_power,     $ ; Starting polynom power
              widget_title=widget_title, $ ; Widget name (if available)
              GROUP_LEADER=wGroup,       $ ; Parent widget if needed
              file=filename                ; Filename prefix

  if(not keyword_set(cs)) then begin
    print,'Usage: b=disp(cs[,line_id=lines[,log=log[,wl_range=wlrange[,poly_power=power[,file=file]]]]])
    print,'where: cs      - 1D vector with comparison spectrum (input, required),'
    print,'       lines   - line identification structure (input/output, optional)'
    print,'       log     - string array with log of comparison spectrum fit (output, optional)'
    print,'       wlrange - 2 element array with the approximate wavelength range used'
    print,'                 only to select the subset of lines from the line list (input, recomended)'
    print,'       power   - set this to starting polynom power. Defaults is 2.
    print,'       file    - filename prefix used for printout, log ID etc. (input, optional)'
    return,-1
  endif
;
;  Construct the compound widget
;
  if((!d.flags AND 65536) ne 65536) then begin
    print,'The graphics device is not capable to display widgets'
    return,[0,1]
  endif

  thar=kpno_thar(wrange=wlrange,thar=thar_spec,maxline=80)
  if(thar[0] gt 0) then begin
    wllist =reform(thar(*,0))
    thar_list=thar
    thar=0
  endif
  if(not keyword_set(thar_list)) then return,0

  if(not keyword_set(filename)) then filename='disp.fit' ; Default filename

  if(n_elements(lines) gt 0) then begin   ; We have an initial solution
    if(n_tags(lines[0]) gt 0) then begin  ; It is a structure
      if(tag_names(lines[0],/STRUCTURE_NAME) eq 'CSLINE') then begin
        j=where(lines.wll gt 0., nj)
        if(nj gt 0) then begin
          line=lines[j]
          line=line(sort(line.posm))
        endif
      endif
    endif
  endif

  dsp_defaults,font,nfont,default_device,box_clr
  log=strcompress('Comparison spectrum: '+filename)
  disp_poly=dsp_widget(cs,wllist,font,nfont,default_device,box_clr, $
                       line_id=line,log=log,poly_power=poly_power, $
                       thar_spec=thar_spec,thar_list=thar_list, $
                       widget_title=widget_title,GROUP_LEADER=wGroup)
  if(n_elements(line) gt 1) then lines=line[0:n_elements(line)-2]
  return,disp_poly
end

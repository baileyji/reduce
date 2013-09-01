;=========================================================================================
;Modification history:
;
;  19/02/2008 (NP) Written
;=========================================================================================
pro refit,line,cs_spec,error=error
  nlines=n_elements(line)
  if(arg_present(error)) then error=replicate(1.d0,nlines)
  for i=0L,nlines-1 do begin
    if((size(cs_spec))(0) gt 1) then cs=cs_spec[*,line[i].order] $
    else                             cs=cs_spec
    case line[i].approx of
    'G':begin
          xxx=line[i].xfirst+dindgen(line[i].xlast-line[i].xfirst+1)
          if(line[i].xfirst ge 0 and $
             line[i].xlast lt n_elements(cs)) then begin
            yapp=disp_gaussfit(xxx,cs[line[i].xfirst:line[i].xlast],a,xxx)
            line[i].height=abs(a[0])
            line[i].posm=a[1]
            line[i].width=abs(a[2])
            if(arg_present(error)) then $
              error[i]=total(abs(yapp-cs[xxx]))/total(cs[xxx])+max(a[3]-cs[xxx]/a[0])
          endif
        end
    'C':begin
          xxx=line[i].xfirst+dindgen(line[i].xlast-line[i].xfirst+1)
          if(line[i].xfirst ge 0 and $
             line[i].xlast lt n_elements(cs)) then begin
            yapp=disp_gaussboxfit(xxx,cs[line[i].xfirst:line[i].xlast],a,xxx)
            line[i].height=abs(a[0])
            line[i].posm=a[1]
            line[i].width=abs(a[2])
            if(arg_present(error)) then $
              error[i]=total(abs(yapp-cs[xxx]))/total(cs[xxx])+max(a[3]-cs[xxx]/a[0])
          endif
        end
    'P':begin
          xxx=line[i].xfirst+dindgen(line[i].xlast-line[i].xfirst+1)
          a=poly_fit(dindgen(line[i].xlast-line[i].xfirst+1), $
                     cs[line[i].xfirst:line[i].xlast],2,/DOUBLE)
          line[i].posm=-a[1]/(2*a[2])+line[i].xfirst
          line[i].height=poly(line[i].posm-line[i].xfirst,a)-a[0]
          line[i].width=(line[i].xlast-line[i].xfirst)*0.5
          if(arg_present(error)) then $
            error[i]=total(abs(poly(xxx-line[i].xfirst,a)-cs[xxx]))/total(cs[xxx]) $
                    +max(a[0]-cs[xxx]/line[i].height)
        end
    'L':begin
          xxx=line[i].xfirst+dindgen(line[i].xlast-line[i].xfirst+1)
          if(line[i].xfirst ge 0 and $
             line[i].xlast lt n_elements(cs)) then begin
            yapp=disp_lorentzfit(xxx,cs[line[i].xfirst:line[i].xlast],a,xxx)
            line[i].height=a[0]/(a[2]*a[2])
            line[i].posm=a[1]
            line[i].width=abs(a[2])
            if(arg_present(error)) then $
              error[i]=total(abs(yapp-cs[xxx]))/total(cs[xxx])+max(a[3]-cs[xxx]/a[0])
          endif
        end
    endcase
  endfor
  return
end

;=========================================================================================
;Modification history:
;
;  09/01/2008 (NP) Added a generation of graphics log file (PostScript) with the line statistics,
;                  resolving power, plots of PSF and PSF variations accross the detector
;=========================================================================================
Pro make_log,cs_spec,cs_lines,oincr,obase,solution_2D,file_prefix $
            ,GROUP_LEADER=wGroup
;
; Dimensions
  npixels=n_elements(cs_spec[*,0])
  norders=n_elements(cs_spec[0,*])
  nlines =n_elements(cs_lines)
;
; Average line shape
  i = where(cs_lines.approx eq 'G', ngauss)
  i = where(cs_lines.approx eq 'P', npara)
  i = where(cs_lines.approx eq 'C', ngbox)
  i = where(cs_lines.approx eq 'L', nloren)
;
; Variables to hold the resolution, dispersion, quadrant limits and polynom power
  resolution = dblarr(nlines)
  dsp        = dblarr(nlines)
  res = dblarr(3,3)
  xlim = round(dindgen(4)/3*npixels)
  ylim = round(dindgen(4)/3*norders)
  polynom_power = round(solution_2D[8])
;
; Loop though spectral orders and compute resolution at each ThAr line
  for iord=0,norders-1 do begin
    ii=where(cs_lines.order eq iord, nii)
    if(nii gt 0) then begin
      abs_order = oincr*iord+obase                       ; Absolute order number
      mkwave, w, solution_2D, abs_order                  ; Produce wavelength scale for selected order
      disp_poly = poly_fit(dindgen(npixels), w $         ; Reconstruct 1D solution for order iord
                         , polynom_power, /DOUBLE)
      disp_poly = disp_poly*dindgen(polynom_power+1)
      dsp[ii] = poly(cs_lines[ii].posm, disp_poly[1:*])  ;local dispersion in Angstroems per pixel
      resolution[ii] = cs_lines[ii].wll/abs(cs_lines[ii].width*dsp[ii]) ;Resolution of individual lines
    endif
  endfor
;
; Loop though spectral orders and compute resolution at each ThAr line
  xx=0.d0
  yy=0.d0
  i00=0L
  i01=0L
  i02=0L
  i10=0L
  i11=0L
  i12=0L
  i20=0L
  i21=0L
  i22=0L
  xx=0.d0
  yy=0.d0
;
; Find dispersion at the center
  x=min((npixels*0.5d0-cs_lines.posm)^2+100*(norders*0.5d0-cs_lines.order)^2,i)
  dsp0=dsp[i]
  wl0 =cs_lines[i].wll
  nn=0L

  for iord=0,norders-1 do begin
    ii=where(cs_lines.order eq iord, nii)
    if(nii gt 0) then begin
      for i=0,nii-1 do begin
        j=ii[i]
        n=cs_lines[j].xlast-cs_lines[j].xfirst+1
        x=cs_lines[j].xfirst+dindgen(n)
        y=cs_spec[x,iord]
        if(cs_lines[j].posm  ge xlim[0] and cs_lines[j].posm  lt xlim[1] and $
           cs_lines[j].order ge ylim[0] and cs_lines[j].order lt ylim[1]) then begin
          i00=[i00,nn+lindgen(n)]
          nn=nn+n
        endif else if(cs_lines[j].posm  ge xlim[1] and cs_lines[j].posm  lt xlim[2] and $
                      cs_lines[j].order ge ylim[0] and cs_lines[j].order lt ylim[1]) then begin
          i10=[i10,nn+lindgen(n)]
          nn=nn+n
        endif else if(cs_lines[j].posm  ge xlim[2] and cs_lines[j].posm  lt xlim[3] and $
           cs_lines[j].order ge ylim[0] and cs_lines[j].order lt ylim[1]) then begin
          i20=[i20,nn+lindgen(n)]
          nn=nn+n
        endif else if(cs_lines[j].posm  ge xlim[0] and cs_lines[j].posm  lt xlim[1] and $
           cs_lines[j].order ge ylim[1] and cs_lines[j].order lt ylim[2]) then begin
          i01=[i01,nn+lindgen(n)]
          nn=nn+n
        endif else if(cs_lines[j].posm  ge xlim[1] and cs_lines[j].posm  lt xlim[2] and $
           cs_lines[j].order ge ylim[1] and cs_lines[j].order lt ylim[2]) then begin
          i11=[i11,nn+lindgen(n)]
          nn=nn+n
        endif else if(cs_lines[j].posm  ge xlim[2] and cs_lines[j].posm  lt xlim[3] and $
           cs_lines[j].order ge ylim[1] and cs_lines[j].order lt ylim[2]) then begin
          i21=[i21,nn+lindgen(n)]
          nn=nn+n
        endif else if(cs_lines[j].posm  ge xlim[0] and cs_lines[j].posm  lt xlim[1] and $
           cs_lines[j].order ge ylim[2] and cs_lines[j].order lt ylim[3]) then begin
          i02=[i02,nn+lindgen(n)]
          nn=nn+n
        endif else if(cs_lines[j].posm  ge xlim[1] and cs_lines[j].posm  lt xlim[2] and $
           cs_lines[j].order ge ylim[2] and cs_lines[j].order lt ylim[3]) then begin
          i12=[i12,nn+lindgen(n)]
          nn=nn+n
        endif else if(cs_lines[j].posm  ge xlim[2] and cs_lines[j].posm  lt xlim[3] and $
           cs_lines[j].order ge ylim[2] and cs_lines[j].order lt ylim[3]) then begin
          i22=[i22,nn+lindgen(n)]
          nn=nn+n
        endif
        yapp=disp_gaussfit(x-cs_lines[j].posm,y,a,x-cs_lines[j].posm)
        
;jjj=where((y-a[3])/cs_lines[j].height lt -0.2, njjj)
;if(njjj gt 0) then stop
;if(max((y-a[3])/a[0]) gt 2.) then stop
        yy=[yy,(y-a[3])/cs_lines[j].height]
        xx=[xx,(x-cs_lines[j].posm)]
      endfor
    endif
  endfor
  xx=xx[1:*]
  yy=yy[1:*]
  if(n_elements(i00) gt 1) then i00=i00[1:*]
  if(n_elements(i10) gt 1) then i10=i10[1:*]
  if(n_elements(i20) gt 1) then i20=i20[1:*]
  if(n_elements(i01) gt 1) then i01=i01[1:*]
  if(n_elements(i11) gt 1) then i11=i11[1:*]
  if(n_elements(i21) gt 1) then i21=i21[1:*]
  if(n_elements(i02) gt 1) then i02=i02[1:*]
  if(n_elements(i12) gt 1) then i12=i12[1:*]
  if(n_elements(i22) gt 1) then i22=i22[1:*]
  napp=n_elements(xx)>1000L
  xapp=(abs(max(xx))>abs(min(xx)))*(dindgen(napp)/(napp-1.d0)-0.5d0)
  yapp=disp_gaussfit(xx,yy,a,xapp)
;
; Construct PS file name
  i=strpos(file_prefix,'.fit')
  if(i gt 0) then file=strmid(file_prefix, 0, i)+'.report.ps' $
  else            file=file_prefix+'.report.ps'

  old_device=!d.name
  set_plot,'ps'
  device,file=file,/COLOR,/TIMES,/BOLD,/ISOLATIN,xs=19,xoff=1,ys=23,yoff=3
  !p.font=0
  colors
  plot,xapp,yapp,xs=3,ys=1,yr=minmax(yy),/NODATA,tit=file_prefix
  oplot,xx,yy,psym=3
  oplot,xapp,yapp,col=4,thick=3
  xyouts,!x.crange[0]+1,0.9,'# of lines '+strtrim(nlines,2),charsiz=1
  xyouts,!x.crange[0]+1,0.8,'Median resolving power ='+string(median(resolution),'(F7.0)'),charsiz=1
  xyouts,!x.crange[0]+1,0.7,'Median dispersion ['+string(197B)+'/pixel] ='+string(median(dsp),'(f7.4)'),charsiz=1
  !p.position=0 & !p.region=0
  !p.multi=[0,3,3]
  erase
  for iy=2,0,-1 do begin
    for ix=0,2 do begin
      ii=where(cs_lines.posm  ge xlim[ix] and cs_lines.posm  lt xlim[ix+1] $
           and cs_lines.order ge ylim[iy] and cs_lines.order lt ylim[iy+1], nii)
      if(nii gt 0) then begin
        if(     ix eq 0 and iy eq 0) then iii=i00 $
        else if(ix eq 1 and iy eq 0) then iii=i10 $
        else if(ix eq 2 and iy eq 0) then iii=i20 $
        else if(ix eq 0 and iy eq 1) then iii=i01 $
        else if(ix eq 1 and iy eq 1) then iii=i11 $
        else if(ix eq 2 and iy eq 1) then iii=i21 $
        else if(ix eq 0 and iy eq 2) then iii=i02 $
        else if(ix eq 1 and iy eq 2) then iii=i12 $
        else if(ix eq 2 and iy eq 2) then iii=i22
        yapp=disp_gaussfit(xx[iii],yy[iii],a,xapp)
        plot,xapp,yapp,xs=3,ys=1,yr=minmax(yy),/NODATA,tit=(ix eq 1 and iy eq 1)?file_prefix:' '
        oplot,xx[iii],yy[iii],psym=3
        oplot,xapp,yapp,col=4,thick=3
        xyouts,!x.crange[0]+1,1.0,'# of lines '+strtrim(nii,2),charsiz=0.4
        xyouts,!x.crange[0]+1,0.9,'Resolving power '+string(median(resolution[ii]),'(F7.0)'),charsiz=0.4
        xyouts,!x.crange[0]+1,0.8,'Dispersion ['+string(197B)+'/pxl]'+string(median(dsp[ii]),'(f7.4)'),charsiz=0.4
        xyouts,!x.crange[1]-1,1.0,'Pixels: '+strtrim(xlim[ix],2)+'-'+strtrim(xlim[ix+1],2),charsiz=0.4,align=1
        xyouts,!x.crange[1]-1,0.9,'Orders: '+strtrim(ylim[iy],2)+'-'+strtrim(ylim[iy+1],2),charsiz=0.4,align=1
      endif else begin
        plot,[0,1],/NODATA,xs=-1,ys=-1
        xyouts,0.5,0.5,'No ThAr lines identified in this quandrant!',align=0.5
      endelse
    endfor
  endfor
  device,/close
  !p.font=-1
  !p.multi=0
  set_plot,old_device
  loadct,0,/SILENT
;
  return
end

;=========================================================================================
;Modification history:
;
; 02-Oct-2003 NP Wrote 2D wavecal based on disp2d algorithm by J. Valenti
; 15-Feb-2004 NP Fixed gauss box fit to the line profiles: removed all debugging
;                plot and switched to using the MFIT minimization package
; 18-Mar-2004 NP Fixed obase change in case the alignment procedure shifts the
;                orders. obase change depends on oincr while cs_lines.order
;                which is a relative order number should change in the opposit
;                sense from oincr*delta_obase. This simply means that no matter
;                what obase and oincr are relative order numbers should start
;                from 0
; 23-Mar-2004 NP Fixed contrast sliders behavior so that the upper cut never
;                gets to the level of the lower cut
; 17-Oct-2005 NP Added automatic clipping of the outlayers with rebuilding of
;                the 2D solution
; 18-Feb-2008 NP Switch to use refit subroutin for re-dermining line centers while
;                keeping the approximation
;=========================================================================================
Function disp_align,cs_spec,cs_lines,obase,oincr,GROUP_LEADER=wGroup
  common align_display_size, base_size

  npixels=n_elements(cs_spec[*,0])
  norders=n_elements(cs_spec[0,*])

;
; Check if our solution is compatible with data
  jj=where(cs_lines.xlast lt npixels, njj)
  if(njj gt 0) then begin
    cs_lines=cs_lines[jj]
  endif else begin
    loadct,0,/SILENT
    return,[0,0]
  endelse

  if(n_elements(base_size) ne 2) then begin
    device, get_screen_size = scr_sz
    base_size=[(scr_sz[0]<1280),scr_sz[1]]
  endif

  BASE_0 = Widget_Base(GROUP_LEADER=wGroup $
      , XOFFSET=5, YOFFSET=5 $
      , /TLB_SIZE_EVENTS, UVALUE='RESIZE' $
      , /COLUMN, TITLE='Model alignment', SPACE=3, XPAD=3 ,YPAD=3)

  BASE_1 = Widget_Base(BASE_0, UVALUE='BASE_1', /ROW)

  BUTTON_CANCEL = Widget_Button(BASE_1, UVALUE='CANCEL'  $
                  , /ALIGN_CENTER, TOOLTIP='Discards alignment results', VALUE='DISCARD')

  BUTTON_EXIT   = Widget_Button(BASE_1, UVALUE='EXIT'  $
                  , /ALIGN_CENTER, TOOLTIP='Feed alignment results to the main program' $
                  , VALUE='EXIT')
  BASE_1a = Widget_Base(BASE_1, UVALUE='BASE_1a', /COLUMN)

  cs_model=cs_spec*0.

  for iord=0,norders-1 do begin
    ii=where(cs_lines.order eq iord, nii)
    if(nii gt 0) then begin
      for i=0,nii-1 do begin
        j=ii[i]
        x=cs_lines[j].xfirst+dindgen(cs_lines[j].xlast-cs_lines[j].xfirst+1)
        x=(x-cs_lines[j].posm)/(cs_lines[j].width*0.5)
        cs_model[cs_lines[j].xfirst:cs_lines[j].xlast,iord]= $
           cs_model[cs_lines[j].xfirst:cs_lines[j].xlast,iord]+exp(-x*x)*cs_lines[j].height
      endfor
    endif
  endfor

  cs_data=cs_spec
  for iord=0,norders-1 do cs_data[*,iord]=cs_data[*,iord]-median(cs_data[*,iord])
  cs_data=cs_data>0.

  csdata =hist_equal(cs_data ,binsize=(max(cs_data )-min(cs_data ))/127.,top=127)
  csmodel=hist_equal(cs_model,binsize=(max(cs_model)-min(cs_model))/127.,top=127)
  r=where(csmodel gt 0B)
  csmodel[r]=csmodel[r]+127B

  max_value=max(round(cs_data))>max(round(cs_model))
  min_value=min(round(cs_data))<min(round(cs_model))
  max_value=max_value*0.25
  max_min_value=max_value
  min_max_value=min_value

  MAX_VAL = Widget_Slider(BASE_1a, MIN=min_value+(max_value-min_value)/799. $
            , MAX=max_value, SCROLL=1, UVALUE='MAX_VAL' $
            , SCR_XSIZE=base_size[0]-200 $
            , VALUE=max_value)
  MIN_VAL = Widget_Slider(BASE_1a, MIN=min_value $
            , MAX=max_value-(max_value-min_value)/799. $
            , SCROLL=1, UVALUE='MIN_VAL' $
            , SCR_XSIZE=base_size[0]-200 $
            , VALUE=min_value)

;=========== Alignment sliders and plot window =================================
  X_SLIDER = Widget_Slider(BASE_0, MIN=-npixels/2, MAX=npixels/2, SCROLL=1 $
            , UVALUE='X_SHIFT' $
            , SCR_XSIZE=base_size[0]-200 $
            , VALUE=0)

  X_SCALE  = Widget_Slider(BASE_0, MIN=-100, MAX=100, SCROLL=1 $
            , UVALUE='X_SCALE' $
            , SCR_XSIZE=base_size[0]-200 $
            , VALUE=0)
  xscale=1.d0

  BASE_2 = Widget_Base(BASE_0, UVALUE='BASE_2', /ROW)

  DRAW = Widget_Draw(BASE_2, UVALUE='ZOOM', /SCROLL $
        , XSIZE=npixels, YSIZE=(norders*20L)>(base_size[1]-180) $
        , X_SCROLL_SIZE=base_size[0]-160 $
        , Y_SCROLL_SIZE=base_size[1]-180  $
        , RETAIN=2, /BUTTON_EVENTS)

  Y_SLIDER = Widget_Slider(BASE_2, MIN=-norders/2,MAX=norders/2, SCROLL=1 $
            , UVALUE='Y_SHIFT' $
            , SCR_YSIZE=base_size[1]-150 $
            , VALUE=0, /VERTICAL)

  xscl=(base_size[1]-160.)/((norders*20.)>(base_size[1]-100.))
  yscl=((base_size[0]<1280)-120.)/double(npixels)

  widget_control, BASE_0, /REALIZE ; Realize the widget
  widget_control, BASE_0, TLB_GET_SIZE=base_size_old

  r=[bindgen(128)*2B,bytarr(127)    ,255B]
  g=[bytarr(128)    ,bytarr(127)    ,255B]
  b=[bytarr(128)    ,bindgen(127)*2B,255B]
  tvlct,r,g,b

  index=indgen(norders)*2

  im=bytarr(npixels,2*norders)
  im[*,index  ]=csmodel
  im[*,index+1]=csdata
  !x.ticklen=0.01*xscl
  !y.ticklen=0.01*yscl

  display,im
  for i=0,norders do oplot,!x.crange,2*[i,i]-0.5
  x_shift=0
  y_shift=0
  mark_ON=0

; Workaround for erratical resize events coming through network
; in client-server use
  resize=0 

  while(1) do begin
    if(resize eq 1) then begin
      repeat begin
        event=widget_event(BASE_0, BAD_ID=bad, /NOWAIT)
      endrep until(event.id eq 0L)
      event = ev_resize
    endif else begin
      event=widget_event(BASE_0, BAD_ID=bad)
    endelse
    if(bad ne 0) then begin
      !x.ticklen=0
      !y.ticklen=0
      loadct,0,/SILENT
      return,[0,0]
    endif
;    if(event.id eq event.top) then begin ; Resizeing the drawing window and the sliders
;      Widget_Control, DRAW, DRAW_XSIZE = (event.x- 40)>800 $
;                          , DRAW_YSIZE = (event.y-140)>600
;      display,im
;     for i=0,norders do oplot,!x.crange,2*[i,i]-0.5
;      if(mark_ON eq 2) then begin
;        oplot,x_model+[-2,2,2,-2,-2],y_model+[0.5,0.5,-0.5,-0.5,0.5]
;      endif else if(mark_ON eq 1) then begin
;        oplot,x_data +[-2,2,2,-2,-2],y_data +[0.5,0.5,-0.5,-0.5,0.5]
;      endif
;      continue
;    endif
    widget_control,event.id,GET_UVALUE=userid
    case userid of
   'RESIZE':if(resize eq 0) then begin
              resize = 1
              ev_resize = event
            endif else if(resize eq 1) then begin
              resize=0
              widget_control, BASE_0,TLB_GET_SIZE=base_size,UPDATE=0
              if(base_size[0] ne base_size_old[0] or $
                 base_size[1] ne base_size_old[1]) then begin
                base_size[0]=base_size[0]>800
                base_size[1]=base_size[1]>600
                widget_control, Y_SLIDER, /DESTROY
                widget_control, DRAW, /DESTROY
                widget_control, BASE_2, /DESTROY
                widget_control, X_SCALE, /DESTROY
                widget_control, X_SLIDER, /DESTROY
                widget_control, MIN_VAL, /DESTROY
                widget_control, MAX_VAL, /DESTROY
                widget_control, BASE_1a, /DESTROY
                widget_control, BUTTON_EXIT, /DESTROY
                widget_control, BUTTON_CANCEL, /DESTROY
                widget_control, BASE_1, /DESTROY
                BASE_1 = Widget_Base(BASE_0, UVALUE='BASE_1', /ROW)
                BUTTON_CANCEL = Widget_Button(BASE_1, UVALUE='CANCEL'  $
                              , /ALIGN_CENTER, TOOLTIP='Discards alignment results', VALUE='DISCARD')
                BUTTON_EXIT   = Widget_Button(BASE_1, UVALUE='EXIT'  $
                              , /ALIGN_CENTER, TOOLTIP='Feed alignment results to the main program' $
                              , VALUE='EXIT')
                BASE_1a = Widget_Base(BASE_1, UVALUE='BASE_1a', /COLUMN)
                MAX_VAL = Widget_Slider(BASE_1a, MIN=min_value+(max_value-min_value)/799. $
                              , MAX=max_value, SCROLL=1, UVALUE='MAX_VAL' $
                              , SCR_XSIZE=base_size[0]-200 $
                              , VALUE=max_value)
                MIN_VAL = Widget_Slider(BASE_1a, MIN=min_value $
                              , MAX=max_value-(max_value-min_value)/799. $
                              , SCROLL=1, UVALUE='MIN_VAL' $
                              , SCR_XSIZE=base_size[0]-200 $
                              , VALUE=min_value)
                X_SLIDER = Widget_Slider(BASE_0, MIN=-npixels/2, MAX=npixels/2, SCROLL=1 $
                              , UVALUE='X_SHIFT' $
                              , SCR_XSIZE=base_size[0]-200 $
                              , VALUE=x_shift)
                X_SCALE  = Widget_Slider(BASE_0, MIN=-100, MAX=100, SCROLL=1 $
                              , UVALUE='X_SCALE' $
                              , SCR_XSIZE=base_size[0]-200 $
                              , VALUE=(xscale-1.)*2.d4)
                BASE_2 = Widget_Base(BASE_0, UVALUE='BASE_2', /ROW)
                DRAW = Widget_Draw(BASE_2, UVALUE='ZOOM', /SCROLL $
                              , XSIZE=npixels, YSIZE=(norders*20L)>(base_size[1]-180) $
                              , X_SCROLL_SIZE=base_size[0]-160 $
                              , Y_SCROLL_SIZE=base_size[1]-180  $
                              , RETAIN=2, /BUTTON_EVENTS)
                Y_SLIDER = Widget_Slider(BASE_2, MIN=-norders/2,MAX=norders/2, SCROLL=1 $
                              , UVALUE='Y_SHIFT' $
                              , SCR_YSIZE=base_size[1]-150 $
                              , VALUE=y_shift, /VERTICAL)
                widget_control, BASE_0, TLB_GET_SIZE=base_size_old, /UPDATE
                display,im
                for i=0,norders do oplot,!x.crange,2*[i,i]-0.5
                if(mark_ON eq 2) then begin
                  oplot,x_model+[-2,2,2,-2,-2],y_model+[0.5,0.5,-0.5,-0.5,0.5]
                endif else if(mark_ON eq 1) then begin
                  oplot,x_data +[-2,2,2,-2,-2],y_data +[0.5,0.5,-0.5,-0.5,0.5]
                endif
              endif
            endif
     'ZOOM':begin
              WIDGET_CONTROL,event.id,get_value=a  ; Get cursor coordinates
              xy=convert_coord(event.x,event.y,/TO_DATA,/DEVICE)
              xy=[round(xy[0]),round(xy[1])]
              if(xy[0] gt !x.crange[0]+2 and xy[0] lt !x.crange[1]-2 and $
                 xy[1] gt !y.crange[0]   and xy[1] lt !y.crange[1]) then begin
                dummy=max(im[xy[0]-2:xy[0]+2,xy[1]],i)
                xy[0]=i+xy[0]-2
                if(mark_ON eq 0 and event.type eq 0) then begin ; Now button is pressed
                  oplot,xy[0]+[-2,2,2,-2,-2],xy[1]+[0.5,0.5,-0.5,-0.5,0.5]
                  if(xy[1] mod 2 eq 0) then begin
                    x_model=xy[0]
                    y_model=xy[1]/2
                    mark_ON=2
                  endif else begin
                    x_data=xy[0]
                    y_data=xy[1]/2
                    mark_ON=1
                  endelse
                endif else if(mark_ON ne 0 and event.type eq 0) then begin ; Now button is pressed again
                  if((xy[1] mod 2 eq 0 and mark_ON eq 1) or $
                     (xy[1] mod 2 eq 1 and mark_ON eq 2)) then begin
                    mark_ON=0                         ; Clear flag
                    oplot,xy[0]+[-2,2,2,-2,-2],xy[1]+[0.5,0.5,-0.5,-0.5,0.5]
                    if(xy[1] mod 2 eq 0) then begin
                      x_model=xy[0]
                      y_model=xy[1]/2
                    endif else begin
                      x_data=xy[0]
                      y_data=xy[1]/2
                    endelse
                    wait,0.5
                    x_shift=x_shift-x_model+x_data
                    y_shift=y_shift-y_model+y_data
                    im=bytarr(npixels,2*norders)
                    im[*,index+1]=csdata
                    if(x_shift gt 0 and y_shift gt 0) then begin
                      im[x_shift:*,index[y_shift:*]]=csmodel[0:npixels-x_shift-1,0:norders-y_shift-1]
                    endif else if(x_shift lt 0 and y_shift gt 0) then begin
                      im[0:npixels+x_shift-1,index[y_shift:*]]=csmodel[-x_shift:*,0:norders-y_shift-1]
                    endif else if(x_shift gt 0 and y_shift lt 0) then begin
                      im[x_shift:*,index[0:norders+y_shift-1]]=csmodel[0:npixels-x_shift-1,-y_shift:*]
                    endif else if(x_shift lt 0 and y_shift lt 0) then begin
                      im[0:npixels+x_shift-1,index[0:norders+y_shift-1]]=csmodel[-x_shift:*,-y_shift:*]
                    endif else if(x_shift gt 0 and y_shift eq 0) then begin
                      im[x_shift:*,index]=csmodel[0:npixels-x_shift-1,*]
                    endif else if(x_shift lt 0 and y_shift eq 0) then begin
                      im[0:npixels+x_shift-1,index]=csmodel[-x_shift:*,*]
                    endif else if(x_shift eq 0 and y_shift gt 0) then begin
                      im[*,index[y_shift:*]]=csmodel[*,0:norders-y_shift-1]
                    endif else if(x_shift eq 0 and y_shift lt 0) then begin
                      im[*,index[0:norders+y_shift-1]]=csmodel[*,-y_shift:*]
                    endif else im[*,index]=csmodel
                    display,im
                    for i=0,norders do oplot,!x.crange,2*[i,i]-0.5
                    widget_control,X_SLIDER,set_value=x_shift
                    widget_control,Y_SLIDER,set_value=y_shift
                  endif
                endif
              endif
            end
  'MAX_VAL':begin
              widget_control,MAX_VAL,get_value=r
              max_min_value=floor(r[0])
;              print,'max)',min_max_value,max_min_value
              csdata =hist_equal(cs_data,min=min_max_value $
                                ,binsize=(max_min_value-min_max_value)/127. $
                                ,top=127)
              csmodel=hist_equal(cs_model,min=min_max_value $
                                ,binsize=(max_min_value-min_max_value)/127. $
                                ,top=127)
              r=where(csmodel gt 0B)
              csmodel[r]=csmodel[r]+127B
              im=bytarr(npixels,2*norders)
              im[*,index+1]=csdata
              if(x_shift gt 0 and y_shift gt 0) then begin
                im[x_shift:*,index[y_shift:*]]=csmodel[0:npixels-x_shift-1,0:norders-y_shift-1]
              endif else if(x_shift lt 0 and y_shift gt 0) then begin
                im[0:npixels+x_shift-1,index[y_shift:*]]=csmodel[-x_shift:*,0:norders-y_shift-1]
              endif else if(x_shift gt 0 and y_shift lt 0) then begin
                im[x_shift:*,index[0:norders+y_shift-1]]=csmodel[0:npixels-x_shift-1,-y_shift:*]
              endif else if(x_shift lt 0 and y_shift lt 0) then begin
                im[0:npixels+x_shift-1,index[0:norders+y_shift-1]]=csmodel[-x_shift:*,-y_shift:*]
              endif else if(x_shift gt 0 and y_shift eq 0) then begin
                im[x_shift:*,index]=csmodel[0:npixels-x_shift-1,*]
              endif else if(x_shift lt 0 and y_shift eq 0) then begin
                im[0:npixels+x_shift-1,index]=csmodel[-x_shift:*,*]
              endif else if(x_shift eq 0 and y_shift gt 0) then begin
                im[*,index[y_shift:*]]=csmodel[*,0:norders-y_shift-1]
              endif else if(x_shift eq 0 and y_shift lt 0) then begin
                im[*,index[0:norders+y_shift-1]]=csmodel[*,-y_shift:*]
              endif else im[*,index]=csmodel
              widget_control,MIN_VAL,set_slider_max=max_min_value- $
                                     (max_min_value-min_max_value)/799.
              widget_control,MIN_VAL,set_value=min_max_value
              display,im
              for i=0,norders do oplot,!x.crange,2*[i,i]-0.5
              if(mark_ON eq 2) then begin
                oplot,x_model+[-2,2,2,-2,-2],y_model+[0.5,0.5,-0.5,-0.5,0.5]
              endif else if(mark_ON eq 1) then begin
                oplot,x_data +[-2,2,2,-2,-2],y_data +[0.5,0.5,-0.5,-0.5,0.5]
              endif
            end
  'MIN_VAL':begin
              widget_control,MIN_VAL,get_value=r
              min_max_value=ceil(r[0])
;              print,'min)',min_max_value,max_min_value
              csdata =hist_equal(cs_data,min=min_max_value $
                                ,binsize=(max_min_value-min_max_value)/127. $
                                ,top=127)
              csmodel=hist_equal(cs_model,min=min_max_value $
                                ,binsize=(max_min_value-min_max_value)/127. $
                                ,top=127)
              r=where(csmodel gt 0B)
              csmodel[r]=csmodel[r]+127B
              im=bytarr(npixels,2*norders)
              im[*,index+1]=csdata
              if(x_shift gt 0 and y_shift gt 0) then begin
                im[x_shift:*,index[y_shift:*]]=csmodel[0:npixels-x_shift-1,0:norders-y_shift-1]
              endif else if(x_shift lt 0 and y_shift gt 0) then begin
                im[0:npixels+x_shift-1,index[y_shift:*]]=csmodel[-x_shift:*,0:norders-y_shift-1]
              endif else if(x_shift gt 0 and y_shift lt 0) then begin
                im[x_shift:*,index[0:norders+y_shift-1]]=csmodel[0:npixels-x_shift-1,-y_shift:*]
              endif else if(x_shift lt 0 and y_shift lt 0) then begin
                im[0:npixels+x_shift-1,index[0:norders+y_shift-1]]=csmodel[-x_shift:*,-y_shift:*]
              endif else if(x_shift gt 0 and y_shift eq 0) then begin
                im[x_shift:*,index]=csmodel[0:npixels-x_shift-1,*]
              endif else if(x_shift lt 0 and y_shift eq 0) then begin
                im[0:npixels+x_shift-1,index]=csmodel[-x_shift:*,*]
              endif else if(x_shift eq 0 and y_shift gt 0) then begin
                im[*,index[y_shift:*]]=csmodel[*,0:norders-y_shift-1]
              endif else if(x_shift eq 0 and y_shift lt 0) then begin
                im[*,index[0:norders+y_shift-1]]=csmodel[*,-y_shift:*]
              endif else im[*,index]=csmodel
              widget_control,MAX_VAL,set_slider_min=min_max_value+ $
                                     (max_min_value-min_max_value)/799.
              widget_control,MAX_VAL,set_value=max_min_value
              display,im
              for i=0,norders do oplot,!x.crange,2*[i,i]-0.5
              if(mark_ON eq 2) then begin
                oplot,x_model+[-2,2,2,-2,-2],y_model+[0.5,0.5,-0.5,-0.5,0.5]
              endif else if(mark_ON eq 1) then begin
                oplot,x_data +[-2,2,2,-2,-2],y_data +[0.5,0.5,-0.5,-0.5,0.5]
              endif
            end
  'X_SHIFT':begin
              widget_control,X_SLIDER,get_value=r
              r=round(r[0])
              im=bytarr(npixels,2*norders)
              im[*,index+1]=csdata
              x_shift=r
              if(x_shift gt 0 and y_shift gt 0) then begin
                im[x_shift:*,index[y_shift:*]]=csmodel[0:npixels-x_shift-1,0:norders-y_shift-1]
              endif else if(x_shift lt 0 and y_shift gt 0) then begin
                im[0:npixels+x_shift-1,index[y_shift:*]]=csmodel[-x_shift:*,0:norders-y_shift-1]
              endif else if(x_shift gt 0 and y_shift lt 0) then begin
                im[x_shift:*,index[0:norders+y_shift-1]]=csmodel[0:npixels-x_shift-1,-y_shift:*]
              endif else if(x_shift lt 0 and y_shift lt 0) then begin
                im[0:npixels+x_shift-1,index[0:norders+y_shift-1]]=csmodel[-x_shift:*,-y_shift:*]
              endif else if(x_shift gt 0 and y_shift eq 0) then begin
                im[x_shift:*,index]=csmodel[0:npixels-x_shift-1,*]
              endif else if(x_shift lt 0 and y_shift eq 0) then begin
                im[0:npixels+x_shift-1,index]=csmodel[-x_shift:*,*]
              endif else if(x_shift eq 0 and y_shift gt 0) then begin
                im[*,index[y_shift:*]]=csmodel[*,0:norders-y_shift-1]
              endif else if(x_shift eq 0 and y_shift lt 0) then begin
                im[*,index[0:norders+y_shift-1]]=csmodel[*,-y_shift:*]
              endif else im[*,index]=csmodel
              display,im
              for i=0,norders do oplot,!x.crange,2*[i,i]-0.5
              if(mark_ON eq 2) then begin
                oplot,x_model+[-2,2,2,-2,-2],y_model+[0.5,0.5,-0.5,-0.5,0.5]
              endif else if(mark_ON eq 1) then begin
                oplot,x_data +[-2,2,2,-2,-2],y_data +[0.5,0.5,-0.5,-0.5,0.5]
              endif
            end
  'X_SCALE':begin
              widget_control,X_SCALE,get_value=r
              im=bytarr(npixels,2*norders)
              im[*,index+1]=csdata
              xscale=r[0]*5.d-5+1.d0
              cs_model=cs_spec*0.

              for iord=0,norders-1 do begin
                ii=where(cs_lines.order eq iord, nii)
                if(nii gt 0) then begin
                  for i=0,nii-1 do begin
                    j=ii[i]
                    x=cs_lines[j].xfirst+dindgen(cs_lines[j].xlast-cs_lines[j].xfirst+1)
                    x=(x-cs_lines[j].posm*xscale)/(cs_lines[j].width*0.5)
                    cs_model[cs_lines[j].xfirst:cs_lines[j].xlast,iord]= $
                    cs_model[cs_lines[j].xfirst:cs_lines[j].xlast,iord]+exp(-x*x)*cs_lines[j].height
                  endfor
                endif
              endfor
              csmodel=hist_equal(cs_model,min=min_max_value $
                                ,binsize=(max_min_value-min_max_value)/127. $
                                ,top=127)
              r=where(csmodel gt 0B)
              csmodel[r]=csmodel[r]+127B

              if(x_shift gt 0 and y_shift gt 0) then begin
                im[x_shift:*,index[y_shift:*]]=csmodel[0:npixels-x_shift-1,0:norders-y_shift-1]
              endif else if(x_shift lt 0 and y_shift gt 0) then begin
                im[0:npixels+x_shift-1,index[y_shift:*]]=csmodel[-x_shift:*,0:norders-y_shift-1]
              endif else if(x_shift gt 0 and y_shift lt 0) then begin
                im[x_shift:*,index[0:norders+y_shift-1]]=csmodel[0:npixels-x_shift-1,-y_shift:*]
              endif else if(x_shift lt 0 and y_shift lt 0) then begin
                im[0:npixels+x_shift-1,index[0:norders+y_shift-1]]=csmodel[-x_shift:*,-y_shift:*]
              endif else if(x_shift gt 0 and y_shift eq 0) then begin
                im[x_shift:*,index]=csmodel[0:npixels-x_shift-1,*]
              endif else if(x_shift lt 0 and y_shift eq 0) then begin
                im[0:npixels+x_shift-1,index]=csmodel[-x_shift:*,*]
              endif else if(x_shift eq 0 and y_shift gt 0) then begin
                im[*,index[y_shift:*]]=csmodel[*,0:norders-y_shift-1]
              endif else if(x_shift eq 0 and y_shift lt 0) then begin
                im[*,index[0:norders+y_shift-1]]=csmodel[*,-y_shift:*]
              endif else im[*,index]=csmodel
              display,im
              for i=0,norders do oplot,!x.crange,2*[i,i]-0.5
              if(mark_ON eq 2) then begin
                oplot,x_model+[-2,2,2,-2,-2],y_model+[0.5,0.5,-0.5,-0.5,0.5]
              endif else if(mark_ON eq 1) then begin
                oplot,x_data +[-2,2,2,-2,-2],y_data +[0.5,0.5,-0.5,-0.5,0.5]
              endif
            end
  'Y_SHIFT':begin
              widget_control,Y_SLIDER,get_value=r
              r=round(r[0])
              im=bytarr(npixels,2*norders)
              im[*,index+1]=csdata
              y_shift=r
              if(x_shift gt 0 and y_shift gt 0) then begin
                im[x_shift:*,index[y_shift:*]]=csmodel[0:npixels-x_shift-1,0:norders-y_shift-1]
              endif else if(x_shift lt 0 and y_shift gt 0) then begin
                im[0:npixels+x_shift-1,index[y_shift:*]]=csmodel[-x_shift:*,0:norders-y_shift-1]
              endif else if(x_shift gt 0 and y_shift lt 0) then begin
                im[x_shift:*,index[0:norders+y_shift-1]]=csmodel[0:npixels-x_shift-1,-y_shift:*]
              endif else if(x_shift lt 0 and y_shift lt 0) then begin
                im[0:npixels+x_shift-1,index[0:norders+y_shift-1]]=csmodel[-x_shift:*,-y_shift:*]
              endif else if(x_shift gt 0 and y_shift eq 0) then begin
                im[x_shift:*,index]=csmodel[0:npixels-x_shift-1,*]
              endif else if(x_shift lt 0 and y_shift eq 0) then begin
                im[0:npixels+x_shift-1,index]=csmodel[-x_shift:*,*]
              endif else if(x_shift eq 0 and y_shift gt 0) then begin
                im[*,index[y_shift:*]]=csmodel[*,0:norders-y_shift-1]
              endif else if(x_shift eq 0 and y_shift lt 0) then begin
                im[*,index[0:norders+y_shift-1]]=csmodel[*,-y_shift:*]
              endif else im[*,index]=csmodel
              display,im
              for i=0,norders do oplot,!x.crange,2*[i,i]-0.5
              if(mark_ON eq 2) then begin
                oplot,x_model+[-2,2,2,-2,-2],y_model+[0.5,0.5,-0.5,-0.5,0.5]
              endif else if(mark_ON eq 1) then begin
                oplot,x_data +[-2,2,2,-2,-2],y_data +[0.5,0.5,-0.5,-0.5,0.5]
              endif
            end
     'EXIT':begin
              widget_control, BASE_0, /DESTROY
              !x.ticklen=0
              !y.ticklen=0
              loadct,0,/SILENT
              widget_control,/HOURGLASS

              csmodel=cs_lines

              i=where(cs_lines.order+y_shift lt norders and $ ; All orders
                      cs_lines.order+y_shift ge 0, nlines)    ; fitting new data
              if(nlines gt 0) then csmodel=csmodel[i] else return,[0,0]
              csmodel.order=csmodel.order+y_shift
              obase=obase-oincr*y_shift

              i=where(cs_lines.xlast*xscale +x_shift lt npixels and $ ; All pixels
                      cs_lines.xfirst*xscale+x_shift ge 0, nlines)   ; fitting new data
              if(nlines gt 0) then csmodel=csmodel[i] else return,[0,0]
              csmodel.xfirst=csmodel.xfirst*xscale+x_shift
              csmodel.posm  =csmodel.posm*xscale  +x_shift
              csmodel.xlast =csmodel.xlast*xscale +x_shift
              refit,csmodel,cs_spec;[xxx,csmodel[i].order]
;              for i=0,nlines-1 do begin ; Re-determine line centers
;                xxx=dindgen(csmodel[i].xlast-csmodel[i].xfirst+1)+csmodel[i].xfirst
;                yapp=disp_gaussfit(xxx,reform(cs_spec[xxx,csmodel[i].order]),a,xxx)  ; Do Gaussian fit
;                posm=a[1]                                ; Measured position (Center of Gaussian)
;                csmodel[i].posc=posm                     ; Store temporary
;                csmodel[i].height=abs(a[0])
;                csmodel[i].width=abs(a[2])
;                csmodel[i].approx='G'
;              endfor
              i=where(csmodel.height gt 0 and $
                      csmodel.width gt 1. and csmodel.width lt 10. and $
                      csmodel.posm gt csmodel.xfirst and $
                      csmodel.posm lt csmodel.xlast, nlines)
              if(nlines eq 0) then return,[0,0]
              csmodel=csmodel[i]

              x_shift=x_shift+mean(csmodel.posc-csmodel.posm)  ; Determine accurate shift
; BUG killed by Eric Stempels 091208. Shift actually goes the other way!
;             csmodel.posm=csmodel.posc                        ; Replace measured positions
              csmodel.posc=csmodel.posm                        ; Replace computed positions

              cs_lines=csmodel                                 ; Replace line list
              return,[x_shift,y_shift]
            end
   'CANCEL':begin
              widget_control, BASE_0, /DESTROY
              !x.ticklen=0
              !y.ticklen=0
              loadct,0,/SILENT
              return,[0,0]
            end
    endcase
  endwhile
  return,0
end

;=========================================================================================
;Function analyze_residuals,m_pix,m_ord,m_flag,npixel,nord,obase,resid_pix,resid_m_s $
;                          ,good=jgood,GROUP_LEADER=wGroup

Function analyze_residuals,m_wave,m_pix,m_ord,m_flag,maxordx $
                          ,maxordy,npixel,nord,obase,oincr,resid_pix,resid_m_s $
                          ,good=jgood,GROUP_LEADER=wGroup,CLIP=res_clip

  device, get_screen_size = scr_sz

  RES_BASE_0 = Widget_Base(GROUP_LEADER=wGroup, UVALUE='RES_BASE_0' $
      , XOFFSET=5, YOFFSET=5, /TLB_SIZE_EVENTS $
;      , SCR_XSIZE=(scr_sz[0]-20)<1024, SCR_YSIZE=(scr_sz[1]-20)<900 $
      , /COLUMN , TITLE='Analysis of the residuals', SPACE=3, XPAD=3 ,YPAD=3)

;=========== Control buttons =================================
  RES_BASE_1 = Widget_Base(RES_BASE_0, UVALUE='RES_BASE_1', /ROW)

  RES_BASE_2 = Widget_Base(RES_BASE_1, UVALUE='RES_BASE_2', /ROW, /EXCLUSIVE)

  RES_BUTTON_0 = Widget_Button(RES_BASE_2, UVALUE='BY_ORDERS'  $
      , /ALIGN_CENTER, TOOLTIP='Show residuals for different spectral orders' $
      , VALUE='By orders')

  RES_BUTTON_1 = Widget_Button(RES_BASE_2, UVALUE='BY_PIXELS'  $
      , /ALIGN_CENTER, TOOLTIP='Show residuals for different column swaths' $
      , VALUE='By pixels')

  RES_BASE_3 = Widget_Base(RES_BASE_1, UVALUE='RES_BASE_3', /ROW, /EXCLUSIVE)

  RES_BUTTON_2 = Widget_Button(RES_BASE_3, UVALUE='RES_PIXEL'  $
      , /ALIGN_CENTER, TOOLTIP='Show pixel residuals' $
      , VALUE='Pixels')

  RES_BUTTON_3 = Widget_Button(RES_BASE_3, UVALUE='RES_WAVE'  $
      , /ALIGN_CENTER, TOOLTIP='Show residuals in m/s' $
      , VALUE='m/s')

  RES_LABEL_1 = Widget_Label(RES_BASE_1, UVALUE='WID_LABEL_7'  $
      , /ALIGN_LEFT, VALUE='Clipping at [m/s]', XOFFSET=10)

; Auto-clipping
  if(not keyword_set(res_clip)) then res_clip=200.
  nclip=5

  RES_TEXT_1 = Widget_Text(RES_BASE_1, UVALUE='RMS_M_S'  $
      , VALUE=[ string(res_clip,form='(f6.0)') ], XSIZE=6, YSIZE=1 $
      , /EDITABLE, /KBRD_FOCUS_EVENTS)

  RES_BUTTON_4 = Widget_Button(RES_BASE_1, UVALUE='AUTO'  $
      , /ALIGN_CENTER, TOOLTIP='Automatically reject outlier', VALUE='AUTO')

  RES_BUTTON_5 = Widget_Button(RES_BASE_1, UVALUE='CANCEL'  $
      , /ALIGN_CENTER, TOOLTIP='Discard outlier rejection', VALUE='DISCARD')

  RES_BUTTON_6 = Widget_Button(RES_BASE_1, UVALUE='EXIT'  $
      , /ALIGN_CENTER, TOOLTIP='Complete outlier rejection', VALUE='EXIT')

;=========== Draw window =================================
  RES_DRAW = Widget_Draw(RES_BASE_0, UVALUE='ZOOM' $
      , SCR_XSIZE=((scr_sz[0]-40)<1024)-10, SCR_YSIZE=(scr_sz[1]-110)<820 $
      , XSIZE=((scr_sz[0]-40)<1024)-10, YSIZE=(scr_sz[1]-110)<820  $
      , RETAIN=1, /BUTTON_EVENTS, /MOTION_EVENTS)

  ngood_on_entry=n_elements(resid_pix)
  jgood=lindgen(ngood_on_entry)

; Column bin information for "by pixels" mode
  nbin=20
  bin=round(npixel/float(nbin))
  bins=[indgen(20)*bin,npixel]

  widget_control, RES_BUTTON_0, /SET_BUTTON
  widget_control, RES_BUTTON_3, /SET_BUTTON
  widget_control, RES_BASE_0, /REALIZE ; Realize the widget

; Residuals for different orders
  devtot=stddev(resid_pix)
  nord=max(m_ord)-min(m_ord)+1
  dev=fltarr(nord)
  med=fltarr(nord)
  res_pixel=0
  resid=resid_m_s
  ytit='Residuals (m/s)'
  for i=0,nord-1 do begin
    ii=where(m_ord eq i*oincr+obase, nii)
    if(nii gt 3) then begin
      dev[i]=stddev(resid[ii])
      med[i]=median(resid[ii])
    endif else begin
      dev[i]=devtot
      med[i]=0.
    endelse
  endfor
  jj=where(m_flag eq 1, njj, complement=ii) ; Use different symbols for flagged lines
  plot,m_ord[ii],resid[ii],psym=4,xs=3,ys=3,yr=minmax(resid) $
      ,xtitle='Spectral orders',ytitle=ytit,xr=minmax(m_ord)
  if(njj gt 0) then oplot,m_ord[jj],resid[jj],psym=1
  oplot,indgen(nord)*oincr+obase, 3*dev,col=200,psym=10
  oplot,indgen(nord)*oincr+obase,-3*dev,col=160,psym=10
  oplot,indgen(nord)*oincr+obase,med,line=3

  zoom_ON=0
  by_orders=1
  box_clr=127

  while(1) do begin
    event=widget_event(RES_BASE_0, BAD_ID=bad)
    if(bad ne 0) then return,0
    if(event.id eq event.top) then begin
      Widget_Control, RES_DRAW, DRAW_XSIZE = (event.x- 40)>800 $
                              , DRAW_YSIZE = (event.y-110)>600
      if(by_orders eq 1) then begin     ; Selected view is by orders
        jj=where(m_flag eq 1, njj, complement=ii) ; Use different symbols for flagges lines
        plot,m_ord[ii],resid[ii],psym=4,xs=3,ys=3,yr=minmax(resid) $
            ,xtitle='Spectral orders',ytitle=ytit
        if(njj gt 0) then oplot,m_ord[jj],resid[jj],psym=1
        oplot,indgen(nord)*oincr+obase, 3*dev,col=160,psym=10
        oplot,indgen(nord)*oincr+obase,-3*dev,col=160,psym=10
        oplot,indgen(nord)*oincr+obase,med,line=3
      endif else begin ; Residuals for different groups of column along orders
        jj=where(m_flag eq 1, njj, complement=ii) ; Use different symbols for flagges lines
        plot,m_pix[ii],resid[ii],psym=4,xs=3,ys=3,yr=minmax(resid) $
            ,xtitle='Columns',ytitle=ytit
        if(njj gt 0) then oplot,m_pix[jj],resid[jj],psym=1
        oplot,indgen(nbin)*bin, 3*dev,col=160,psym=10
        oplot,indgen(nbin)*bin,-3*dev,col=160,psym=10
        oplot,indgen(nbin)*bin,med,line=3
      endelse
      continue
    endif
    widget_control,event.id,GET_UVALUE=userid
    case userid of
     'ZOOM':begin
              WIDGET_CONTROL,event.id,get_value=a  ; Get cursor coordinates
              xy=convert_coord(event.x,event.y,/TO_DATA,/DEVICE)
              if(xy[0] ge !x.crange[0] and xy[0] le !x.crange[1] and $
                 xy[1] ge !y.crange[0] and xy[1] le !y.crange[1]) then begin
                if(zoom_ON eq 0 and event.type eq 0) then begin ; Now button was pressed
                  zoom_ON=1                         ; Set flag
                  device,set_graphics=6             ; Set XOR
                  xx=xy(0) & yy=xy(1)               ; Select one corner
                  box_px=replicate(xx,5)
                  box_py=replicate(yy,5)
                  oplot,box_px,box_py,thick=1, $    ; Plot new box
                        /noclip,lines=0,col=box_clr
                endif else if((zoom_ON eq 1) and $  ; Mouse moved while button
                   (event.type eq 2)) then begin    ; is still pressed
                  oplot,box_px,box_py,thick=1, $    ; Erase old box
                        /noclip,lines=0,col=box_clr
                  box_px(1)=xy(0) & box_px(2)=xy(0) ; Update box coordinates
                  box_py(2)=xy(1) & box_py(3)=xy(1)
                  oplot,box_px,box_py,thick=1, $    ; Plot new box
                        /noclip,lines=0,col=box_clr
                endif else if((zoom_ON eq 1) and $  ; Button was released
                   (event.type eq 1)) then begin
                  zoom_ON=0
                  device,set_graphics=3             ; Return to the normal mode
                  x0=min(box_px) & y0=min(box_py)   ; Determine LL corner
                  x1=max(box_px) & y1=max(box_py)   ; Determine UR corner

                  if(by_orders eq 1) then begin     ; Selected view is by orders
                    ii=where(m_ord gt x0 and m_ord lt x1 and $
                             resid gt y0 and resid lt y1, nii, complement=jj)
                    if(nii gt 0) then begin
                      jgood=jgood[jj]
                      m_pix=m_pix[jj]
                      m_ord=m_ord[jj]
                      m_flag=m_flag[jj]
                      m_wave=m_wave[jj]
                      resid=resid[jj]
                      resid_pix=resid_pix[jj]
                      resid_m_s=resid_m_s[jj]
                    endif
                    dev=fltarr(nord)
                    med=fltarr(nord)
                    for i=0,nord-1 do begin
                      ii=where(m_ord eq i*oincr+obase, nii)
                      if(nii gt 3) then dev[i]=stddev(resid[ii]) else dev[i]=devtot
                      if(nii gt 0) then med[i]=median(resid[ii]) else med[i]=0.
                    endfor
                    jj=where(m_flag eq 1, njj, complement=ii) ; Use different symbols for flagges lines
                    plot,m_ord[ii],resid[ii],psym=4,xs=3,ys=3,yr=minmax(resid) $
                        ,xtitle='Spectral orders',ytitle=ytit
                    if(njj gt 0) then oplot,m_ord[jj],resid[jj],psym=1
                    oplot,indgen(nord)*oincr+obase, 3*dev,col=160,psym=10
                    oplot,indgen(nord)*oincr+obase,-3*dev,col=160,psym=10
                    oplot,indgen(nord)*oincr+obase,med,line=3
                  endif else begin ; Residuals for different groups of column along orders
                    ii=where(m_pix gt x0 and m_pix lt x1 and $
                             resid gt y0 and resid lt y1, nii, complement=jj)
                    if(nii gt 0) then begin
                      jgood=jgood[jj]
                      m_pix=m_pix[jj]
                      m_ord=m_ord[jj]
                      m_flag=m_flag[jj]
                      m_wave=m_wave[jj]
                      resid=resid[jj]
                      resid_pix=resid_pix[jj]
                      resid_m_s=resid_m_s[jj]
                    endif
                    dev=fltarr(nbin)
                    med=fltarr(nbin)
                    for i=0,nbin-1 do begin
                      ii=where(m_pix ge bins[i] and m_pix lt bins[i+1], nii)
                      if(nii gt 3) then dev[i]=stddev(resid[ii]) else dev[i]=devtot
                      if(nii gt 0) then med[i]=median(resid[ii]) else med[i]=0.
                    endfor
                    jj=where(m_flag eq 1, njj, complement=ii) ; Use different symbols for flagges lines
                    plot,m_pix[ii],resid[ii],psym=4,xs=3,ys=3,yr=minmax(resid) $
                        ,xtitle='Columns',ytitle=ytit
                    if(njj gt 0) then oplot,m_pix[jj],resid[jj],psym=1
                    oplot,indgen(nbin)*bin, 3*dev,col=160,psym=10
                    oplot,indgen(nbin)*bin,-3*dev,col=160,psym=10
                    oplot,indgen(nbin)*bin,med,line=3
                  endelse
                endif
              endif
            end
'BY_ORDERS':begin
              by_orders=1
              dev=fltarr(nord)
              med=fltarr(nord)
              for i=0,nord-1 do begin
                ii=where(m_ord eq i*oincr+obase, nii)
                if(nii gt 3) then dev[i]=stddev(resid[ii]) else dev[i]=devtot
                if(nii gt 0) then med[i]=median(resid[ii]) else med[i]=0.
              endfor
              jj=where(m_flag eq 1, njj, complement=ii) ; Use different symbols for flagges lines
              plot,m_ord[ii],resid[ii],psym=4,xs=3,ys=3,yr=minmax(resid) $
                  ,xtitle='Spectral orders',ytitle=ytit
              if(njj gt 0) then oplot,m_ord[jj],resid[jj],psym=1
              oplot,indgen(nord)*oincr+obase, 3*dev,col=160,psym=10
              oplot,indgen(nord)*oincr+obase,-3*dev,col=160,psym=10
              oplot,indgen(nord)*oincr+obase,med,line=3
            end
'BY_PIXELS':begin
              by_orders=0
              dev=fltarr(nbin)
              med=fltarr(nbin)
              for i=0,nbin-1 do begin
                ii=where(m_pix ge bins[i] and m_pix lt bins[i+1], nii)
                if(nii gt 3) then dev[i]=stddev(resid[ii]) else dev[i]=devtot
                if(nii gt 0) then med[i]=median(resid[ii]) else med[i]=0.
              endfor
              jj=where(m_flag eq 1, njj, complement=ii) ; Use different symbols for flagges lines
              plot,m_pix[ii],resid[ii],psym=4,xs=3,ys=3,yr=minmax(resid) $
                  ,xtitle='Columns',ytitle=ytit
              if(njj gt 0) then oplot,m_pix[jj],resid[jj],psym=1
              oplot,indgen(nbin)*bin, 3*dev,col=160,psym=10
              oplot,indgen(nbin)*bin,-3*dev,col=160,psym=10
              oplot,indgen(nbin)*bin,med,line=3
            end
'RES_PIXEL':begin
              res_pixel=1
              resid=resid_pix
              ytit='Residuals (pixels)'
              devtot=stddev(resid_pix)
              if(by_orders eq 1) then begin
                dev=fltarr(nord)
                med=fltarr(nord)
                for i=0,nord-1 do begin
                  ii=where(m_ord eq i*oincr+obase, nii)
                  if(nii gt 3) then dev[i]=stddev(resid[ii]) else dev[i]=devtot
                  if(nii gt 0) then med[i]=median(resid[ii]) else med[i]=0.
                endfor
                jj=where(m_flag eq 1, njj, complement=ii) ; Use different symbols for flagges lines
                plot,m_ord[ii],resid[ii],psym=4,xs=3,ys=3,yr=minmax(resid) $
                    ,xtitle='Spectral orders',ytitle=ytit
                if(njj gt 0) then oplot,m_ord[jj],resid[jj],psym=1
                oplot,indgen(nord)*oincr+obase, 3*dev,col=160,psym=10
                oplot,indgen(nord)*oincr+obase,-3*dev,col=160,psym=10
                oplot,indgen(nord)*oincr+obase,med,line=3
              endif else begin
                dev=fltarr(nbin)
                med=fltarr(nbin)
                for i=0,nbin-1 do begin
                  ii=where(m_pix ge bins[i] and m_pix lt bins[i+1], nii)
                  if(nii gt 3) then dev[i]=stddev(resid[ii]) else dev[i]=devtot
                  if(nii gt 0) then med[i]=median(resid[ii]) else med[i]=0.
                endfor
                jj=where(m_flag eq 1, njj, complement=ii) ; Use different symbols for flagges lines
                plot,m_pix[ii],resid[ii],psym=4,xs=3,ys=3,yr=minmax(resid) $
                    ,xtitle='Columns',ytitle=ytit
                if(njj gt 0) then oplot,m_pix[jj],resid[jj],psym=1
                oplot,indgen(nbin)*bin, 3*dev,col=160,psym=10
                oplot,indgen(nbin)*bin,-3*dev,col=160,psym=10
                oplot,indgen(nbin)*bin,med,line=3
              endelse
            end
 'RES_WAVE':begin
              res_pixel=0
              resid=resid_m_s
              ytit='Residuals (m/s)'
              devtot=stddev(resid_m_s)
              if(by_orders eq 1) then begin
                dev=fltarr(nord)
                med=fltarr(nord)
                for i=0,nord-1 do begin
                  ii=where(m_ord eq i*oincr+obase, nii)
                  if(nii gt 3) then dev[i]=stddev(resid[ii]) else dev[i]=devtot
                  if(nii gt 0) then med[i]=median(resid[ii]) else med[i]=0.
                endfor
                jj=where(m_flag eq 1, njj, complement=ii) ; Use different symbols for flagges lines
                plot,m_ord[ii],resid[ii],psym=4,xs=3,ys=3,yr=minmax(resid) $
                    ,xtitle='Spectral orders',ytitle=ytit
                if(njj gt 0) then oplot,m_ord[jj],resid[jj],psym=1
                oplot,indgen(nord)*oincr+obase, 3*dev,col=160,psym=10
                oplot,indgen(nord)*oincr+obase,-3*dev,col=160,psym=10
                oplot,indgen(nord)*oincr+obase,med,line=3
              endif else begin
                dev=fltarr(nbin)
                med=fltarr(nbin)
                for i=0,nbin-1 do begin
                  ii=where(m_pix ge bins[i] and m_pix lt bins[i+1], nii)
                  if(nii gt 3) then dev[i]=stddev(resid[ii]) else dev[i]=devtot
                  if(nii gt 0) then med[i]=median(resid[ii]) else med[i]=0.
                endfor
                jj=where(m_flag eq 1, njj, complement=ii) ; Use different symbols for flagges lines
                plot,m_pix[ii],resid[ii],psym=4,xs=3,ys=3,yr=minmax(resid) $
                    ,xtitle='Columns',ytitle=ytit
                if(njj gt 0) then oplot,m_pix[jj],resid[jj],psym=1
                oplot,indgen(nbin)*bin, 3*dev,col=160,psym=10
                oplot,indgen(nbin)*bin,-3*dev,col=160,psym=10
                oplot,indgen(nbin)*bin,med,line=3
              endelse
            end
  'RMS_M_S':begin
              widget_control, RES_TEXT_1, get_value=dummy
              dummy=float(dummy[0])
              if(dummy gt 0. and dummy lt 1000.) then begin
                res_clip=dummy
              endif else begin
                widget_control, RES_TEXT_1 $
                              , set_value=[string(res_clip,form='(f5.0)')]
              endelse
            end
     'AUTO':begin
              nord=round(max(m_ord)-min(m_ord)+1)         ; Range of orders for which we have solutions
              obase_tmp=obase<(obase+oincr*(nord-1))
              if(nord gt 2) then disp_2d,m_wave,m_pix,m_ord,m_flag,maxordx $ ; 2D solver
                                        ,maxordy,npixel,nord,obase_tmp,wvc $
                                        ,resid_pix=resid_pix,resid_m_s=resid_m_s
              resid_pix=resid_pix
              resid_m_s=resid_m_s
              if(res_pixel) then resid=resid_pix else resid=resid_m_s
              repeat begin
                jj=lindgen(n_elements(resid_m_s))
                clip='CLIP_OFF'
                for iclip=1,nclip do begin
                  ii=where(jj ge 0, nii)
                  if(max(abs(resid_m_s[ii]), imax) gt res_clip) then begin
                    clip='CLIP_ON'
                    jj[ii[imax]]=-1
                  endif
                endfor
                ii=where(jj ge 0, nii)
                if(nii eq 0) then clip='CLIP_OFF' $
                else jj=jj[ii]
                if(clip eq 'CLIP_ON') then begin
                  jgood=jgood[jj]
                  m_pix=m_pix[jj]
                  m_ord=m_ord[jj]
                  m_flag=m_flag[jj]
                  m_wave=m_wave[jj]
                  nord=round(max(m_ord)-min(m_ord)+1)         ; Range of orders for which we have solutions
                  obase_tmp=obase<(obase+oincr*(nord-1))
                  if(nord gt 2) then disp_2d,m_wave,m_pix,m_ord,m_flag,maxordx $ ; 2D solver
                                            ,maxordy,npixel,nord,obase_tmp,wvc $
                                            ,resid_pix=resid_pix,resid_m_s=resid_m_s
                  resid_pix=resid_pix
                  resid_m_s=resid_m_s
                  if(res_pixel) then resid=resid_pix else resid=resid_m_s

                  if(by_orders eq 1) then begin     ; Selected view is by orders
                    dev=fltarr(nord)
                    med=fltarr(nord)
                    for i=0,nord-1 do begin
                      ii=where(m_ord eq i*oincr+obase, nii)
                      if(nii gt 3) then dev[i]=stddev(resid[ii]) else dev[i]=devtot
                      if(nii gt 0) then med[i]=median(resid[ii]) else med[i]=0.
                    endfor
                    jj=where(m_flag eq 1, njj, complement=ii) ; Use different symbols for flagges lines
                    plot,m_ord[ii],resid[ii],psym=4,xs=3,ys=3,yr=minmax(resid) $
                        ,xtitle='Spectral orders',ytitle=ytit
                    if(njj gt 0) then oplot,m_ord[jj],resid[jj],psym=1
                    oplot,indgen(nord)*oincr+obase, 3*dev,col=160,psym=10
                    oplot,indgen(nord)*oincr+obase,-3*dev,col=160,psym=10
                    oplot,indgen(nord)*oincr+obase,med,line=3
                    empty
                  endif else begin ; Residuals for groups of column along orders
                    dev=fltarr(nbin)
                    med=fltarr(nbin)
                    for i=0,nbin-1 do begin
                      ii=where(m_pix ge bins[i] and m_pix lt bins[i+1], nii)
                      if(nii gt 3) then dev[i]=stddev(resid[ii]) else dev[i]=devtot
                      if(nii gt 0) then med[i]=median(resid[ii]) else med[i]=0.
                    endfor
                    jj=where(m_flag eq 1, njj, complement=ii) ; Use different symbols for flagges lines
                    plot,m_pix[ii],resid[ii],psym=4,xs=3,ys=3,yr=minmax(resid) $
                        ,xtitle='Columns',ytitle=ytit
                    if(njj gt 0) then oplot,m_pix[jj],resid[jj],psym=1
                    oplot,indgen(nbin)*bin, 3*dev,col=160,psym=10
                    oplot,indgen(nbin)*bin,-3*dev,col=160,psym=10
                    oplot,indgen(nbin)*bin,med,line=3
                    empty
                  endelse
                endif
              endrep until(clip eq 'CLIP_OFF')
            end
     'EXIT':begin
              widget_control, RES_BASE_0, /DESTROY
              return,n_elements(jgood)
            end
   'CANCEL':begin
              widget_control, RES_BASE_0, /DESTROY
              jgood=lindgen(ngood_on_entry)
              return,0
            end
    endcase
  endwhile
  return,0

;  SDtoFWHM=2.d0*sqrt(2.d0*alog(2.d0))
;  !p.region=0 & !p.position=0
;  !p.multi=0
;  device,decomposed=0 & colors
;  iter=0
;  xmin=-0.5>min(resid)
;  xmax= 0.5<max(resid)
;next:
;
;disp_2d,m_wave,m_pix,m_ord,m_flag,maxordx,maxordy,npix, $
;        nord,obase,wvc,residual=resid,wldev=wldev
;
;  dx=(xmax-xmin)*0.01d0
;  i=where(m_flag eq 0)
;  h=histogram(resid(i),bin=dx,min=xmin,max=xmax)
;  x=xmin+dindgen(n_elements(h))*dx
;  xx=rebin(x,n_elements(x)*10)
;  yy=gaussfit(x,h,a,xx)
;  plot,x,h,psym=10,xr=[-1,1]*0.5
;  oplot,xx,yy,col=2
;  i=where(abs(resid) gt a(2)*4./SDtoFWHM, ni)
;  if(ni gt 0) then m_flag(i)=1
;  xmin=-6.*a(2)/SDtoFWHM>min(resid)
;  xmax= 6.*a(2)/SDtoFWHM<max(resid)
;  iter=iter+1
;  if(iter lt 4) then goto,next

end

;=========================================================================================
; This function opens a small separate dialog asking to enter/confirm absolute
; order number and wavelength range
Function setWLrange, rel_order, abs_order, wlrange_in, GROUP_LEADER=wGroup

; Main base
  base=widget_base(title='Order #'+strtrim(rel_order,2), /COLUMN, /FRAME, GROUP_LEADER=wGroup)

; Base for buttons
  bbase=widget_base(base, /ROW, GROUP_LEADER=wGroup)
  b1=widget_button(bbase, value='CANCEL', uvalue='CANCEL', GROUP_LEADER=wGroup)
  b2=widget_button(bbase, value='Continue', uvalue='EXIT', GROUP_LEADER=wGroup)

; Base for input windows
  wlrange=wlrange_in
  tbase0=widget_base(base, /ROW, GROUP_LEADER=wGroup)
  WID_LABEL_0 = Widget_Label(tbase0, UVALUE='WID_LABEL_0'  $
      , /ALIGN_LEFT, VALUE='Absolute order number', XOFFSET=10)
  if(not keyword_set(abs_order)) then abs_order= 80
  WID_SLIDE = Widget_Slider(tbase0, MINIMUM=1,MAXIMUM=500 $
      , YSIZE=48, SCROLL=1, VALUE=abs_order, UVALUE='ABS_ORD', /VERTICAL)
  tbase1=widget_base(base, /ROW, GROUP_LEADER=wGroup)
  WID_LABEL_1 = Widget_Label(tbase1, UNAME='WLmin', /ALIGN_LEFT, VALUE='Minimum wavelength' $
                                  , XOFFSET=10)
  WID_TEXT_1 = Widget_Text(tbase1, UVALUE='WLMIN', VALUE=string(wlrange[0],form='(F11.2)') $
                                , XSIZE=11, YSIZE=1, /EDITABLE, /KBRD_FOCUS_EVENTS)
  tbase2=widget_base(base, /ROW, GROUP_LEADER=wGroup)
  WID_LABEL_2 = Widget_Label(tbase2, UVALUE='WLmax', /ALIGN_LEFT, VALUE='Maximum wavelength' $
                                  , XOFFSET=10)
  WID_TEXT_2 = Widget_Text(tbase2, UVALUE='WLMAX', VALUE=string(wlrange[1],form='(F11.2)') $
                                , XSIZE=11, YSIZE=1, /EDITABLE, /KBRD_FOCUS_EVENTS)

  widget_control, base,/REALIZE ; Realize the widget

  while(1) do begin
    event=widget_event(base, BAD_ID=bad) & if(bad ne 0) then return, [-1, -1]
    widget_control,event.id,GET_UVALUE=userid
    if(Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_KBRD_FOCUS') then userid='REFRESH'
    case userid of
   'REFRESH':if(event.enter eq 0) then begin
               widget_control, WID_TEXT_1, get_value=dummy
               wlrange[0]=dummy[0]
               widget_control, WID_TEXT_2, get_value=dummy
               wlrange[1]=dummy[0]
             endif
   'ABS_ORD':begin
               widget_control, WID_SLIDE, get_value=dummy
               abs_order=dummy[0]
             end
     'WLMIN':begin
               widget_control, WID_TEXT_1, get_value=dummy
               wlrange[0]=dummy[0]
             end
     'WLMAX':begin
               widget_control, WID_TEXT_2, get_value=dummy
               wlrange[1]=dummy[0]
             end
      'EXIT':begin
               widget_control, base,/DESTROY
               return, wlrange
             end
    'CANCEL':begin
               widget_control, base,/DESTROY
               return, [-1, -1]
             end
    endcase
  endwhile
  widget_control,line_base,/DESTROY
  return,[-1, -1]
end

;=========================================================================================
; IDL Widget Interface Procedures. This Code is automatically
;     generated and should not be modified..., but of course, I went in and changed everything ...
;
; Generated on: 09/21/2003 17:23.11
;
; Changes:
;  09/01/2008 (NP) Added a generation of graphics log file (PostScript) with the line statistics,
;                  resolving power, plots of PSF and PSF variations accross the detector
;
pro WID_BASE_0_event, Event

  common ORDER_INFO, rel_ord_currect, rel_order_selected, abs_ord_selected $
                   , cspec, cs_lines, solution_1D, solution_2D, wl_range $
                   , maxordx, maxordy, npixel, nlist, norders, obase, oincr $
                   , file_prefix, polynom_power, is_polariz, pol_suffix, min_disp $
                   , max_disp, bad_order, bad_order_flag, res_clip

  common WIDGET_IDs, WID_BASE_0,   WID_TEXT_0, WID_TEXT_1, WID_TEXT_2, WID_TEXT_3 $
                   , WID_TEXT_4,   WID_TEXT_5, WID_TEXT_6, WID_TEXT_7, WID_TEXT_8 $
                   , WID_DRAW_0,   WID_BUTTON_0, WID_BUTTON_1, WID_BUTTON_3, WID_BUTTON_4 $
                   , WID_BUTTON_5, WID_BUTTON_7, WID_MAX_VAL, WID_MIN_VAL, graph0

  common  BASE_SIZE, base_x_size, base_y_size


  wTarget = (widget_info(Event.id,/NAME) eq 'TREE' ?  $
      widget_info(Event.id, /tree_root) : event.id)


  if(Event.id eq Event.top) then begin ; Resize event or the main base was killed
    widget_control,Event.id,get_uvalue=uvalue
    if(uvalue eq 'RESIZE') then begin
      widget_control,WID_BASE_0,TLB_GET_SIZE=base_size,UPDATE=0
      base_x_size = base_size[0]>1024
      base_y_size = base_size[1]>600
;    print,base_x_size,base_y_size
      Widget_Control, WID_DRAW_0, DRAW_XSIZE = (base_x_size- 20) $
                                , DRAW_YSIZE = (base_y_size-130)
      base_size=[base_x_size,base_y_size]
      widget_control,WID_BASE_0,/UPDATE
      display, cspec, min=min_disp, max=max_disp
      return
    endif
  endif


  wWidget =  Event.top

  case wTarget of

;-----------------------------------------------------------------------------------------
    Widget_Info(wWidget, FIND_BY_UNAME='1D_SOLUTION'): begin
      if(Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON') then begin

        wlrange=reform(wl_range[*,rel_order_selected])       ; Take the range stored
        disp_poly = 0
        if(oincr ne 0) then abs_order = oincr * rel_order_selected + obase $
        else if(not keyword_set(abs_order)) then abs_order = rel_order_selected
        if(n_elements(solution_2D) gt 0) then begin          ; If 2D solution exists, get the range
          mkwave, w, solution_2D, abs_order                  ; Produce wavelength scale for selected order
          polynom_power = round(solution_2D[8])              ; Get the actual polynomial power in X
          disp_poly = poly_fit(dindgen(npixel), w $          ; Estimate of 1D solution
                                              , polynom_power, /DOUBLE)
          wlrange = minmax(w)                                ; Wavelength limits
          w=0
        endif

; Verify order number and wavelength range
        wlrange = setWLrange(rel_order_selected, abs_order, wlrange, GROUP_LEADER=WID_BASE_0)
        if(wlrange[0] lt 0 or wlrange[1] lt 0) then goto,flat_order
        wl_range[*,rel_order_selected]=wlrange

; Do automatic line pre-ID based on lab spectrum
        cs = reform(cspec[*,rel_order_selected])
        if(min(cs) eq max(cs)) then begin
          bad_order[rel_order_selected]=1
          goto,flat_order
        endif
        if(wlrange[0] gt 0.d0 and wlrange[1] gt wlrange[0]) then begin
          nline=0
          if(n_elements(cs_lines) gt 0) then begin           ; Check if have any identifications so far
            iline = where(cs_lines.order eq rel_order_selected, nline) ; Check if we have done this order before
            if(nline gt 0) then begin                        ; We will start from the last run
              line = cs_lines[iline]                         ; Subset of lines for selected order
            endif
          endif

; Check if we did this order before
          current_save_file=file_prefix+'_order_'+strmid(strtrim(rel_order_selected+1000,2),1)+pol_suffix+'.sav'
          if(file_test(current_save_file)) then begin
            restore,current_save_file,/relaxed_structure_assignment ; Restore the structure
            wlrange=poly([0.d0,double(npixel)],b)                   ; Get the wavelength range
            ii=where(line.posm gt 0 and line.posm lt npixel, nline) ; Skip non-identified lines
            if(nline gt 0) then line=line[ii]
            if(n_elements(cs_lines) eq 0) then cs_lines=line
            disp_poly=b
          endif

          if(n_elements(disp_poly) gt 0) then polynom_power=(n_elements(disp_poly)-1)>1
          if(n_elements(wldev) eq 0) then wldev=0.01d0

          wlrange=[wlrange[0]-1.d0, wlrange[1]+1.d0]       ; Stretch limits a bit

          if(n_elements(cs_lines) eq 0 or nline eq 0) then begin
            wllist=reform((kpno_thar(wrange=wlrange))[*,0])  ; Get lab spectrum
            nlist=n_elements(wllist)                         ; Size of the line list
            jmax=where(cs[1:npixel-2] gt cs[0:npixel-3] and $; Locate main local maxima
                       cs[1:npixel-2] gt cs[2:npixel-1], njmax)
            jjj=(reverse(sort(cs[jmax+1])))[0:(njmax<nlist)-1]
            jmax=jmax[jjj]+1
            njmax=n_elements(jmax)
            height=median(cs[jmax])
            jmin=where(cs[1:npixel-2] lt cs[0:npixel-3] and $; Locate all local minima
                       cs[1:npixel-2] lt cs[2:npixel-1], njmin)+1

            line=[{csline,-1.d0,-1.d0,-1.d0,-1.d0,0,0,'G',0.d0,0,0.d0,0}] ; Initialize line structure
            nline=0                                          ; Next line to add
            for j=0,njmax-1 do begin                         ; Loop through peaks
              if(jmax[j] gt 10 and jmax[j] lt npixel-10 and $; Ignore peaks close to the edge
                 jmax[j] gt min(jmin) and jmax[j] lt max(jmin)) then begin
                xx1=max(where(jmax[j] gt jmin)) & xx1=jmin[xx1]>(jmax[j]-10)
                xx2=min(where(jmax[j] lt jmin)) & xx2=jmin[xx2]<(jmax[j]+10)
                xxx=dindgen(xx2-xx1+1)+xx1
                                       ; Check if there is an overlap with earlier identification
                if(nline gt 0) then ii=where(jmax[j] ge line[0:nline-1].xfirst and $
                                             jmax[j] le line[0:nline-1].xlast, nii) $
                else                nii=0
                if(jmax[j]-xx1 ge 3 and $
                   xx2-jmax[j] ge 3 and nii eq 0) then begin ; If the part above the bottom is sufficiently
                  yapp=disp_gaussfit(xxx,cs[xxx],a,xxx)      ; wide and away from the edges do the Gaussian fit
                  posm=a[1]                                  ; Measured position (Center of Gaussian)
                  if(n_elements(disp_poly) gt 1) then wlc=poly(posm,disp_poly) $ ; Computed wavelength
                  else                                wlc=-1.d0
                  if(posm gt 0 and posm lt npixel-1) then begin  ; Computed position
                    if(n_elements(disp_poly) gt 1) then begin    ; If we have a solution, fine
                      if(min(abs(line.posm-posm)) gt 0.5 and $
                         min(abs(wllist-wlc),i) lt 0.02) then begin
                        p=fz_roots([disp_poly[0]-wllist[i],disp_poly[1:polynom_power]],/DOUBLE)
                      endif else p=posm                          ; Fake if we don't
                      posc=p((where(p ge 0 and p lt npixel-1))[0])
                      line[nline].xfirst=min(xxx)
                      line[nline].xlast=max(xxx)
                      line[nline].posm=posm
                      line[nline].posc=posc
                      if(n_elements(disp_poly) gt 1) then line[nline].wll= wllist[i] $
                      else                                line[nline].wll=-1.d0
                      line[nline].wlc =wlc
                      line[nline].height=abs(a[0])
                      line[nline].width=abs(a[2])
                      line[nline].approx='G'
                      line[nline].flag=0
                      line[nline].order=rel_order_selected
; Make sure this line is not a defect or a cosmic
                      if(line[nline].width gt 1 and line[nline].width lt 8 $
;                        and line[nline].height gt height) then begin
                        ) then begin
                        line=[line,{csline,-1.d0,-1.d0,-1.d0,-1.d0,0,0,'G',0.d0,0,0.d0,0}]
                        nline=nline+1
                      endif
                    endif
                  endif
                endif
              endif
            endfor
            if(nline gt 1) then begin
              ii=sort(line[0:nline-1].posm)                    ; Sort identified lines
              line=line[ii]                                    ; Trim the list
            endif

; If we have a first cut on solution, get rid of obvious outliers
            if(n_elements(disp_poly) gt 1 and nline gt polynom_power) then begin
              b=poly(line.posm,disp_poly)-line.wll
              ii=[-1,uniq(line.wll)]
              nii=n_elements(ii)-1
              l=0
              for i=0,nii-1 do begin
                x=min(abs(b[ii[i]+1:ii[i+1]]),ll)
                l=[l,ll+ii[i]+1]
              endfor
              if(nii gt 1) then begin
                l=l[1:nii]
                line=line[l]
                b=(poly(line.posm,disp_poly)-line.wll)/line.wll*2.9979246d8
                ii=where(abs(b) lt 2000., nline)
              endif else nline=0
              if(nline gt 0) then line=line[ii] else line=0
            endif
          endif
        endif

        disp_poly=disp(cs,log=log,line=line,wl_range=wlrange $
                         ,file=file,poly_power=polynom_power, GROUP_LEADER=wGroup $
                         ,widget_title='Fit wavelength scale to order '+strtrim(rel_order_selected,2))

        if(n_elements(disp_poly) gt 1) then begin
          line.order=rel_order_selected
          abs_order_selected=abs_order
          wl_range[*,rel_order_selected]=poly([0.d0,npixel],disp_poly); wlrange
          if(n_elements(cs_lines) eq 0) then begin           ; New identified line list
            cs_lines=line
          endif else if(n_elements(cs_lines) gt 0) then begin; Existing identified lines
            iline = where(cs_lines.order ne rel_order_selected, nline); Check if we have done this order before
            if(nline gt 0) then begin                        ; We will start from the last run
              cs_lines=[cs_lines[iline],line]                ; Subset of lines for selected order
            endif else cs_lines=line                         ; Only this order was done before
            if(oincr gt 0) then begin                        ; Sort in increasing absolute order
              cs_lines=cs_lines[        sort(cs_lines.order) ]
            endif else begin
              cs_lines=cs_lines[reverse(sort(cs_lines.order))]
            endelse
          endif
          b=disp_poly                                        ; Fit will be saved under historical name
          save,line,b,file=current_save_file                 ; Save the fit to this order

; Check if we can do 2D solution already
          m_wave=cs_lines.wll
          m_pix =cs_lines.posm
          m_flag=cs_lines.flag
          m_ord =cs_lines.order*oincr+obase
          nord=round(max(m_ord)-min(m_ord)+1)           ; Range of orders for which we have solutions
          ii=sort(m_ord)
          m_wave=m_wave[ii]
          m_pix = m_pix[ii]
          m_flag=m_flag[ii]
          m_ord = m_ord[ii]
          if(nord gt 1) then begin
            obase_tmp=obase<(obase+oincr*(norders-1)) ; min(m_ord)
            maxordx=5                                   ; Polynomial power in dispersion direction
            maxordy=5                                   ; Polynomial power in order direction
            disp_2d,m_wave,m_pix,m_ord,m_flag,maxordx $ ; 2D solver
                   ,maxordy,npixel,nord,obase_tmp,wvc,resid_pix=resid_pix,resid_m_s=resid_m_s
            wldev=stddev(resid_m_s)
            widget_control,WID_TEXT_6,set_value=string(mean(abs(resid_pix)),FORM='(F8.5)')
            widget_control,WID_TEXT_7,set_value=string(mean(abs(resid_m_s)),FORM='(F6.0)')
            if(n_elements(wvc) gt 0) then begin
              solution_2D=wvc
              widget_control, WID_BUTTON_1,sensitive=1
              widget_control, WID_BUTTON_3,sensitive=1
              widget_control, WID_BUTTON_4,sensitive=1
            endif
;            i=analyze_residuals(m_pix,m_ord,m_flag,npixel,nord,obase,resid_pix,resid_m_s $
;                               ,good=jgood,group_leader=wWidget)
            i=analyze_residuals(m_wave,m_pix,m_ord,m_flag,maxordx $
                               ,maxordy,npixel,nord,obase,oincr,resid_pix,resid_m_s $
                               ,good=jgood,GROUP_LEADER=wWidget,CLIP=res_clip)
            if(i gt 0) then begin
              cs_lines=cs_lines[ii[jgood]]
              cs_lines.flag=0
;              resid_m_s=resid_m_s[jgood]
;              resid_pix=resid_pix[jgood]
              m_wave=cs_lines.wll
              m_pix =cs_lines.posm
              m_flag=cs_lines.flag
              m_ord =cs_lines.order*oincr+obase
              ii=sort(m_ord)
              m_wave=m_wave[ii]
              m_pix = m_pix[ii]
              m_flag=m_flag[ii]
              m_ord = m_ord[ii]
              nord=round(max(m_ord)-min(m_ord)+1)         ; Range of orders for which we have solutions
              obase_tmp=obase<(obase+oincr*(norders-1)); min(m_ord)
              if(nord gt 1) then disp_2d,m_wave,m_pix,m_ord,m_flag,maxordx $ ; 2D solver
                                        ,maxordy,npixel,nord,obase_tmp,wvc $
                                        ,resid_pix=resid_pix,resid_m_s=resid_m_s
              wldev=stddev(resid_m_s)
              widget_control,WID_TEXT_6,set_value=string(mean(abs(resid_pix)),FORM='(F8.5)')
              widget_control,WID_TEXT_7,set_value=string(mean(abs(resid_m_s)),FORM='(F6.0)')
            endif
            widget_control,WID_TEXT_8,set_value=strtrim(n_elements(cs_lines),2)
            if(n_elements(uniq(m_ord)) le 6) then  widget_control, WID_BUTTON_5,sensitive=0 $; Adjust the status
            else                                     widget_control, WID_BUTTON_5,sensitive=1  ; of the report button
          endif
        endif

; Refresh order plot
        wset,graph0
        display, cspec, min=min_disp,max=max_disp       ; Plot spectrum
        if(n_elements(cs_lines) gt 0) then begin
          for iord=0,norders-1 do begin
            i=where(cs_lines.order eq iord,ni)
            if(ni gt 0) then oplot, !x.crange, [iord,iord], col=160 $
            else             oplot, !x.crange, [iord,iord], col= 80
          endfor
        endif else begin
          for iord = 0, norders-1 do begin
            oplot, !x.crange, [iord,iord], col= 80
          endfor
        endelse

; Verify if the relation between relative and absolute order numbers is not set and
; have more than two orders processed, we can define this relation here
        if(n_elements(disp_poly) gt 1) then begin
          omin=min(cs_lines.order,i1)
          omax=max(cs_lines.order,i2)
          if(oincr eq 0 and omax ne omin) then begin
            if(i1 lt i2) then oincr=1 else if(i1 gt i2) then oincr=-1
            obase=abs_order_selected-oincr*rel_order_selected
            widget_control,WID_TEXT_4,set_value=strtrim(obase,2)
            widget_control,WID_TEXT_5,set_value=strtrim(oincr,2)
          endif
        endif
      endif
flat_order:
      wset,graph0
    end

;-----------------------------------------------------------------------------------------
; Rebuild 2D solution
    Widget_Info(wWidget, FIND_BY_UNAME='2D_SOLUTION'): begin
      if(Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON') then begin
        m_wave=cs_lines.wll
        m_pix =cs_lines.posm
        m_flag=cs_lines.flag
        m_ord =cs_lines.order*oincr+obase
        ii=sort(m_ord)
        m_wave=m_wave[ii]
        m_pix = m_pix[ii]
        m_flag=m_flag[ii]
        m_ord = m_ord[ii]
        nord=round(max(m_ord)-min(m_ord)+1)         ; Range of orders for which we have solutions
        obase_tmp=obase<(obase+oincr*(norders-1)); min(m_ord)
        disp_2d,m_wave,m_pix,m_ord,m_flag,maxordx $ ; 2D solver
               ,maxordy,npixel,nord,obase_tmp,wvc,resid_pix=resid_pix,resid_m_s=resid_m_s
        wldev=stddev(resid_m_s)
        widget_control,WID_TEXT_6,set_value=string(mean(abs(resid_pix)),FORM='(F8.5)')
        widget_control,WID_TEXT_7,set_value=string(mean(abs(resid_m_s)),FORM='(F6.0)')
        if(n_elements(wvc) gt 0) then begin
          solution_2D=wvc
          widget_control, WID_BUTTON_1,sensitive=1
          widget_control, WID_BUTTON_3,sensitive=1
          widget_control, WID_BUTTON_4,sensitive=1
        endif
;        i=analyze_residuals(m_pix,m_ord,m_flag,npixel,nord,obase,resid_pix,resid_m_s $
;                                 ,good=jgood,group_leader=wWidget)
        i=analyze_residuals(m_wave,m_pix,m_ord,m_flag,maxordx $
                           ,maxordy,npixel,nord,obase,oincr,resid_pix,resid_m_s $
                           ,good=jgood,GROUP_LEADER=wWidget,CLIP=res_clip)
        if(i gt 0) then begin
          cs_lines=cs_lines[ii[jgood]]
          cs_lines.flag=0
;          resid_m_s=resid_m_s[jgood]
;          resid_pix=resid_pix[jgood]
          m_wave=cs_lines.wll
          m_pix =cs_lines.posm
          m_flag=cs_lines.flag
          m_ord =cs_lines.order*oincr+obase
          ii=sort(m_ord)
          m_wave=m_wave[ii]
          m_pix = m_pix[ii]
          m_flag=m_flag[ii]
          m_ord = m_ord[ii]
          nord=round(max(m_ord)-min(m_ord)+1)         ; Range of orders for which we have solutions
          obase_tmp=obase<(obase+oincr*(norders-1)); min(m_ord)
          if(nord gt 2) then disp_2d,m_wave,m_pix,m_ord,m_flag,maxordx $ ; 2D solver
                                    ,maxordy,npixel,nord,obase_tmp,wvc $
                                    ,resid_pix=resid_pix,resid_m_s=resid_m_s
          wldev=stddev(resid_m_s)
          widget_control,WID_TEXT_6,set_value=string(mean(abs(resid_pix)),FORM='(F8.5)')
          widget_control,WID_TEXT_7,set_value=string(mean(abs(resid_m_s)),FORM='(F6.0)')
        endif
        widget_control,WID_TEXT_8,set_value=strtrim(n_elements(cs_lines),2)
        if(n_elements(uniq(m_ord)) le 6) then  widget_control, WID_BUTTON_5,sensitive=0 $; Adjust the status
        else                                     widget_control, WID_BUTTON_5,sensitive=1  ; of the report button
; Refresh order plot
        wset,graph0
        display, cspec, min=min_disp,max=max_disp       ; Plot spectrum
        if(n_elements(cs_lines) gt 0) then begin
          for iord=0,norders-1 do begin
            i=where(cs_lines.order eq iord,ni)
            if(ni gt 0) then oplot, !x.crange, [iord,iord], col=160 $
            else             oplot, !x.crange, [iord,iord], col= 80
          endfor
        endif else begin
          for iord = 0, norders-1 do begin
            oplot, !x.crange, [iord,iord], col= 80
          endfor
        endelse
      endif
    end

;-----------------------------------------------------------------------------------------
; AutoID all relevant lines using 2D solution
    Widget_Info(wWidget, FIND_BY_UNAME='AUTO_ID'): begin
      if(Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON') then begin
        widget_control,/HOURGLASS


        for iord=0,norders-1 do begin
          if(not bad_order[iord]) then begin
            abs_order=oincr*iord+obase                         ; Absolute order number
            mkwave, w, solution_2D, abs_order                  ; Produce wavelength polynomial
            polynom_power = round(solution_2D[8])              ; for the selected order
            disp_poly = poly_fit(dindgen(npixel), w $          ; Construct 1D wavelength scale
                                                , polynom_power, /DOUBLE)
            wlrange = minmax(w)                                ; Wavelength limits
            w=0
            wlrange=[wlrange[0]-1.d0, wlrange[1]+1.d0]         ; Stretch limits a bit
            cs = reform(cspec[*,iord])                         ; ThAr spectrum for a single order

            wllist=reform((kpno_thar(wrange=wlrange))[*,0])    ; Get lab spectrum
            nlist=n_elements(wllist)                           ; Size of the line list
            jmax=where(cs[1:npixel-2] gt cs[0:npixel-3] and $  ; Locate primary local maxima
                       cs[1:npixel-2] gt cs[2:npixel-1], njmax)
            if(njmax eq 0) then goto,skip_order

; Restrict the lines to significant signal
            limit=median(cs[jmax+1])
            jjj=where(cs[jmax+1] gt limit, njmax)
            jmax=jmax[jjj]
            if(njmax eq 0) then goto,skip_order

; Reject maxima which are too close:
; We now may have continuous (in pixel order) groups of lines which are too close.
; All are marked as tight except one which is believed to be OK ("which" statement that sets jjj
; looks only at the neighour on the right). We have to rarefy those groups by dropping the weaker lines
; until no tight pair is left. This why iterations are needed.
tight:
            jgood=lindgen(njmax)
            jjj=where(jmax[1:njmax-1]-jmax[0:njmax-2] gt 6 $  ; Find tight maxima
                      , njmax, COMPLEMENT=jcompl, NCOMPLEMENT=njcompl) 
            if(njmax eq 0) then goto,skip_order                ; Bad order: all lines are one blob
            jjj=[jjj,n_elements(jmax)-1L]                      ; Mark the last line as good
            njmax=njmax+1                                      ; Update number of lines marked as good

            if(njcompl gt 0) then begin
              iii=uniq(jcompl-lindgen(njcompl))    ; rightmost boundaries of continuous index stretches
              iii=[-1L,iii]                        ; first boundary
              for j=0L,n_elements(iii)-2L do begin ; Loop through tight groups, jcompl refers to the left peak
                ifirst=jcompl[iii[j]+1L]           ; Starting number of the group
                ilast =jcompl[iii[j+1]]+1          ; Ending number of the group
                mx=min(cs[jmax[ifirst:ilast]+1],imin)
                imin=imin+ifirst
                jgood[imin]=-1
              endfor
              jgood=jgood[where(jgood ge 0)]
              jmax=jmax[jgood]
              njmax=n_elements(jmax)
              goto,tight
            endif
            jmax=jmax[jjj]+1
            jjj=(reverse(sort(cs[jmax])));[0:(njmax<nlist)-1]  ; Sort lines in decreasing strength
            njmax=n_elements(jmax)
            height=median(cs[jmax])
            jmin=where(cs[1:npixel-2] lt cs[0:npixel-3] and $  ; Locate all local minima
                       cs[1:npixel-2] lt cs[2:npixel-1], njmin)+1
            jmin=[0,jmin,npixel-1]
            ii=where(cs_lines.order eq iord, nline)
            if(nline gt 0) then begin
              refit,cs_lines[ii],cs,error=error ; Quickly check if manually identified lines are OK
              jgood=where(error lt 0.1d0, ngood)
              if(ngood gt 0) then $ ; Have some good lines
                line=[cs_lines[ii[jgood]], {csline,-1.d0,-1.d0,-1.d0,-1.d0,0,0,'G',0.d0,0,0.d0,0}] $
              else $                ; No good lines; start with an empty structure
                line=[{csline,-1.d0,-1.d0,-1.d0,-1.d0,0,0,'G',0.d0,0,0.d0,0}]
              nline=ngood
            endif else begin
              line=[{csline,-1.d0,-1.d0,-1.d0,-1.d0,0,0,'G',0.d0,0,0.d0,0}] ; Initialize line structure
              nline=0                                          ; Next line to add
            endelse
            for j=0L,njmax-1 do begin                          ; Loop through peaks
              if(jmax[j] gt 10 and jmax[j] lt npixel-10) then begin ; Ignore peaks close to the edge
                xx1=max(where(jmax[j] gt jmin)) & xx1=jmin[xx1]>(jmax[j]-10)
                xx2=min(where(jmax[j] lt jmin)) & xx2=jmin[xx2]<(jmax[j]+10)
                xxx=dindgen(xx2-xx1+1)+xx1
                                     ; Check if there is an overlap with earlier identification
                if(nline gt 0) then ii=where(jmax[j] ge line[0:nline-1].xfirst and $
                                             jmax[j] le line[0:nline-1].xlast, nii) $
                else                nii=0
                if(jmax[j]-xx1 ge 3 and $
                   xx2-jmax[j] ge 3 and nii eq 0) then begin   ; If the peak is in the center and
                  yapp=disp_gaussfit(xxx,cs[xxx],a,xxx)        ; sufficiently wide do Gaussian fit
                  error=total(abs(yapp-cs[xxx]))/total(cs[xxx])
                  slop=abs(a[4]/(max(cs[xxx])-min(cs[xxx])))
                  posm=a[1]                                    ; Measured position (Center of Gaussian)
                  if(n_elements(disp_poly) gt 1) then wlc=poly(posm,disp_poly) $ ; Computed wavelength
                  else                                wlc=-1.d0
                  if(posm gt 0 and posm lt npixel-1) then begin; Computed position
                    if(min(abs(line.posm-posm)) gt 0.5 and $
                       min(abs(wllist-wlc),i) lt 0.02) then begin
                      p=fz_roots([disp_poly[0]-wllist[i],disp_poly[1:polynom_power]],/DOUBLE)
                    endif else p=posm                          ; Fake if we don't
                    posc=p((where(p ge 0 and p lt npixel-1))[0])
                    line[nline].xfirst=min(xxx)
                    line[nline].xlast=max(xxx)
                    line[nline].posm=posm
                    line[nline].posc=posc
                    line[nline].wll= wllist[i]
                    line[nline].wlc =wlc
                    line[nline].height=abs(a[0])
                    line[nline].width=abs(a[2])
                    line[nline].approx='G'
                    line[nline].flag=1
                    line[nline].order=iord
; Make sure this line is not a defect or a cosmic
                    if($
;                       a[3] le min(cs[xxx]) and a[1] gt xx1+3 and a[1] lt xx2-3 and $
;                       error lt 0.1d0 and slop lt 0.005d0 and $
                       line[nline].width gt 1 and $
                       line[nline].width lt 8) then begin
;                       line[nline].width lt 8 and $
;                       line[nline].height gt height) then begin
                      line=[line,{csline,-1.d0,-1.d0,-1.d0,-1.d0,0,0,'G',0.d0,1,0.d0,0}]
                      nline=nline+1
                    endif
                  endif
                endif
              endif
            endfor
            if(nline gt 1) then begin
              ii=sort(line[0:nline-1].posm)                    ; Sort identified lines
              line=line[ii]                                    ; Trim the list
            endif
            if(iord eq 0) then all_lines=line else all_lines=[all_lines,line]

; If we have a first cut on solution, get rid of obvious outliers
            if(n_elements(disp_poly) gt 1 and nline gt polynom_power) then begin
              b=line.wlc-line.wll
              ii=[-1,uniq(line.wll)]
              nii=n_elements(ii)-1
              l=0
              for i=0,nii-1 do begin
                x=min(abs(b[ii[i]+1:ii[i+1]]),ll)
                l=[l,ll+ii[i]+1]
              endfor
              l=l[1:nii]
              line=line[l]
              b=(line.wlc-line.wll)/line.wll*2.9979246d8
;              dev=stddev(b)
              ii=where(abs(b) lt 2000., nii)
              if(nii gt 0) then line=line[ii]
              nline=n_elements(line)
            endif
            if(n_elements(cs_lines) eq 0) then begin           ; New identified line list
              cs_lines=line
            endif else if(n_elements(cs_lines) gt 0) then begin; Existing identified lines
              iline = where(cs_lines.order ne iord, nline)     ; Check if we have done this order before
              if(nline gt 0) then begin                        ; We will start from the last run
                cs_lines=[cs_lines[iline],line]                ; Subset of lines for selected order
              endif else cs_lines=line                         ; Only this order was done before
              if(oincr gt 0) then begin                        ; Sort in increasing absolute order
                cs_lines=cs_lines[        sort(cs_lines.order) ]
              endif else begin
                cs_lines=cs_lines[reverse(sort(cs_lines.order))]
              endelse
            endif
          endif
skip_order:
        endfor
        save,all_lines,file=file_prefix+pol_suffix+'_posm.sav'

        m_wave=cs_lines.wll
        m_pix =cs_lines.posm
        m_flag=cs_lines.flag
        m_ord =cs_lines.order*oincr+obase
        nord=round(max(m_ord)-min(m_ord)+1)         ; Range of orders for which we have solutions
        resid_pix=cs_lines.posm-cs_lines.posc
        resid_m_s=(cs_lines.wll-cs_lines.wlc)/cs_lines.wll*2.9979246d8
;        i=analyze_residuals(m_pix,m_ord,m_flag,npixel,nord,obase,resid_pix,resid_m_s $
;                                 ,good=jgood,group_leader=wWidget)
        i=analyze_residuals(m_wave,m_pix,m_ord,m_flag,maxordx $
                           ,maxordy,npixel,nord,obase,oincr,resid_pix,resid_m_s $
                           ,good=jgood,GROUP_LEADER=wWidget,CLIP=res_clip)
        if(i gt 0) then begin
          cs_lines=cs_lines[jgood]
          cs_lines.flag=0
          resid_m_s=resid_m_s[jgood]
          resid_pix=resid_pix[jgood]
          m_wave=cs_lines.wll
          m_pix =cs_lines.posm
          m_flag=cs_lines.flag
          m_ord =cs_lines.order*oincr+obase
          ii=sort(m_ord)
          m_wave=m_wave[ii]
          m_pix = m_pix[ii]
          m_flag=m_flag[ii]
          m_ord = m_ord[ii]
          nord=round(max(m_ord)-min(m_ord)+1)         ; Range of orders for which we have solutions
          obase_tmp=obase<(obase+oincr*(norders-1)); min(m_ord)
          if(nord gt 1) then disp_2d,m_wave,m_pix,m_ord,m_flag,maxordx $ ; 2D solver
                                    ,maxordy,npixel,nord,obase_tmp,wvc $
                                    ,resid_pix=resid_pix,resid_m_s=resid_m_s
          wldev=stddev(resid_m_s)
          widget_control,WID_TEXT_6,set_value=string(mean(abs(resid_pix)),FORM='(F8.5)')
          widget_control,WID_TEXT_7,set_value=string(mean(abs(resid_m_s)),FORM='(F6.0)')
        endif
        widget_control,WID_TEXT_8,set_value=strtrim(n_elements(cs_lines),2)
        if(n_elements(uniq(m_ord)) le 6) then  widget_control, WID_BUTTON_5,sensitive=0 $; Adjust the status
        else                                   widget_control, WID_BUTTON_5,sensitive=1  ; of the report button

        wset,graph0
        display, cspec, min=min_disp,max=max_disp       ; Plot spectrum
        if(n_elements(cs_lines) gt 0) then begin
          for iord=0,norders-1 do begin
            i=where(cs_lines.order eq iord,ni)
            if(ni gt 0) then oplot, !x.crange, [iord,iord], col=160 $
            else             oplot, !x.crange, [iord,iord], col= 80
          endfor
        endif else begin
          for iord = 0, norders-1 do begin
            oplot, !x.crange, [iord,iord], col= 80
          endfor
        endelse
      endif
    end

;-----------------------------------------------------------------------------------------
; Add 2D wavelength solution to the selected ECH files
    Widget_Info(wWidget, FIND_BY_UNAME='MOD_ECH'): begin
      if(Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON') then begin
        widget_control,/HOURGLASS
        SDtoFWHM=2.d0*sqrt(2.d0*alog(2.d0))
        resolution=dblarr(norders)
        for iord=0,norders-1 do begin
          abs_order=oincr*iord+obase                         ; Absolute order number
          mkwave, w, solution_2D, abs_order                  ; Produce wavelength scale for selected order
          disp_poly = poly_fit(dindgen(npixel), w $          ; Get 1D solution
                                              , polynom_power, /DOUBLE)
          w=mean(w)
          i=where(cs_lines.order eq iord and cs_lines.flag eq 0, n)
          if(n gt 0) then begin
            line=cs_lines[i]
            mean_width=total(line[i].width)/n
            delta_mean_width=sqrt(total((mean_width-line[i].width)^2))/n
            ii=where(abs(mean_width-line[i].width) lt 5*delta_mean_width, n)
            if(n gt 0) then line=line[ii]
            n=n_elements(line)
            mean_width=total(line.width)/n
            resolution[iord]=w/(abs(disp_poly[1])*mean_width)
          endif
        endfor

        files=dialog_pickfile(dialog_parent=wWIDGET,/MULTIPLE_FILES,/MUST_EXIST $
                             ,title='Select ECH file(s)',filter='*.ech')
        nfiles=n_elements(files)

        if(is_polariz) then orders=(indgen(norders*2)/2)*oincr+obase $
        else                orders= indgen(norders)     *oincr+obase

        solution_2D[2]=n_elements(orders)
        for i=0,nfiles-1 do begin
          if(file_test(files[i])) then begin
            rdech,e,files[i]
            modech,files[i],orders=orders,wave=solution_2D $
                    ,resolut=round(mean(resolution))
            rdech,e,files[i]                                ;check for heliocentric
            helcor=hierarch(e.head,'E_HELCOR',COUNT=cnt)    ;velocity in header
            if cnt gt 0 then begin                          ;and 'copy' it to
              modech,files[i],barycorr=helcor               ;BARYCORR keyword
            endif
          endif
        endfor
      endif
    end

;-----------------------------------------------------------------------------------------
; Selection window for spectral order
    Widget_Info(wWidget, FIND_BY_UNAME='WID_DRAW_0'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_DRAW' )then begin
        if( Event.type eq 0 )then begin ; Zero means "button pressed"
          col= 80
          if(n_elements(cs_lines) gt 0) then begin
            ii=where(cs_lines.order eq rel_order_selected, nii)
            if(nii gt 0) then col=160 else col= 80
          endif
          if(rel_order_selected ge 0) then oplot, !x.crange $
                                                , [rel_order_selected, rel_order_selected] $
                                                , col=col
          x = Event.X
          y = Event.Y
          y = convert_coord([x, y, 0], /DEVICE, /TO_DATA)
          y = y[1]
          rel_order_selected = 0>round(y)<(norders-1)
          abs_order_selected = oincr*rel_order_selected+obase
          widget_control, WID_BUTTON_0, sensitive=1
          widget_control, WID_BUTTON_7, sensitive=1
          if(bad_order[rel_order_selected]) then begin
            if(not bad_order_flag) then begin
              bad_order_flag=1
              widget_control, WID_BUTTON_0,sensitive=0
              widget_control, WID_BUTTON_7,set_value='Bad order ', /set_button
            endif
          endif else begin
            if(bad_order_flag) then begin
              bad_order_flag=0
              widget_control, WID_BUTTON_0,sensitive=1
              widget_control, WID_BUTTON_7,set_value='Good order', set_button=0
            endif
          endelse
          Widget_Control, WID_TEXT_2, set_value = strtrim(rel_order_selected, 2)
          Widget_Control, WID_TEXT_3, set_value = strtrim(abs_order_selected, 2)
          oplot,!x.crange,[rel_order_selected, rel_order_selected]
        endif else if(Event.type eq 2) then begin   ; Just dragging the cursor
          x = Event.X
          y = Event.Y
          y = (convert_coord([x, y, 0], /DEVICE, /TO_DATA))[1]
          rel_order_current = 0>round(y)<(norders-1)
          if(bad_order[rel_order_current]) then begin
            if(not bad_order_flag) then begin
              bad_order_flag=1
              widget_control, WID_BUTTON_0,sensitive=0
              widget_control, WID_BUTTON_7,set_value='Bad order ', /set_button
            endif
          endif else begin
            if(bad_order_flag) then begin
              bad_order_flag=0
              widget_control, WID_BUTTON_0,sensitive=1
              widget_control, WID_BUTTON_7,set_value='Good order', set_button=0
            endif
          endelse
          abs_order_current = oincr*rel_order_current+obase
          Widget_Control, WID_TEXT_0, set_value = strtrim(rel_order_current, 2)
          Widget_Control, WID_TEXT_1, set_value = strtrim(abs_order_current, 2)
        endif
      endif else if(Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_TRACKING') then begin
        if(Event.Enter eq 0 and rel_order_selected ge 0)then begin ; Cursor leaves the drawing area
          if(bad_order[rel_order_selected]) then begin
            if(not bad_order_flag) then begin
              bad_order_flag=1
              widget_control, WID_BUTTON_0,sensitive=0
              widget_control, WID_BUTTON_7,set_value='Bad order ', /set_button
            endif
          endif else begin
            if(bad_order_flag) then begin
              bad_order_flag=0
              widget_control, WID_BUTTON_0,sensitive=1
              widget_control, WID_BUTTON_7,set_value='Good order', set_button=0
            endif
          endelse
        endif
      endif
    end

;-----------------------------------------------------------------------------------------
; User may have edited the selected absolute order number
    Widget_Info(wWidget, FIND_BY_UNAME='SELECTED_ABSOLUTE_ORDER'): begin
      if(Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_KBRD_FOCUS') then begin
        widget_control, WID_TEXT_3, get_value=dummy
        dummy=fix(dummy[0])
        if(n_elements(abs_order_selected) gt 0) then begin
          if(abs_order_selected ne dummy) then begin
            ii=where(cs_lines.order eq abs_order_selected, nii)
            if(nii gt 0) then cs_lines[ii].order=abs_order_selected
            abs_order_selected=dummy
          endif
        endif
      endif
    end

;-----------------------------------------------------------------------------------------
; User may have explicitely set the base order
    Widget_Info(wWidget, FIND_BY_UNAME='BASE_ORDER'): begin
      if(Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_KBRD_FOCUS') then begin
        widget_control, WID_TEXT_4, get_value=dummy
        dummy=fix(dummy[0])
        if(obase ne dummy) then obase=dummy
      endif
    end

;-----------------------------------------------------------------------------------------
; User may have explicitely set the order increment
    Widget_Info(wWidget, FIND_BY_UNAME='ORDER_INCREMENT'): begin
      if(Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_KBRD_FOCUS') then begin
        widget_control, WID_TEXT_5, get_value=dummy
        dummy=fix(dummy[0])
        if(oincr ne dummy) then oincr=dummy
      endif
    end

;-----------------------------------------------------------------------------------------
; Mark/unmark current order as bad
    Widget_Info(wWidget, FIND_BY_UNAME='BAD_ORDER'): begin
      if(Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON') then begin
        if(rel_order_selected ge 0 and rel_order_selected lt norders) then begin
          bad_order[rel_order_selected]=(bad_order[rel_order_selected]+1B) mod 2B
          bad_order_flag=1
          if(bad_order[rel_order_selected]) then begin
            widget_control, WID_BUTTON_7, set_value='Bad order '
            widget_control, WID_BUTTON_0, SENSITIVE=0
          endif else begin
            widget_control, WID_BUTTON_7, set_value='Good order'
            widget_control, WID_BUTTON_0, SENSITIVE=1
          endelse
        endif
      endif
    end

;-----------------------------------------------------------------------------------------
; Align model
    Widget_Info(wWidget, FIND_BY_UNAME='ALIGN'): begin
      if(Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON') then begin
        offset=disp_align(cspec,cs_lines,obase,oincr,GROUP_LEADER=wGroup)
;        print,offset
        wset,graph0
        display, cspec, min=min_disp,max=max_disp       ; Plot spectrum
        widget_control,WID_TEXT_4,set_value=strtrim(obase,2)
        if(n_elements(cs_lines) gt 0) then begin
          for iord=0,norders-1 do begin
            i=where(cs_lines.order eq iord,ni)
            if(ni gt 0) then oplot, !x.crange, [iord,iord], col=160 $
            else             oplot, !x.crange, [iord,iord], col= 80
          endfor
        endif else begin
          for iord = 0, norders-1 do begin
            oplot, !x.crange, [iord,iord], col= 80
          endfor
        endelse
      endif
    end

;-----------------------------------------------------------------------------------------
; Adjust top display cut
    Widget_Info(wWidget, FIND_BY_UNAME='MAX_VAL'): begin
      if(Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_SLIDER') then begin
        widget_control,WID_MAX_VAL,get_value=r
        max_disp=round(r[0])
        display, cspec, min=min_disp,max=max_disp       ; Plot spectrum
        if(n_elements(cs_lines) gt 0) then begin
          for iord=0,norders-1 do begin
            i=where(cs_lines.order eq iord,ni)
            if(ni gt 0) then oplot, !x.crange, [iord,iord], col=160 $
            else             oplot, !x.crange, [iord,iord], col= 80
          endfor
        endif else begin
          for iord = 0, norders-1 do begin
            oplot, !x.crange, [iord,iord], col= 80
          endfor
        endelse
        widget_control,WID_MIN_VAL,set_slider_max=max_disp-(max_disp-min_disp)/79.
        widget_control,WID_MIN_VAL,set_value=min_disp
      endif
    end

;-----------------------------------------------------------------------------------------
; Adjust bottom display cut
    Widget_Info(wWidget, FIND_BY_UNAME='MIN_VAL'): begin
      if(Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_SLIDER') then begin
        widget_control,WID_MIN_VAL,get_value=r
        min_disp=round(r[0])
        display, cspec, min=min_disp,max=max_disp       ; Plot spectrum
        if(n_elements(cs_lines) gt 0) then begin
          for iord=0,norders-1 do begin
            i=where(cs_lines.order eq iord,ni)
            if(ni gt 0) then oplot, !x.crange, [iord,iord], col=160 $
            else             oplot, !x.crange, [iord,iord], col= 80
          endfor
        endif else begin
          for iord = 0, norders-1 do begin
            oplot, !x.crange, [iord,iord], col= 80
          endfor
        endelse
        widget_control,WID_MAX_VAL,set_slider_min=min_disp+(max_disp-min_disp)/79.
        widget_control,WID_MAX_VAL,set_value=max_disp
      endif
    end

;-----------------------------------------------------------------------------------------
; Generate a graphic report
    Widget_Info(wWidget, FIND_BY_UNAME='REPORT'): begin
      if(Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON') then begin
        widget_control,/HOURGLASS
        make_log,cspec,cs_lines,oincr,obase,solution_2D,file_prefix $
                ,GROUP_LEADER=wWidget
      endif
    end

;-----------------------------------------------------------------------------------------
; Normal way out
    Widget_Info(wWidget, FIND_BY_UNAME='EXIT'): begin
      if(Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON') then begin
        answ=dialog_message('Do you really want to exit?',TITLE='Wavecal', $
                            /QUESTION,/CANCEL,DIALOG_PARENT=wWIDGET)
        if(answ eq 'Yes') then begin ; Normal exit, save the results
          widget_control,wWidget,/DESTROY
          save,solution_2D,cs_lines,obase,oincr,bad_order,file=file_prefix+pol_suffix+'_2D.sav'
          loadct,0
          !p.position=0
          !p.background=0
          return
        endif else if(answ eq 'Cancel') then begin
          widget_control,wWidget,/DESTROY
          loadct,0
          !p.position=0
          !p.background=0
          return
        endif
      endif
    end
    else:
  endcase
end


;=========================================================================================
pro WID_BASE_0, CS_struct, CS_ech_file, prefix, save_file, POLARIZ=polariz     $
              , REVERSE=flip_wl, UPSIDE_DOWN=flip_orders, GROUP_LEADER=wGroup $
              , _EXTRA=_VWBExtra_
  common ORDER_INFO, rel_ord_currect, rel_order_selected, abs_ord_selected $
                   , cspec, cs_lines, solution_1D, solution_2D, wl_range $
                   , maxordx, maxordy, npixel, nlist, norders, obase, oincr $
                   , file_prefix, polynom_power, is_polariz, pol_suffix, min_disp $
                   , max_disp, bad_order, bad_order_flag, res_clip

  common WIDGET_IDs, WID_BASE_0,   WID_TEXT_0, WID_TEXT_1, WID_TEXT_2, WID_TEXT_3 $
                   , WID_TEXT_4,   WID_TEXT_5, WID_TEXT_6, WID_TEXT_7, WID_TEXT_8 $
                   , WID_DRAW_0,   WID_BUTTON_0, WID_BUTTON_1, WID_BUTTON_3, WID_BUTTON_4 $
                   , WID_BUTTON_5, WID_BUTTON_7, WID_MAX_VAL, WID_MIN_VAL, graph0

  common  BASE_SIZE, base_x_size, base_y_size

  file_prefix = prefix
  rel_order_selected = -1
  oincr = 0
  obase = 0
  device, get_screen_size = scr_sz

  if(keyword_set(polariz)) then begin     ; Use only half of the orders if polariz is set 
    nord=n_elements(CS_struct.spec[0,*])  ; Selecting polariz to be odd of even allows to choose same polarization orders
    cspec=CS_struct.spec[*,indgen(nord/2)*2+(polariz mod 2)]
    pol_suffix=((polariz mod 2) eq 1)?'_pol1':'_pol2'
    is_polariz = 1 eq 1
  endif else begin
    cspec=CS_struct.spec
    pol_suffix=''
    is_polariz= 1 ne 1
  endelse

  polynom_power=3                                     ; Polynom power for single order fit
  maxordx=5                                           ; Polynomial power in dispersion direction
  maxordy=5                                           ; Polynomial power in order direction

  if(keyword_set(save_file)) then begin
    f=save_file
    if(not file_test(f)) then begin
      print,'Save file "'+f+'" was not found'
      return
    endif
    restore,f
    npix=n_elements(cspec[*,0])
    nord=n_elements(cspec[0,*])
    if(keyword_set(flip_wl)) then begin
      x1             =cs_lines.xfirst
      x2             =cs_lines.xlast
      cs_lines.xfirst=npix-1L-x2
      cs_lines.xlast =npix-1L-x1
      cs_lines.posc  =npix-1L-cs_lines.posc
      cs_lines.posm  =npix-1L-cs_lines.posm
    endif
    if(keyword_set(flip_orders)) then begin
      cs_lines.order =nord-1L-cs_lines.order
;      i=sort(cs_lines.order)
;      cs_lines=cs_lines[i]
;      i=0
      obase          =obase+oincr*(nord-1L)
      oincr          =-oincr
    endif
  endif else begin
    f=file_prefix+pol_suffix+'_2D.sav'
    if(file_test(f)) then restore,f
  endelse

  base_x_size = (scr_sz[0]-20)<1024
  base_y_size = (scr_sz[1]-30)<950
  WID_BASE_0 = Widget_Base( GROUP_LEADER=wGroup, UNAME='WID_BASE_0' $
      , XOFFSET=5, YOFFSET=5 $
;      , SCR_XSIZE=base_x_size, SCR_YSIZE=base_y_size $
      , /COLUMN, /TLB_Size_Events, UVALUE='RESIZE' $
      , TITLE='Wavelength calibration '+CS_ech_file, SPACE=3, XPAD=3 ,YPAD=3)

;=========== Information boxes =================================

  WID_BASE_1 = Widget_Base(WID_BASE_0, UNAME='WID_BASE_1', /ROW)

  WID_LABEL_0 = Widget_Label(WID_BASE_1, UNAME='WID_LABEL_0'  $
      , /ALIGN_LEFT, VALUE='Current relative order', XOFFSET=10)

  WID_TEXT_0 = Widget_Text(WID_BASE_1, UNAME='CURRENT_RELATIVE_ORDER' $
      , VALUE=[ '  0' ], XSIZE=3, YSIZE=1)

  WID_LABEL_1 = Widget_Label(WID_BASE_1, UNAME='WID_LABEL_1'  $
      , /ALIGN_LEFT, VALUE='Current absolute order', XOFFSET=10)

  WID_TEXT_1 = Widget_Text(WID_BASE_1, UNAME='CURRENT_ABSOLUTE_ORDER'  $
      , VALUE=[ '  0' ], XSIZE=3, YSIZE=1)

  WID_LABEL_2 = Widget_Label(WID_BASE_1, UNAME='WID_LABEL_2'  $
      , /ALIGN_LEFT, VALUE='Selected relative order', XOFFSET=10)

  WID_TEXT_2 = Widget_Text(WID_BASE_1, UNAME='SELECTED_RELATIVE_ORDER'  $
      , VALUE=[ ' -1' ], XSIZE=3, YSIZE=1)

  WID_LABEL_3 = Widget_Label(WID_BASE_1, UNAME='WID_LABEL_3'  $
      , /ALIGN_LEFT, VALUE='Selected absolute order', XOFFSET=10)

  WID_TEXT_3 = Widget_Text(WID_BASE_1, UNAME='SELECTED_ABSOLUTE_ORDER'  $
      , VALUE=[ ' -1' ], XSIZE=3, YSIZE=1, /EDITABLE, /KBRD_FOCUS_EVENTS)

  WID_LABEL_4 = Widget_Label(WID_BASE_1, UNAME='WID_LABEL_4'  $
      , /ALIGN_LEFT, VALUE='Base order number', XOFFSET=10)

  WID_TEXT_4 = Widget_Text(WID_BASE_1, UNAME='BASE_ORDER'  $
      , VALUE=[ ' -1' ], XSIZE=3, YSIZE=1, /EDITABLE, /KBRD_FOCUS_EVENTS)

  WID_LABEL_5 = Widget_Label(WID_BASE_1, UNAME='WID_LABEL_5'  $
      , /ALIGN_LEFT, VALUE='Order increment', XOFFSET=10)

  WID_TEXT_5 = Widget_Text(WID_BASE_1, UNAME='ORDER_INCREMENT'  $
      , VALUE=[ '  0' ], XSIZE=3 ,YSIZE=1, /EDITABLE, /KBRD_FOCUS_EVENTS)

;===================== Controls =============================
  WID_BASE_0a=Widget_Base(WID_BASE_0, /ROW)
  WID_BASE_0b=Widget_Base(WID_BASE_0a, /COLUMN)

;================ Contrast sliders ==========================
  WID_BASE_3 = Widget_Base(WID_BASE_0a, /COLUMN)

  min_disp=min(round(cspec>0))
  max_disp=max(round((min_disp*9+max(cspec>0))/10.))

  WID_MAX_VAL = Widget_Slider(WID_BASE_3, MAX=max_disp $
            , MIN=min_disp+(max_disp-min_disp)/79., SCROLL=1, UNAME='MAX_VAL' $
            , SCR_XSIZE=80, VALUE=max_disp)
  WID_MIN_VAL = Widget_Slider(WID_BASE_3, MAX=max_disp-(max_disp-min_disp)/79. $
            , MIN=min_disp, SCROLL=1, UNAME='MIN_VAL' $
            , SCR_XSIZE=80, VALUE=min_disp)

;=========== Control buttons =================================
  WID_BASE_2c = Widget_Base(WID_BASE_0b, /ROW)
  WID_LABEL_6 = Widget_Label(WID_BASE_2c, UNAME='WID_LABEL_6'  $
      , /ALIGN_LEFT, VALUE='RMS pixels', XOFFSET=10)

  WID_TEXT_6 = Widget_Text(WID_BASE_2c, UNAME='RMS_PIX'  $
      , VALUE=[ '  -1.000' ], XSIZE=8, YSIZE=1)

  WID_LABEL_7 = Widget_Label(WID_BASE_2c, UNAME='WID_LABEL_7'  $
      , /ALIGN_LEFT, VALUE='RMS m/s', XOFFSET=10)

  WID_TEXT_7 = Widget_Text(WID_BASE_2c, UNAME='RMS_M_S'  $
      , VALUE=[ '   -1.' ], XSIZE=6, YSIZE=1)

  WID_LABEL_8 = Widget_Label(WID_BASE_2c, UNAME='WID_LABEL_8'  $
      , /ALIGN_LEFT, VALUE='Number of lines used', XOFFSET=10)

  WID_TEXT_8 = Widget_Text(WID_BASE_2c, UNAME='NLINES'  $
      , VALUE=[ '   ' ], XSIZE=6, YSIZE=1)

  WID_BASE_2a = Widget_Base(WID_BASE_2c, UNAME='WID_BASE_2a', /ROW, /NONEXCLUSIVE)
  WID_BUTTON_7 = Widget_Button(WID_BASE_2a, UNAME='BAD_ORDER'  $
      , /ALIGN_CENTER, TOOLTIP='Include/exclude selected spectral order', VALUE='Good order')

  WID_BASE_2b = Widget_Base(WID_BASE_0b, /ROW)

  WID_BUTTON_0 = Widget_Button(WID_BASE_2b, UNAME='1D_SOLUTION'  $
      , /ALIGN_CENTER, TOOLTIP='Construct wavelength solution for'+ $
      ' selected order', VALUE='1D solution')

  WID_BUTTON_2 = Widget_Button(WID_BASE_2b, UNAME='ALIGN'  $
      , /ALIGN_CENTER, TOOLTIP='Align 2D model with the actual data', VALUE='Align')

  WID_BUTTON_1 = Widget_Button(WID_BASE_2b, UNAME='2D_SOLUTION'  $
      , /ALIGN_CENTER, TOOLTIP='Construct 2D wavelength solution using'+ $
      ' all processed orders', VALUE='2D solution')

  WID_BUTTON_3 = Widget_Button(WID_BASE_2b, UNAME='AUTO_ID'  $
      , /ALIGN_CENTER, TOOLTIP='Automatically identify lines in all'+ $
      ' orders', VALUE='AutoID')

  WID_BUTTON_4 = Widget_Button(WID_BASE_2b, UNAME='MOD_ECH'  $
      , /ALIGN_CENTER, TOOLTIP='Add 2D wavelength solution to the ECH'+ $
      ' file with the comparison spectrum', VALUE='Add 2D to ECH')

  WID_BUTTON_5 = Widget_Button(WID_BASE_2b, UNAME='REPORT'  $
      , /ALIGN_CENTER, TOOLTIP='Produces graphics report in a PS file'+ $
      ' including spectral resolution and the PSF', VALUE='Report')

  WID_BUTTON_6 = Widget_Button(WID_BASE_2b, UNAME='EXIT'  $
      , /ALIGN_CENTER, VALUE='EXIT')

;=========== Draw window =================================
  WID_DRAW_0 = Widget_Draw(WID_BASE_0, UNAME='WID_DRAW_0' $
      , SCR_XSIZE=((base_x_size-20)<1024)-10, SCR_YSIZE=(base_y_size-130)<800 $
      , XSIZE=((base_x_size-20)<1024)-10, YSIZE=(base_y_size-130)<800  $
      , RETAIN=1, /BUTTON_EVENTS, /MOTION_EVENTS, /TRACKING_EVENTS)

  device,decomposed=0
  loadct,0,/SILENT


  Widget_Control, /REALIZE, WID_BASE_0
  Widget_Control, WID_DRAW_0, get_value=graph0


  !x.margin=[8,5]
  display, cspec, min=min_disp,max=max_disp       ; Plot spectrum
  npixel = n_elements(cspec[*,0])
  norders = n_elements(cspec[0,*])
  wl_range = dblarr(2, norders)-1
  solution_1D=dblarr(7,norders)
  bad_order_flag=0

  line={csline, wlc:-1.d0,  $ ; computed wavelength
                wll:-1.d0,  $ ; laboratory wavelength
                posc:-1.d0, $ ; computed position
                posm:-1.d0, $ ; measured position
                xfirst:0,   $ ; starting pixel
                xlast:0,    $ ; ending pixel
                approx:'G', $ ; approximation type
                width:0.d0, $ ; Line width in pixels
                flag:0,     $ ; flag status
                height:0.d0,$ ; Line strength
                order:0}      ; Absolute spectral order

  bad_order=bytarr(norders)   ; 1 - bad order, 0 - order is fine

;  if(keyword_set(save_file)) then f=save_file else f=file_prefix+'_2D.sav'
;  if(file_test(f)) then restore,f

  widget_control, WID_BUTTON_0,sensitive=0
  widget_control, WID_BUTTON_5,sensitive=0
  if(n_elements(cs_lines) gt 0) then begin
    for iord = 0, norders-1 do begin
      i=where(cs_lines.order eq iord,ni)
      if(ni gt 0) then oplot, !x.crange, [iord,iord], col=160 $
      else             oplot, !x.crange, [iord,iord], col= 80
    endfor
    if(min(cs_lines.order) eq max(cs_lines.order)) then begin
      widget_control, WID_BUTTON_1,sensitive=0
      widget_control, WID_BUTTON_3,sensitive=0
      widget_control, WID_BUTTON_4,sensitive=0
    endif
    widget_control, WID_TEXT_8, set_value=strtrim(n_elements(cs_lines),2)
  endif else begin
    for iord = 0, norders-1 do oplot, !x.crange, [iord,iord], col= 80
    widget_control, WID_BUTTON_1,sensitive=0
    widget_control, WID_BUTTON_3,sensitive=0
    widget_control, WID_BUTTON_4,sensitive=0
    widget_control, WID_TEXT_8, set_value='  0'
  endelse
  widget_control, WID_BUTTON_7,sensitive=0
  if(obase gt 0) then widget_control,WID_TEXT_4,set_value=strtrim(obase,2)
  if(oincr ne 0) then widget_control,WID_TEXT_5,set_value=strtrim(oincr,2)

  XManager, 'WID_BASE_0', WID_BASE_0;, /NO_BLOCK

end


;=========================================================================================
; Empty stub procedure used for autoloading.
;
pro wavecal, CS_ech_file, save=save_file, POLARIZ=polariz $
           , REVERSE=flip_wl, UPSIDE_DOWN=flip_orders    $
           , GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_

  if(not keyword_set(CS_ech_file)) then begin
    print,'Usage: wavecal,<ech_file_with_ThAr_spectrum>[, SAVE=<2D_solution>[, /REVERSE'
    print,'                                            [, /UPSIDE_DOWN[, POLARIZ=ord]]]]'
    print,'       SAVE    an existing 2D solution. If the solution was the ech file was already'
    print,'               obtained earlier, it will be found automatically with no need for'
    print,'               using this option.'
    print,'      /REVERSE flips pixel numbering in the suggested model given by "save"'
    print,'      /UPSIDE_DOWN flips order numbering in the suggested model given by "save"'
    print,'       POLARIZ is set when a polarization spectrum (2 spectra per sp. order) is used'
    print,'               If ord is 1 than the lower spectrum in each pair will be used, while'
    print,'               2 refers to the upper spectrum. The output files will receive the'
    print,'               corresponding suffix.' 
    return
  endif
  if(not file_test(CS_ech_file)) then begin
    print,'ECH file "'+CS_ech_file+'" was not found'
    return
  endif
  prefix = strmid(CS_ech_file, 0, strpos(CS_ech_file, '.', /REVERSE_SEARCH))
  rdech, CS_struct, CS_ech_file,/NOCONT
  CS_struct.spec = CS_struct.spec > 0.
  if(max(CS_struct.spec) lt 1.) then $
    CS_struct.spec=(CS_struct.spec-min(CS_struct.spec))*10000

;  if(keyword_set(flip_wl)) then begin
;    nord=n_elements(CS_struct.spec[0,*])
;    npix=n_elements(CS_struct.spec[*,0])
;    ind=npix-1L-lindgen(npix)
;    for iord=0,nord-1 do CS_struct.spec[*,iord]=CS_struct.spec[ind,iord]
;  endif

  loadct,0
  !p.region=0
  !p.position=0
  WID_BASE_0, CS_struct, CS_ech_file, prefix, save_file, POLARIZ=polariz $
            , REVERSE=flip_wl, UPSIDE_DOWN=flip_orders, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_
end

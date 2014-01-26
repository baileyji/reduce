#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits

##First call
#hamspec, im, head, orders, xwd, sp $
#    , sig=sig, colrange=colrange, SF_SMOOTH=10., OSAMPLE=12 $
#    , swath_width=200, MASK=mask, FILENAME=file_basename(imfile);, $
#
##which calls
#getspec, im, orc, xwd, spec, sunc $         ;extract spectrum
#    , gain, dark, readn $
#    , colrange=colrange, sf_smooth=sf_smooth, sp_smooth=sp_smooth $
#    , osample=osample, SWATH_WIDTH=swath_width $
#    , MASK=mask, PLOT=iplot, ORDER_RANGE=order_range $
#    , TELLURIC=telluric, FILENAME=filename $
#    , SLIT_TILT_FILE=slit_tilt_file $
#    , WING_SMOOTH_FACTOR=wing_smooth_factor $
#    , SWATH_BOUNDARIES=swath_boundaries
#
##or for thar
#getspec,im,orc,xwd,spec,0,gain,dark,readn $
#    ,colrange=colrange, noopt=True

def getspec( im, orc, xwd, spec, sunc, gain, dark ,readn  $
           , NOOPT=False, colrange=None, sf_smooth=6
           , SP_SMOOTH=sp_smooth, OSAMPLE=osample, SWATH_WIDTH=swath_width $
           , MASK=mask, PLOT=iplot, process_orders=None,
           , TELLURIC=telluric, FILENAME=filename $
           , use_2d_slit=slit_tilt_file $
           , SWATH_BOUNDARIES=swath_boundaries $
           , WING_SMOOTH_FACTOR=wing_smooth_factor):
"""
im - the image
orc - nord x n_ord_coeff array
xwd - nord x 2 array
colrange - nord x 2 array first col, last col +1 to extract
process_orders - list of order numbers to process, None for all. default None. 1 indexed
use_2d_slit - True or file containing info, passed to load_2d_slit_info Doesnt work
"""

nrow, ncol = im.shape
nord,ncoef = orc.shape


#Default order range span the whole image
if not colrange:
    colrange=np.tile(np.array([0,ncol+1]),(nord,1))

    
ix = np.arange(ncol)                        #column indicies
spec = np.zeros((ncol,nord),dtype=np.float)				#init spectrum
sunc = spec.copy()                          #init uncertainties
if noopt:
     sunc +=1.0                             #why? thar unc = 1

#;GETARC needs order location coefficients (orc) on both sides of arc swath to
#;  be extracted. In order to extract the lowest and highest orders specified
#;  in orc, we need to extend orc one extra order on each end. We shall do so
#;  by linearly extrapolating the last two orc on each end.
#;Extend orc on the low end. Check that requested swath lies on image.
if nord >1:
    ixx=ix[slice(*colrange[0])]
    orclo = 2*orc[0] - orc[1]               #extrapolate orc
    if xwd[0,0] > 1.5:                      #Extraction width in pixels
      coeff=orc[0]
      coeff[0]=coeff[0]-xwd[0,0]
    else:                                   #Fraction extraction width
      coeff = 0.5 * ((2+xwd[0,0])*orc[0] - xwd[0,0]*orc[1])
    yoff, = np.where(np.polyval(coeff, ixx) < 0)		#pixels off low edge
else:
    yoff=[]
    orclo = np.zeros(ncoeff,dtype=np.float)

if len(yoff) > 0:				#check if on image
#   GETARC will reference im(j) where j<0. These array elements do not exist.
    print 'GETSPEC: Top order off image in columns [{},{}]'.format(
        yoff[0], yoff[noff-1])

#Extend orc on the high end. Check that requested swath lies on image.
if nord > 1:
    ixx=ix[slice(*colrange[nord-1])]
    orchi = 2*orc[nord-1] - orc[nord-2]      #extrapolate orc
    if xwd[nord-1,1] > 1.5:                  #Extraction width in pixels
      coeff=orc[nord-1]
      coeff[0]=coeff[0]+xwd[nord-1,1]
    else:                                    #Fraction extraction width
      coeff = 0.5 * ((2+xwd[nord-1,1])*orc[nord-1] - xwd[nord-1,1]*orc[nord-2])
    yoff, = np.where(np.polyval(coeff, ixx) > nrow-1)	#pixels off high edge
  else:
    yoff = []
    orchi = np.zeros(ncoeff,dtype=np.float)
    orchi[0] = nrow-1

if len(yoff) > 0:
# GETARC will reference im(j) where j > ncol*nrow-1. These elements do not exist.
    print 'GETSPEC: Bottom order off image in columns [{},{}].'.format(
            yoff[0],yoff[noff-1])

#Define an order set (orcend) extended one extra order on either end.
orcend = np.zeros((nord+2,ncoef),dtype=np.float)
orcend[1:-2,:]=orc
orcend[0] = orclo
orcend[nord+1] = orchi


#Orders to process
if not process_orders:
    process_orders=range(1,nord+1)


#Doing optimal extraction?
if noopt:                           #no optimal extraction
    if xwd.max() <= 1:              #fractional orders
        print ('GETSPEC: Interpreting extraction width '
                'as fraction of each order to mash.')
    else:                           #fixed pixel widths
       print ('GETSPEC: Interpreting extraction width '
              'as number of pixels to mash.')
else:                               #optimal extraction
    print,'GETSPEC: Using optimal extraction to produce spectrum.'

if not noopt:                       #Optimal extraction of the spectrum

    if use_2d_slit:                 #2D slit function is requested
        pass
        #get slit tilt info (shear_x, ..., ???)

    swath_boundaries=intarr(2,1)

    for onum in process_orders:                    #loop thru orders
    
        #Columns for extraction
        cole0, cole1 = colrange[onum-1]         #first & last column to extract
        ncole =cole1 -cole0 +1                  #number of columns to extract
    
        #Status reporting
        if nord <= 10:
            print 'GETSPEC: extracting relative order {} of {}.'.format(
                    onum, nord)
        elif (onum-1) % 5 == 0:
            print'GETSPEC: extracting relative order {}-{} of {}.'.format(
                    onum, min(nord, onum+4), nord)


        #First & last column to extract
        cole0, cole1 = colrange[onum-1]

        #Column indices to extract
        ixx = ix[slice(*colrange[onum-1])]

        #Compute extraction region (ycen, ycenn, ymin, ymax)
        ycen = np.polyval(orcend[onum], ix)			#row at order center

        ycenn = ycen[slice(*colrange[onum-1])]

        if xwd[onum-1,0] > 1.5:                #Extraction width in pixels
            ymin = ycenn - xwd[onum-1,0]
        else:                                  #Fractional extraction width
            ymin = (ycenn - xwd[onum-1,0] *
                    (ycenn - np.polyval(orcend[onum-1],ixx)))   #trough below

        ymin = floor(ymin)

        if ymin.min() < 1:
            ymin = ymin - ymin.min() + 1        #help for orders at edge

        if xwd[onum-1,1] > 1.5:                #Extraction width in pixels
            ymax = ycenn + xwd[onum-1,1]
        else:                                  #Fractional extraction width
            ymax = (ycenn + xwd[onum-1,1] *
                    (np.polyval(orcend[onum-1],ixx)-ycenn))   #trough above

        ymax = ceil(ymax)

        if ymax.max() >= nrow -1:
            ymax = ymax - ymax.max() + nrow - 2         #help for orders at edge


#        #Define a fixed height area containing one spectral order
#      y_lower_lim = fix(min(ycen(cole0:cole1)-ymin))              ;Pixels below center line
#      y_upper_lim = fix(min(ymax-ycen(cole0:cole1)))              ;Pixels above center line
#      if(y_lower_lim <= 0 or y_upper_lim le 0) then begin
#        message, 'Bad order trace or something, skiping'+string(onum),/info
#        continue
#      endif
#
#
#      unc=replicate(-2.,ncol)
#      mkslitf, im,
#             , ycen, y_lower_lim, y_upper_lim $
#             , yslitf, slitf, binc, onum $
#             , x_left_lim=cole0, x_right_lim=cole1-1, PLOT=iplot $ x_right_lim is column number not a range spec
#             , LAMBDA_SF=sf_smooth, LAMBDA_SP=sp_smooth, blz=blz $
#             , OSAMPLE=osample, SWATH_WIDTH=swath_width, MASK=mask $
#             , GAIN=gain, READN=readn, /NO_SCATTER $
#             , TELLURIC=telluric, FILENAME=filename $
#             , slit_tilt=shear_x $
#             , SWATH_BOUNDARIES=sw_b $
#             , WING_SMOOTH_FACTOR=wing_smooth_factor $
#             , UNCERTAINTY=unc
#
#
#      spec(cole0:cole1,onum-1) = blz(cole0:cole1)	 ;save in output array
#      getwindow,4
#      axplot(spec[onum-1,cole0:cole1])
#      ,title='Spectrum, order '+strtrim(onum,2)
#      swath_boundaries=[[swath_boundaries],[sw_b]]
#      if(max(unc) lt -1) then begin
#        sunc(cole0:cole1,onum-1) = sqrt(abs(blz(cole0:cole1)))
#      endif else begin
#        sunc(cole0:cole1,onum-1) = unc(cole0:cole1)
#      endelse
#    endfor
#    swath_boundaries=swath_boundaries(*,1:*)
#
#  endif else begin					;else: getarc style
#    for onum=1,nord do begin				;loop thru orders
#      cole0  = colrange[0,onum-1]                  ;First column to extract
#      cole1 = colrange[1,onum-1]                  ;Last column to extract
#      awid = xwd(0,onum-1) + xwd(1,onum-1)
#      getarc,im,orcend,onum,awid,arc,pix $
#      	    ,x_left_lim=cole0,x_right_lim=cole1	;extract counts/pixel
#;      if subback then arc = arc - back(*,onum-1)	        ;subtract backgd/pixel
#      spec(cole0:cole1,onum-1) = arc * pix		;store total counts
#      sunc(cole0:cole1,onum-1) = sqrt(abs(arc*pix*gain + dark + $
#                            pix*readn^2)) / gain          ;estimate uncertainty
#    endfor
#  endelse
#
#return


def mkslitf(im, varim, $
           , ycen, y_lower_lim, y_upper_lim, yslitf $
           , slitf, bincen, ord_num, PLOT=iplot $
           , X_LIM=None $
           , lambda_sf=1.0, lambda_sp=0, swath_width=None $
           , BLZ=blz, osample=10, MASK=mask, gain=1.0, readn=0.0 $
           , no_scatter=True, telluric=0, FILENAME=filename $
           , slit_tilt=None, normalize=False, NORM_IMAGE=im_norm  $
           , ORDER_IMAGE=im_ordr, THRESHOLD=threshold $
           , return_swath_bounds=return_swath_bounds $
           , WING_SMOOTH_FACTOR=wing_smooth_factor $
           , UNCERTAINTY=sunc):
"""
x_lim - 2 element iterable, first & last column to extract, def: 1st & last col
gain - ccd gain in e-/dn
readnoise - in e-
varim - image of the variance in e-, if set im must be in e-
shearx - variable describing the 2d slit
swath_width - ??? def: estimate it
ycen - ncol array of order centers
y_lower_lim/y_upper_lim - scalar, min number of pixels above/below the the center line
x_left_lim/x_right_lim - the first and last column indices for extraction
return_swath_bounds - T/F if set return swath boundaries
osample - slitf pixels / real pixel default = 10
lambda_sf - ??? def =1.0
lambda_sp - ??? def =0
no_scatter - whether to skip worrying about scattered light. Only False is supported
telluric - If set to specific integer larger than 1 is used as the offset from 
    the order center line. The sky is then estimated by computing median signal 
    between this offset and the upper/lower limit of the extraction window.
    Support for this is unfinished in this version as M2FS can not use it.
normalize - def False. Set this to ... ????
"""
#;Determines slit function along echelle order
#;Input:
#; im (array(ncol,nrow)) image containing echelle spectrum
#; back (vector(ncol) vector containing background under the current order
#; ymin (vector(ncol)) row numbers along bottom of region to map
#; ycen (vector(ncol)) row numbers of zero point for slit function
#; ymax (vector(ncol)) row numbers along top of region to map
#;Output:
#; yslitf (vector(nslitf)) subpixel row offsets for slitf
#; sflit (array(nslitf,nbin)) subpixel slit functions
#; bincen (vector(nbin)) column of bin centers
#;


    ret={}  #The return dictionary

    nrow, ncol = im.shape                #number of rows & columns

    #Columns for extraction
    if not x_lim:
        x_left_lim,x_right_lim= 0,ncol-1        #first & last column to extract
    else:
        x_left_lim,x_right_lim=x_lim

    #TODO: Sort gain & readnoise out
    CCD_gain  = gain  #CCD gain e- per 1 ADU
    CCD_readn = readn #Readout noise
    noise=CCD_readn/CCD_gain

    #TODO: Understand the slit tilt
    if slit_tilt:
        shear=-np.polyval(slit_tilt[*,ord_num-1], np.arange(ncol,dtype=np.float))
        delta_x=ceil(np.abs(max(y_lower_lim, y_upper_lim) *np.abs(shear).max()))

    #Vet the mask
    msk = 0
    imask = False   #Have mask to use?
    if mask and mask.shape != im.shape:
        raise ValueError('MKSLITF: Your mask is '
                         'defined but does not match your image')
    elif mask:
        imask=True
        #TODO: Get rid of this flag


    #Compute a swath width if not specified, as ncol/avg_ncol_btwn_crossings
    min_swaths=3 #legacy min and max from NP, don't know the reasoning
    max_swaths=20
    if not swath_width:
        #TODO this isn't true np.unique meaning is different
        _, i = np.unique(ycen.round(), return_index=True)  #Points of row crossing
        ni = len(i)                          #This is how many times this order
                                             #crosses to the next row
        if ni > 1:                              #Curved order crosses rows
            i = (i[1:ni]-i[0:ni-1]).sum()/(ni-1) #average columns between row crossings
            #number of swaths along the order
            nbin = min(max( round(ncol/i) / 3 , min_swaths) , max_swaths)
        else:                            #Perfectly aligned orders
            nbin = max(ncol/400 , 3)     #Still follow the changes in PSF

        #Adjust for the true order length
        nbin = nbin * (x_right_lim - x_left_lim + 1) / ncol
                 
    else:                                  #If swath width is preset
        nbin = max(round((x_right_lim - x_left_lim + 1.) / swath_width),1)

    #integer trace centers
    yc = ycen.astype(np.int)

    #Find rows that contain spectrum.
    imin = yc - y_lower_lim                       #bottom row of order
    imax = yc + y_upper_lim                       #top row of order

    #Calculate boundaries of distinct slit function regions.
    nhalf=2*nbin-1
    if nbin > 1:
        #TODO: could just do this by tweaking the params to linspace
        # later
        ibound=np.linspace(x_left_lim, x_right_lim, nbin+1) #boundaries of bins
        ibeg = ceil(ibound[:-1])               #beginning of each bin
        iend = np.concatenate([ibeg[1:] - 1, x_right_lim])
        ibeg_half = (ibeg[1:] + ibeg[:-1]) / 2 #boundaries of overalapping steps
        iend_half = (iend[1:] + iend[:-1]) / 2

        ibeg_half = np.concatenate((ibeg_half, ibeg))
        ibeg_half.sort()

        iend_half = np.concatenate((iend_half, iend))
        iend_half.sort()
    else:
        ibeg_half=np.array([x_left_lim ,x_right_lim])
        iend_half=np.array([x_right_lim,x_right_lim])

    #Store swath bounaries for return if requested
    if return_swath_bounds:
        ret['swath_boundaries']=np.concatenate((ibeg_half[np.newaxis,:-1],
                                                iend_half[np.newaxis,:-1])).T

    #Swath centers
    bincen = 0.5*(ibeg_half + iend_half)

    #Initialize more stuff
    irow = np.arange(nrow)                      #indices of all rows
    nysf = y_upper_lim + y_lower_lim + 1        #subpixel range required
    yslitf0 = -y_lower_lim                      #minimum value for yslitf
    yslitf1 =  y_upper_lim                      #maximum value for yslitf

    blz = np.zeros(ncol, dtype=np.float)        #???



    #Perform slit decomposition within each swath
    #   stepping through the order with half swath width. Spectra for each
    #   decomposition are combined with linear weights.
    for ihalf in range(nhalf):               #loop thru swaths
        ib = ibeg_half[ihalf]                         #left column
        if slit_tilt and ihalf > 0:
            ib = ib - delta_x
    
        ie = iend_half[ihalf]                         #right column
        if slit_tilt and ihalf < nhalf-1:
            ie = ie + delta_x

        nc = ie - ib + 1                        #number of columns in the swath

        #Load slit function data into vectors.
        j0 = lonarr(nc)                     #???
        j1 = lonarr(nc)         #???
        sf = np.zeros((nc,nysf), dtype=np.float)    #The slit function ???
        sfpnt = fltarr(nc * nysf )
        ysfpnt = fltarr(nc * nysf )
        if imask:
            msk = bytarr(nc,nysf)

        #Populate sf, sfpnt, ysfpnt, j0, j1, & msk
        #I think this section is taking the curvy spectrum out of im
        # and putting into into a rectangular array just big enough to hold the
        # data. Its also doing this for the mask data.
        # finally it is generating a flattened version and a indexing arrays
        # all of which are used in an unclear way as of yet!
        #TODO: understand this
        for j in range(nc):                  #loop thru columns in region
            icen = yc[ib+j]                  #row closest to peak
            k0 = icen + yslitf0              #lowest row to consider
            k1 = icen + yslitf1              #highest row to consider
            j0[j] = j*nysf                   #begining of storage area
            j1[j] = j0[j] + k1 - k0          #new logic (works at edge)
            if no_scatter:
                ssf = im[ib+j,k0:k1+1]
            else:
                raise Exception('Scatter not implemented in python version')
                #dy_scatter = (dindgen(nysf) + k0 - yscatter_below[ib+j]) $
                #             / (yscatter_above(ib+j) - yscatter_below(ib+j))
                #scatter = (scatter_above(ib+j)-scatter_below(ib+j))* $ #Interpolate background
                #           dy_scatter+scatter_below(ib+j)
                #ssf = (im(ib+j,k0:k1) - scatter)>0                     #subtract scatter

            sf[j] = ssf
            sfpnt[j0[j]:j1[j]+1] = ssf
            ysfpnt[j0[j]:j1[j]+1] = irow[k0:k1+1] - ycen[ib+j]
            if imask:
                msk[j] = mask[ib+j,k0:k1+1]

        #Do telluric processing of some sort, see the telluric param
        if telluric:
            print 'WARNING: telluric param untested. Not useful for M2FS'
            if telluric < nysf/2:
                tel_lim=telluric
            else:
                tel_lim=min(5, nysf/3)
            tel=sf.sum(axis=1)          #sum over nysf
            itel=np.arange(nysf)
            itel=itel[abs(itel-nysf/2) >= tel_lim]
            tel=sf[:,itel]

            ipdb.set_trace()
            #TODO: verifiy identical functionality and broadcasting
            #sc=np.zeros(nc)
            #for itel in range(nc-1):
            #    sc[itel]=np.median(tel[itel,:]) # in IDL: sc[itel]=median(reform(tel[itel,*]))

            sc=np.median(tel,axis=1)
            sf -= sc # numpy should broadcast this #np.ones(yslitf1-yslitf0+1)

            sfpnt = sf.T.reshape(len(sfpnt))

        #Compute the slit function
        jgood=0
        if slit_tilt:
            raise Exception('slit_func_2d not implemented')
            #sp, sfsm = slit_func_2d(sf, ycen[ib:ie+1]-yc[ib:ie+1] $
            #       ,delta_x,shear,NOISE=noise $
            #       ,OVERSAMPLE=osample,IM_OUT=sfbin,LAMBDA_SF=lambda_sf $
            #       ,LAMBDA_SP=lambda_sp,USE_COL=jgood,BAD=jbad,MASK=msk $
            #       ,WING_SMOOTH_FACTOR=wing_smooth_factor $
            #       ,UNCERTAINTY=unc,DEBUG=iplot)
        else:
            sp, sfsm, sfbin, jbad, unc = slit_func(sf,ycen[ib:ie+1]-yc[ib:ie+1],
                noise=noise, osample=osample, im_out=True, lambda_sf=lambda_sf,
                lambda_sp=lambda_sp, use_col=None, mask=msk,
                wing_smooth_factor=wing_smooth_factor, debug=iplot)



        weight = np.ones(nc,np.float) #weight vector for combining the overalaping part
        if ihalf > 0:
            weight[:nc/2 + 1] = np.arange(nc/2 + 1,dtype=np.float)/(nc/2)
        oweight = 1.0 - weight

        #In case we do FF normalization replace the original image by the
        # ratio of sf/sfbin where number of counts is larger than threshold
        # and with 1 elsewhere
        scale=1
        if normalize:
        #TODO: This is where I left off
            if imask:
                ii = where(sfbin gt threshold/gain and msk eq 1B, nii) $
            else:
                ii = where(sfbin gt threshold/gain, nii)
            sss = replicate(1.,nc,nysf)                 #1     - by default
            ddd = fltarr(nc,nysf)                       #save orders, helps analysing ghosts
            if nii > 0:
                sss[ii] = sf(ii)/sfbin(ii);ratio - when high signal
            ddd = sfbin
            if ihalf > 0:

                #Here is the layout to understand the lines below
                #
                #        1st swath    3rd swath    5th swath      ...
                #     /============|============|============|============|============|
                #
                #               2nd swath    4th swath    6th swath
                #            |------------|------------|------------|------------|
                #            |.....|
                #            overlap
                #
                #            +     ******* 1
                #             +   *
                #              + *
                #               *            weights (+) previous swath, (*) current swath
                #              * +
                #             *   +
                #            *     +++++++ 0


                overlap=iend_half[ihalf-1] - ibeg_half[ihalf] + 1
                #Code was commented in the IDL version
                #scale=median((blz[ibeg_half[ihalf]:iend_half[ihalf-1]]>1)/(sp[0:overlap-1]>1))
                #print,scale
                #stop
                scale=1
                if nii > 0:
                    sss[ii]=sss[ii]/scale
                sp=sp*scale
            else:
                nc_old = nc                         #previous swath width
                sss_old = fltarr(nc,nysf)           #dummy normalized ff from previous swath
                ddd_old = fltarr(nc,nysf)           #dummy order shape ff from previous swath
                overlap = nc_old + 1

            #This loop is complicated because swaths do overlap
            # to ensure continuity of the spectrum.
            if ihalf == 0:
                ncc=ibeg_half[1] - ibeg_half[0]
            else:
                ncc=overlap

            for j=0, ncc-1:                               #loop thru columns in region
                icen = yc(ib+j)                           #row closest to peak
                k0 = icen + yslitf0                       #lowest row to consider
                k1 = icen + yslitf1                       #highest row to consider
                j0(j) = j * nysf                          #begining of storage area
                j1(j) = j0(j) + k1 - k0                   #new logic (works at edge)
                jj = nc_old - ncc + j                     #column number in the previous swath
                im_norm(ib+j,k0:k1) = sss_old(jj,*) * oweight(j) $#replace the original
                                    + sss    (j ,*) *  weight(j)
                im_ordr(ib+j,k0:k1) = ddd_old(jj,*) * oweight(j) $#replace the original
                                    + ddd    (j ,*) *  weight(j)

            if ihalf == 2*nbin-2:                           #finish the very last swath
                for j=ncc,nc-1:                             #loop thru columns in region
                    icen = yc(ib+j)                         #row closest to peak
                    k0 = icen + yslitf0                     #lowest row to consider
                    k1 = icen + yslitf1                     #highest row to consider
                    j0(j) = j * nysf                        #begining of storage area
                    j1(j) = j0(j) + k1 - k0                 #new logic (works at edge)
                    im_norm(ib+j,k0:k1) = sss(j ,*)         #replace the original
                    im_ordr(ib+j,k0:k1) = ddd(j ,*)         #replace the original

            nc_old  = nc
            sss_old = sss
            ddd_old = ddd
        #End of FF normalization part


    nbad=n_elements(jbad)
    if(nbad eq 1) then if(jbad eq -1) then nbad=0

;    sp=total(sf,2)

    for j=0,nc-1 do begin          ; Normalize sfpnt for the plot
      sfpnt(j0(j):j1(j)) = sfpnt(j0(j):j1(j)) / (sp(j)>1.)
      sfbin(j,*) = sfbin(j,*) / (sp(j)>1.)
    endfor
    
    if(slit_tilt) then begin
      if(ihalf gt 0 and ihalf lt 2*nbin-2) then begin   ; Copy the extracted spectrum
        blz(ib+delta_x:ie-delta_x) = blz(ib+delta_x:ie-delta_x) * oweight(delta_x:nc-delta_x-1) $
                                   + sp(delta_x:nc-delta_x-1)   *  weight(delta_x:nc-delta_x-1)
      endif else if(ihalf eq 0) then begin
        blz(ib        :ie-delta_x) = sp(0:nc-delta_x-1)
      endif else if(ihalf eq 2*nbin-2) then begin
        blz(ib+delta_x:ie)         = blz(ib+delta_x:ie) * oweight(delta_x:nc-1) $
                                   + sp(delta_x:nc-1)   *  weight(delta_x:nc-1)
      endif
    endif else begin
      blz(ib:ie) = blz(ib:ie) * oweight + sp * weight
    endelse
    if(arg_present(sunc)) then sunc(ib:ie) = sunc(ib:ie) * oweight + unc * weight

    if(ihalf eq 0) then begin
      nslitf = n_elements(sfsm)
      yslitf = yslitf0 + (findgen(nslitf)-0.5)/osample-1.5 ;final subpixel scale
      slitf = fltarr(nslitf, 2*nbin-1)                   ;init final slit function
    endif

    sfsm2 = reform(transpose(sfbin), nc * nysf )
    j = sort(ysfpnt)
;    if(keyword_set(istop)) then stop
;    if(nbad gt 0) then for k=0L,nbad-1 do jbad(k) = (where(j eq jbad(k)))(0)
; The line below does exactly the same as the line above. I have no clue why. NP
    if(nbad gt 0) then jbad = (sort(j))(jbad)
    ysfpnt = ysfpnt(j)
    sfpnt = sfpnt(j)
    sfsm2 = sfsm2(j)
;    j=uniq(ysfpnt)
;    jj=where(yslitf gt min(ysfpnt(j)) and yslitf le max(ysfpnt(j)))
;    sfsm(jj) = spl_interp(ysfpnt(j), sfsm2(j), spl_init(ysfpnt(j), $
;                          sfsm2(j)), yslitf(jj))
    slitf(*,ihalf) = sfsm/total(sfsm)*osample         ;save slit function

;    plot,blz&oplot,indgen(ie-ib+1)+ib,sp,col=2
;    stop
;    sf0=sf & sfbin0=sfbin & sfsm0=sfsm

;Plot diagnostics.
    if keyword_set(iplot) then begin
      colors
      getwindow,0,/show
      !p.multi=[0,2,2]

      pscale = mean(sp)
      sfplot=smooth(sfsm*pscale,(osample-1)>1)/scale
      plot, ysfpnt, sfpnt*pscale, ps=3, /xsty, ysty=3, yr=minmax(sfpnt*pscale) $
          , /NODATA, title='Order '+strtrim(ord_num,2)+ $
          ', Columns '+strtrim(ib,2)+' through '+strtrim(ie,2)
      oplot, ysfpnt, sfpnt*pscale, ps=3
      if(nbad gt 0) then oplot, [ysfpnt(jbad)], [sfpnt(jbad)*pscale], ps=1, co=3
      ymiddle=0.5*(!y.crange(0)+!y.crange(1))
      oplot, yslitf, sfplot, thick=2, col=2
      oplot, !x.crange, [0.,0.], co=4

      plot, ysfpnt, sfpnt*pscale, ps=3, /xsty, ysty=3, yr=minmax(sfsm*pscale) $
          , /NODATA, title='Order '+strtrim(ord_num,2)+ $
          ', Columns '+strtrim(ib,2)+' through '+strtrim(ie,2)
      oplot, ysfpnt, sfpnt*pscale, ps=3
      if(nbad gt 0) then oplot, [ysfpnt(jbad)], [sfpnt(jbad)*pscale], ps=1, co=3
      ymiddle=0.5*(!y.crange(0)+!y.crange(1))
      oplot, yslitf, sfplot, thick=2, col=2
      oplot, !x.crange, [0.,0.], co=4

      plot, ysfpnt,(sfpnt - sfsm2)*pscale, ps=3, xs=1, title='Data - Fit'
      if(nbad gt 0) then oplot, [ysfpnt(jbad)], [(sfpnt(jbad) - sfsm2(jbad))*pscale] $
                              , ps=1, col=3
      oplot, !x.crange,  [0.,0.], co=4
      if no_scatter:
          poffset = 0
      elif (scatter_above and scatter_below):
          poffset = (scatter_below[ib:ie+1] + scatter_above[ib:ie+1]).mean()*0.5
      else:
          raise Exception('Scatter logic not implemented fully')

      oplot, yslitf,  sqrt(((sfsm*pscale/scale + poffset + CCD_readn^2)/CCD_gain)>0),col=2,thick=2
      oplot, yslitf, -sqrt(((sfsm*pscale/scale + poffset + CCD_readn^2)/CCD_gain)>0),col=2,thick=2
      oplot,  (!x.crange(1)-!x.crange(0))*[0.03,0.1]+!x.crange(0) $
           ,  (!y.crange(1)-!y.crange(0))*[0.06,0.06]+!y.crange(0), col=2,thick=2
      xyouts, (!x.crange(1)-!x.crange(0))*0.1+!x.crange(0) $
            , (!y.crange(1)-!y.crange(0))*[0.05,0.05]+!y.crange(0) $
            , '1!7r!3 + scatter + read noise'

      yr = sqrt(max(sfsm*pscale + poffset + CCD_readn^2)/CCD_gain)*2.
      yr = [-yr,yr]
      plot, ysfpnt,(sfpnt - sfsm2)*pscale, ps=3, xs=1, ys=3 $
;          , yr= minmax(sfpnt(jgood) - sfsm2(jgood))*pscale, title='Fit - Data'
          , yr = yr, title='Data - Fit'
      if(nbad gt 0) then oplot, [ysfpnt(jbad)], [(sfpnt(jbad) - sfsm2(jbad))*pscale] $
                              , ps=1, col=3
      oplot, !x.crange,  [0.,0.], co=4
      oplot, yslitf,  sqrt(((sfsm*pscale/scale + poffset + CCD_readn^2)/CCD_gain)>0),col=2,thick=2
      oplot, yslitf, -sqrt(((sfsm*pscale/scale + poffset + CCD_readn^2)/CCD_gain)>0),col=2,thick=2
      oplot,  (!x.crange(1)-!x.crange(0))*[0.03,0.1]+!x.crange(0) $
           ,  (!y.crange(1)-!y.crange(0))*[0.06,0.06]+!y.crange(0), col=2,thick=2
      xyouts, (!x.crange(1)-!x.crange(0))*0.1+!x.crange(0) $
            , (!y.crange(1)-!y.crange(0))*[0.05,0.05]+!y.crange(0) $
            , '1!7r!3 + scatter + read noise'

      if(keyword_set(filename)) then xyouts,0.5,0.5,filename,align=0.5,/NORM
      empty
      wshow,0
      !p.multi=[0,0,0]
    endif
  endfor

; Combine fitting errors and the Poisson shot noise
  sunc=sqrt((sunc+blz)>1.)



def make_o(y, yycen, i, nrow, omega, sf, osample):
    """Snippet of code just to clean things up"""
    #Shift the positions of the subpixels by the amount the trace is into the central row
    #for instace if the trace center fell at row 12.5, then we'd be shifting by 0.5
    yy=y+yycen[i]     # Offset SF,
    
    
    #omega=fltarr(n,nrow)                    ; weights for ovsersampling
    #                                        ; omega(k,j) is how much
    #                                        ; point k in oversampled SF
    #                                        ; contributes to the image pixel j
    #for j=0,nrow-1 do begin
    #    ind=where(yy gt j and yy lt j+1, nind)
    #    if(nind gt 0) then begin
    #        omega(ind,j)=weight           ; All sub-pixels that fall to image
    #        i1=ind(0)                           ; pixel j have weight 1/osample
    #        i2=ind(nind-1)
    #        omega(i1,j)=yy(i1)-j                ; On the two ends we may have
    #        omega(i2+1,j)=j+1-yy(i2)            ; boundary crossing so weights
    #                                            ; could be less than 1./osample
    #    endif
    #endfor
    #omega=reform(sf#omega)
    
    
    i1,i2=np.where(yy >= 0 and yy < nrow)[[0,-1]]
    
    
    #TODO: Understand. This really doesn't make sense. yy is some sort of offset in units of pixels
    # and omega is a weight
    #I think omega are the weights each sub pixel contributes to the pixel
    omega[0]=yy[i1]
    
    #sf is the oversampled slit function with (i think) some padding on each end
    #I think ssf is supposed to be the oversampled slit function without padding and
    # corrected to be centered on the trace (by shifting y by ycen
    # and reshaped to be osampleXnrow, but I don't really understand what the shift
    # is doing or why it is being done. I suspect this might be because I don't quite
    # understand the output of bandsol
    ssf=sf[i1:i2+1].reshape((osample,nrow)) #reform(sf(i1:i2),osample,nrow)
    
    #o seems to be the slit function in pixel space
    o=np.dot(omega, ssf) #reform(ssf##omega)
    
    #o needs more processing to make, not sure what is happening ???
    yyy=nrow-yy[i2]
    o[:nrow-1]+= ssf[0,1:nrow]*yyy
    o[nrow-1] += sf[i2+1]*yyy
    
    return o



def generate_Akl_Bl(y, yycen, osample, imm, sp, iter,
                    olind, bklind, msk, nrow, lambda_sf, sf, n,
                    wing_smooth_factor=wing_smooth_factor):

    Akl=np.zeros((n,2*osample+1),dtype=float)           # Initialize matrix
    Bl=np.zeros(n,dtype=float)                          # and RHS
    omega=np.zeros(osample+1, dtype=np.float) + weight  # Replicate constant weights
    for i in range(0,ncol):                    # Fill up matrix and RHS
        yy=y+yycen[i]                          # Offset SF

        i1,i2=np.where(yy >= 0 and yy < nrow)[[0,-1]] # Weights are the same within pixel except
                                                      # for the first and the last subpixels.

        omega[[0,osample]]=yy[i1],1.0-yy[i2]    # Fix the first and the last subpixel, here
                                                # the weight is split between the two subpixels

        o=np.outer(omega, omega)
        o[osample,osample]+=o[0,0]

        oo=o[olind] #TODO, is this the right subscripting

        bkl=dblarr(n,2*osample+1)       # Band-diagonal part that will contain omega#omega
        for l in range(0,nrow):
#            for m=osample,2L*osample:     # Explicit and slow way of filling bkl
#                mm=m-osample
#                bkl(l*osample+i1:l*osample+i1+osample-mm,m)=o(oind(0:osample-mm)+mm)
            bkl[l*osample+i1+bklind]=oo*msk[l,i]
            
        for l in range(1,nrow):
            bkl[l*osample+i1, osample]=o[osample,osample]*msk[l,i]
        
        bkl[nrow*osample+i1,osample]=msk[nrow-1,i]*omega[osample]**2 #blk in IDL subscripting order

        for m in range(0,osample):
            bkl[osample-m:n,m]=bkl[:n-osample+m,2*osample-m] #blk in IDL subscripting order
        
        Akl+=sp[i]**2 * bkl

        o=np.zeros(n,dtype=np.float)
        for l in range(0,nrow):
            o[l*osample+i1:l*osample+i1+osample+1]= imm[l,i]*weight*msk[l,i]

        for l in range(1,nrow):
            o[l*osample+i1]= (imm[l-1,i]*omega[osample]*msk[l-1,i] +
                              imm[l,i]*omega[0]*msk[l,i])
        
        o[i1]=imm[0,i]*omega[0]*msk[0,i]
        o[nrow*osample+i1]=imm[nrow-1,i]*omega[osample]*msk[nrow-1,i]
        
        Bl+=sp[i]*o

    #TODO: Deal with moving this into the function
    if debug:
        time1=time()

    lambda_ = lambda__sf * Akl[:,osample].sum()/n
    if wing_smooth_factor and iter > 1:
        #lambda_=lambda_*(1.+wing_smooth_factor*(2.d0*dindgen(n)/(n-1)-1.d0)^2)
        lambda_*=(1.0 + wing_smooth_factor/sf.clip(1e-5))
    else:
        lambda_=np.zeros(n, dtype=np.float) + lambda_


    # 1st order Tikhonov regularization (minimum 1st derivatives)
    # Add the following 3-diagonal matrix * lambda_:
    #  1 -1  0  0  0  0
    # -1  2 -1  0  0  0
    #  0 -1  2 -1  0  0
    #  0  0 -1  2 -1  0
    #      .  .  .
    #
    #Akl(  0,osample)=Akl(  0,osample)+lambda_   # +lambda_ to the upper-left element
    #Akl(n-1,osample)=Akl(n-1,osample)+lambda_   # and to the lower-right
    #Akl(1L:n-2L,osample)=Akl(1L:n-2L,osample)+2.*lambda_    # +2*lambda_ to the rest of the main diagonal
    #Akl(0L:n-2L,osample+1L)=Akl(0L:n-2L,osample+1L)-lambda_ # -lambda_ to the upper sub-diagonal
    #Akl(1L:n-1L,osample-1L)=Akl(1L:n-1L,osample-1L)-lambda_ # -lambda_ to the lower sub-diagonal

    Akl[  0,osample]=Akl[  0,osample]+lambda_[0]# +lambda_ to the upper-left element
    Akl[n-1,osample]=Akl[n-1,osample]+lambda_[n-1]# and to the lower-right
    Akl[1:n-1,osample]=Akl[1:n-1,osample]+2.0*lambda_[1:n-1]# +2*lambda_ to the rest of the main diagonal
    Akl[:n-1,osample+1]=Akl[:n-1,osample+1]-lambda_[:n-1]# -lambda_ to the upper sub-diagonal
    Akl[1:n,osample-1]=Akl[1:n,osample-1]-lambda_[1:n]# -lambda_ to the lower sub-diagonal

    # 2nd order Tikhonov regularization (minimum 2nd derivative)
    # Add the following 5-diagonal matrix * lambda_:
    #  1 -2  1  0  0  0
    # -2  5 -4  1  0  0
    #  1 -4  6 -4  1  0
    #  0  1 -4  6 -4  1
    #      .  .  .
    #
    #lambda_=0.1*lambda_
    #Akl(  0,osample)=Akl(  0,osample)+1.*lambda_ # Main diagonal
    #Akl(n-1,osample)=Akl(n-1,osample)+1.*lambda_
    #Akl(  1,osample)=Akl(  1,osample)+5.*lambda_
    #Akl(n-2,osample)=Akl(n-2,osample)+5.*lambda_
    #Akl(2L:n-3L,osample)=Akl(2L:n-3L,osample)+6.*lambda_
    #Akl(0L,osample+1L)=Akl(0L,osample+1L)-2.*lambda_ # upper sub-diagonal
    #Akl(n-2L,osample+1L)=Akl(n-2L,osample+1L)-2.*lambda_
    #Akl(1L:n-3L,osample+1L)=Akl(1L:n-3L,osample+1L)-4.*lambda_
    #Akl(1L,osample-1L)=Akl(1L,osample-1L)-2.*lambda_ # lower sub-diagonal
    #Akl(n-1L,osample-1L)=Akl(n-1L,osample-1L)-2.*lambda_
    #Akl(2L:n-2L,osample-1L)=Akl(2L:n-2L,osample-1L)-4.*lambda_

    return Akl, Bl




#def bandsol(aa,rr):
#
#  n=n_elements(rr)
#  sz=size(aa)
#  if(n eq 0 or sz(0) ne 2 or sz(1) ne n) then begin
#    print,'bandsol solve a sparse system of linear equations with band-diagonal matrix.'
#    print,'Band is assumed to be symmetrix relative to the main diaginal. Usage:'
#    print,'res=bandsol(a,r[,/DOUBLE])'
#    print,'where a is 2D array [n,m] where n - is the number of equations and m'
#    print,'        is the width of the band (3 for tri-diagonal system),'
#    print,'        m is always an odd number. The main diagonal should be in a(*,m/2)'
#    print,'        The first lower subdiagonal should be in a(1:n-1,m-2-1), the first'
#    print,'        upper subdiagonal is in a(0:n-2,m/2+1) etc. For example:'
#    print,'               / 0 0 X X X \'
#    print,'               | 0 X X X X |'
#    print,'               | X X X X X |'
#    print,'               | X X X X X |'
#    print,'           A = | X X X X X |'
#    print,'               | X X X X X |'
#    print,'               | X X X X X |'
#    print,'               | X X X X 0 |'
#    print,'               \ X X X 0 0 /'
#    print,'      r is the array of RHS of size n.'
#    print,'      /DOUBLE forces the calculations to double precision.'
#    return,0
#  endif
#  nd=sz(2)
#
#  if(keyword_set(dbl)) then begin
#    a=double(aa)
#    r=double(rr)
#  endif else begin
#    a=float(aa)
#    r=double(rr)
#  endelse
#
#  for i=0L,n-2L do begin
#    r(i)=r(i)/a(i,nd/2) & a(i,*)=a(i,*)/a(i,nd/2)
#    for j=1L,(nd/2)<(n-i-1L) do begin
#      r(i+j)=r(i+j)-r(i)*a(i+j,nd/2-j)
#      a(i+j,0:nd-j-1)=a(i+j,0:nd-j-1)-a(i,j:nd-1)*a(i+j,nd/2-j)
#    endfor
#  endfor
#
#  r(n-1L)=r(n-1L)/a(n-1L,nd/2)
#  for i=n-1L,1L,-1L do begin
#    for j=1L,(nd/2)<i do begin
#      r(i-j)=r(i-j)-r(i)*a(i-j,nd/2+j)
#    endfor
#    r(i-1L)=r(i-1L)/a(i-1L,nd/2)
#  endfor
#  r(0L)=r(0L)/a(0L,nd/2)
#
#  return r

def slit_func(im, ycen, osample=1, lambda_sf=0.1,
                lambda_sp=0, im_out=False, use_col=None,
                mask=None, noise=0, model_only=False,
                debug=False,wing_smooth_factor=wing_smooth_factor,
                preset_slit_func=None):
    """
    im - mxn array, don't use nan's, trace is centered on central row
    ycen - float fraction al position of trace within central row as function of column
    mask - none or array of same shape as image 1= good pixels
    osample - the amount of oversampling ??? def: 1 = none
    lambda_sf - ??? default 0.1 I think this is a typo
    lambda_sp - default 0. ???
    use_col - sets columns of im to use Def: None-> use all columns, 
        set to a 1d numpy array for indexing im's columns if needed
    model_only - def false. not used in REDUCE package, so implementation might be buggy
        If desired set to a tuple of (sp, sf). sp & sf should have the same
        format as would normally be returned by this routine. It doesn't make any sense 
        to set this without also setting im_out
    preset_slit_func - Set to the slit function for preset value, the array is copied
        and should have length equal to nrows in im
    debuf - bool do debugging
    bad - set to return jbad
    im_out - set to return the reconstructed swath image
    wing_smooth_factor - ??? controls the amount of wing smoothing

    Returns:
    sp - the spectrum
    sf - the slit function
    optional (returned if im_out is set):
    im_out - the reconstructed swath image
    unc - some sort of uncertainty ???
    bad - the indicies of bad pixels in im, perhaps ??? (formerly jbad)

    """
    import scipy.signal.medfilt
    #Vet oversample
    if osample < 1:
        raise ValueError('oversample must be >= 1.')

    #Setup the mask
    if not mask:
        mmsk=np.ones_like(im, dtype=np.bool)
    elif mask.shape!=im.shape:
        raise ValueError('SLIT_FUNC: Mask must have the same size as the image')
    else:
        mmsk=mask.copy()

    oind=np.arange(osample+1)*(osample+2) #oind ???
    weight=1.0/osample   # nominal weight each subpixel contributes to pixel


    #Pick out the columns to use for the procedure, which is???
    # use_col is called as 0 everytime slit_func is called by mkslitf
    # so normal flow is through option 2
    if use_col:
        imm=im[:,use_col]
        yycen=ycen[use_col]
        msk=mmsk[:,use_col]
    else:
        use_col=np.arange(im.shape[0])
        imm=im.copy()
        yycen=ycen.copy()
        msk=mmsk.copy()

    #Get rows and columns
    nrow,ncol=imm.shape

    #Number of oversampled points incluing a full padding pixel on either side
    n=(nrow + 1)*osample+1
    #y positions of the subpixels (includes padding)
    y=np.arange(n, dtype=np.float)/osample - 1

    #fraction of unmasked pixels
    norm=float(len(msk))/msk.sum() #JIB added float, seems appropriate

    #Build bklind & olind
    # bklind is ???
    # olind is ???
    bklind=np.arange(osample+1) + n*osample
    olind=oind[0:osample+1]
    for m in range(osample+1, 2*osample+1):
        mm=m-osample
        bklind=np.concatenate((bklind, np.arange(osample+1-mm)+n*m))
        olind=np.concatenate((olind,oind[:osample-mm+1]+mm))

    #If we are only doing the model (whatever that means)
    if model_only:
        print ('Warning model only never used in REDUCE. '
               'Implementation untested.')
        sp, sf=model_only
        dev=np.sqrt(sp.mean()+noise**2)
    else:

        #Compute initial guess for spectrum and slit function
        if not preset_slit_func:           #Slit function is unknown
            sf=(imm*msk).sum(axis=1)        #The slit function
            if osample > 2 and len(sf) > 5:
                sf=scipy.signal.medfilt(sf,5)   #Smooth it a bit
          
            if imm.sum(axis=0).mean() < 1e3:  #in case of low S/N robust guess for sf
                sf=np.exp(-((np.arange(nrow,dtype=np.float)-nrow/2.)/
                            (nrow/4.))**2)
            sf=sf/sf.sum()
          
            #Initial guess for the spectrum
            sp=norm * ((imm*msk)*np.tile(sf,(ncol,1))).sum(axis=0)
        else:                     #Slit function is given
            #Initial guess for the spectrum
            sp=(imm*msk).sum(axis=0)/max(msk.sum(axis=0), 1)
            sf=preset_slit_func.copy()

        #Finish of the guess at the spectrum¸
        if osample > 2:
            sp=scipy.signal.medfilt(sp,5)
        sp*=(imm*msk).sum()/sp.sum()


        #rms of differences from model
        spsf_outer_prod=np.outer(sp,sf)
        dev=np.sqrt((msk* (imm-spsf_outer_prod)**2).sum()/msk.sum())

        #mask all with outliers > 3*rms
        msk[np.abs(imm-spsf_outer_prod) > 3.0*dev]=False


        #Iterate and solve for the spectrum and slit function
        iter=0  #Iteration counter
        max_iter=18
        while True:
        
            iter+=1

            # Solve for the slit function if not preset
            if not preset_slit_func:
  
                if debug:
                    time0=time()
                
                Bl, Akl=generate_Akl_Bl(y, yycen, osample, imm, sp, olind, iter,
                                        bklind, msk, nrow, lambda_sf, sf, n,
                                        wing_smooth_factor=wing_smooth_factor)
                sf=bandsol(Akl,Bl)
                #TODO: enforece sf >= 0
                sf=sf/sf.sum()*osample

                if debug:
                    time2=time()

            sp_old=sp.copy()
            r=sp.copy()

            omega=np.zeros(osample, dtype=np.float) + weight

            if debug:
                sssf=np.zeros((nrow,ncol), dtype=np.float)   # ????

            #Evaluate the new spectrum
            dev_new=0.0
            for i in range(ncol):
                
                #Compute o. I think this is some thing like the slit function in
                # in pixel space normalized somehow
                o=make_o(y, yycen, i, nrow, omega, sf, osample)
                
                #Compute the ??? at this column
                # r seems to be the profile-weighted mashed values at this column, masked
                # sp seems to be the mashed profile squared, masked
                # neither seem to take into account normalization effects due to masked pixels
                r[i] = (imm[:,i]*msk[:,i]*o).sum() #IDL: ((imm(i,*)*msk(i,*))#o)  = 1xnrow # nrow
                sp[i]=(o**2 * msk[:,i]).sum()
                
                if sp[i]==0.0: #ok
                    sp[i]=(o**2).sum()

                if debug:
                    sssf[:,i]=r[i]/sp[i]*o

                #Locate and mask outliers
                if iter > 1:
                    norm=r[i]/sp[i]
                    
                    j=np.abs(imm[:,i]-norm*o) > 6.0*dev #j=where(abs((imm[i,*]-norm*o)) > 6.*dev,nj,COMPLEMENT=b)

                    msk[j,i]=False
                    #load back the original mask for anything that isn't an
                    # outlier this iteration
                    msk[~j,i]=mmsk[~j,use_col[i]].copy()
                    
                    #Add this columns squared error from the data ???
                    dev_new+=(msk[:,i] * (imm[:,i] - norm*o)**2).sum()
                    #Different variant, commented out in original code
                    #dev_new+=(msk[:,i] * (imm[:,i]/norm-o)**2).sum()

            ### Finished evaluating the new spectrum ###

            #Do spectral smoothing if requested
            if lambda_sp:
                lambda_ = lambda_sp*sp.sum()/ncol
                
                #in IDL it was:
                #a =[0.,replicate(-lambda_,ncol-1)] #ldiag
                #b =[lambda_+1.,replicate(2.*lambda_+1.,ncol-2),lambda_+1.] #diag
                #c =[replicate(-lambda_,ncol-1),0.]  #udiag
                #sp=trisol(a,b,c,r/sp)
                
                sp=scipy.linalg.solve_banded(
                        (nrows-3, ncols-3),
                        np.array([
                          [-lambda_]*(ncol-1)+[0.0],
                          [lambda_ + 1]+[2*lambda_+1]*(ncol-2)+[lambda_+1],
                          [0.0]+[-lambda_]*(ncol-1)
                        ]), r/sp,
                        overwrite_ab=True,
                        overwrite_b=True,
                        check_finite=True)
            else:
                sp=r/sp


            #Update the RMS error
            if iter gt 1:
                dev=np.sqrt(noise**2 + dev_new/msk.sum())

            #yy  =imm/(sp#replicate(1.d,nrow))
            #dev_new=total(msk*(yy-sssf)^2)
            if iter > 1:
                dev=np.sqrt(dev_new/msk.sum())

            if debug:
                time3=time()
                
                figure(1)
                
                x=np.arange(im.shape[0])+ 0.5 + 0.5/osample
                
                sssf/=np.tile(sp,(nrow,1)) #IDL sssf/(replicate(1.d,nrow)#sp)

                xx=x-ycen[0]
                yy=im[:,0]/sp[0]
                
                #TODO: sort this out
                for i in range(1,ncol):
                    xx=[xx,x-ycen[i]]
                    yy=[yy,im[:,i]/sp[i]]
                
                ii=yy.argsort()[:11:-1]
                
                print ('Not shure why this matters... '
                      'ii/nrow={}, ii mod nrow ={}'.format(ii/nrow, ii % nrow))
                      
                ii=xx.argsort()
                xx=xx[ii]
                yy=yy[ii]
                
                #TODO: init Subplot
                
                plt.plot(xx,sssf(ii),'k',label='sssf[ii]')

                plt.axhline(sf.mean(), linestyle='-.')
                plt.plot(xx,yy,'k.',label='yy')
                bad,=np.where(~(transpose(msk))[ii])
                if len(bad):
                    plt.plot(xx[bad],yy[bad],'r+',label='bad')
                plt.plot(xx,yy-sssf[ii]+sf.mean(),'g.',label='yy-sssf[ii]+sf.mean()')
                plt.xlabel('Shifted Row')
                plt.ylabel('Units Here')
                plt.title('Slit Function. Iteration={}'.format(iter))

                #TODO: init Subplot
                plt.plot(sp)
                plt.xlabel('Swath index')
                plt.ylabel('Counts')
                plt.title('Swath Spectrum')
                print ('Maximum change in the spectrum'
                       ': {}. dev: {}. something {}.'.format(
                            max(np.abs(sp-sp_old)/sp.max()),
                            dev,
                            np.sqrt(((yy-yyy)**2).sum()/len(yy))))
                print 'Runtimes: {} {} {}'.format(time1-time0,
                                                  time2-time1,
                                                  time3-time2)
                time.sleep(.3)

            if iter >= max_iter or (np.abs(sp-sp_old)/sp.max()).max() <= 1e-5:
                break

        ### End of iteration loop ###
    ### End of if model_only else ###


    if not im_out:
        im_out=None
        jbad=None
    else:
        jbad=[]
        if use_col:
            sp=interpol(sp, use_col, np.arange(ncol,dtype=np.float))

        im_out=np.zeros((nrow,ncol), dtype=np.float)
        unc=np.zeros(ncol, dtype=np.float)

        omega=np.zeros(osample, dtype=np.float) + weight

        for i in range(0,ncol):                #Evaluate the new spectrum
        
            #Compute o. I think this is some thing like the slit function in
            # in pixel space normalized somehow
            o=make_o(y, yycen, i, nrow, omega, sf, osample)
   
            #Compute column in model image
            im_out[:,i]=sp[i]*o
   
            #Compute uncertainties
   
            #Pixels in the column < 5 * dev from the model
            good=np.abs(im[:,i]-sp[i]*o) < 5*dev
            j=np.where(good)

            if len(j) < nrow:          #Bad pixels in column i
                #TODO: Depending on how the caller uses this we may want
                # to change the way the bad subscripts are returned
                jbad.append(i*nrow + np.where(~good))

            if len(j) > 2:  #At least 2 good pixels in column
                #TODO: Understand the uncertainty computation
                
                # sum of squares of difference from image and model @ col i
                #  good pix only
                ss=((im[j,i]-sp[i]*o[j])**2).sum()
                #JIB: shouldn't this be
                #ss=(msk[j,i]*(im[j,i]-sp[i]*o[j])**2).sum()

                # ngood*(ngood-2) * variance of o[good]
                xx=((o[j]-o[j].mean())**2).sum()*(len(j)-2)
                unc[i]=ss/xx
            else:
                #TODO: Why 0? I think this is a flag to indicate that the value is bad
                unc[i]=0.0

        if debug:
            plt.figure(2)
            plt.imshow(im_out)
            plt.title('Reconstructed swath')
            plt.xaxis('swath column')
            plt.yaxis('trace row')
            print 'Click with mouse to terminate.'
            if plt.waitforbuttonpress()==False:
                raise Exception('Aborted by user')

        if jbad:
            jbad=np.concatenate(jbad)
        else:
            jbad=np.array([],dtype=np.int)

    j,=np.where(~msk.any(axis=0))
    if len(j) > 0:
        print 'Warning: columns {} is/are fully masked.'.format(j)
        #TODO: get the numpy function for interpol
#      sp[j]=interpol(sp[msk.any(axis=0)], np.where(msk.any(axis=0)), j)

    if im_out != None:
        return sp, sf, im_out, jbad, unc
    else:
        return sp, sf

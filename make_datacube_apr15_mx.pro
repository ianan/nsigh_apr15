pro make_datacube_apr15_mx,mosid=mosid,chuid=chuid,dirout=dirout,maindir=maindir

  ;  Make the maps per each mos pointing
  ;         Do it per FPMA & FPMB (though only FPMA at the moment)
  ;         Filter out so only CHU12, 3, 123 or 23, grade=0 and > 2keV
  ;         Save all these out without gaussian smoothing.
  ;
  ;         For non-IGH use need to change
  ;         maindir - where April data is kept
  ;         dirout - where to put the maps
  ;
  ;         mosid - 10 to 26 for MOS1, 30 to 46 for MOS2
  ;         chuid - 1 to 4 (1=12,2=3,3=123,4=23)

  ; based on make_datacube_Nov14_O4.pro

  ; IGH 17-Jun-2015 - Modified version for one mosiac tile in April 2015
  ; IGH 29-Jun-2015 - Added dirout, maindir options
  ; IGH 08-Apr-2015 - Removes "hot" pixels and corrects bug in engs before histogram
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ; For saving it out
  ; shifts, smoothing, combinations can be done in another file
  if (n_elements(dirout) ne 1) then dirout='~/data/ns_data/obs4_maps/'
  if (n_elements(mosid) ne 1) then mosid=40;20
  mnm=string(mosid,format='(i2)')
  if (n_elements(chuid) ne 1) then chuid=4

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (n_elements(maindir) ne 1) then maindir='~/data/ns_data/obs4_bg/'

  subdir='201100'+mnm+'_Sol_15119_MOS0'+mnm
  ddname='201100'+mnm+'001'

  chm=[5,9,14,13]
  chmn='CHU'+['012','003','123','023']
  chumask=chm[chuid-1]
  chunam=chmn[chuid-1]

  ; Assuming centre of fov is Sun Centre
  ; Shift from this applied later via xshf, yshf
  xc=0.0
  yc=0.0

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Get the *_cl_sunpost.evt files

  cla_file=maindir+subdir+'/'+ddname+'/event_cl/nu'+ddname+'A06_cl_sunpos.evt'
  evta = mrdfits(cla_file, 1,evtah)

  clb_file=maindir+subdir+'/'+ddname+'/event_cl/nu'+ddname+'B06_cl_sunpos.evt'
  evtb = mrdfits(clb_file, 1,evtbh)


  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Before doing anything else need to filter out the "bad"/"hot" pixels in FPMA
  ; these are manually found in test_evtlst_apr15.pro and mostly different from the Nov15 ones

  use = bytarr(n_elements(evta)) + 1
  thisdet = where(evta.det_id eq 2)
  badones = where(evta[thisdet].rawx eq 16 and evta[thisdet].rawy eq 5, nbad)
  if nbad gt 0 then use[thisdet[badones]]=0
  badones = where(evta[thisdet].rawx eq 8 and evta[thisdet].rawy eq 22, nbad)
  if nbad gt 0 then use[thisdet[badones]]=0

  thisdet = where(evta.det_id eq 3)
  badones = where(evta[thisdet].rawx eq 16 and evta[thisdet].rawy eq 11, nbad)
  if nbad gt 0 then use[thisdet[badones]]=0
  badones = where(evta[thisdet].rawx eq 0 and evta[thisdet].rawy eq 16, nbad)
  if nbad gt 0 then use[thisdet[badones]]=0
  badones = where(evta[thisdet].rawx eq 10 and evta[thisdet].rawy eq 1, nbad)
  if nbad gt 0 then use[thisdet[badones]]=0
  badones = where(evta[thisdet].rawx eq 10 and evta[thisdet].rawy eq 7, nbad)
  if nbad gt 0 then use[thisdet[badones]]=0
  
  evta=evta[where(use)]


  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Doing it in just one CHU so can filter out other CHU as well
  ; Based on https://github.com/NuSTAR/nustar_solar/blob/master/solar_mosaic_20141211/solar_mosaic_hk.pro
  chufile = file_search(maindir+subdir+'/'+ddname+'/hk/', '*chu123.fits')
  for chunum= 1, 3 do begin
    chu = mrdfits(chufile, chunum)
    maxres = 20 ;; [arcsec] maximum solution residual
    qind=1 ; From KKM code...
    if chunum eq 1 then begin
      mask = (chu.valid EQ 1 AND $          ;; Valid solution from CHU
        chu.residual LT maxres AND $  ;; CHU solution has low residuals
        chu.starsfail LT chu.objects AND $ ;; Tracking enough objects
        chu.(qind)(3) NE 1)*chunum^2       ;; Not the "default" solution
    endif else begin
      mask += (chu.valid EQ 1 AND $            ;; Valid solution from CHU
        chu.residual LT maxres AND $    ;; CHU solution has low residuals
        chu.starsfail LT chu.objects AND $ ;; Tracking enough objects
        chu.(qind)(3) NE 1)*chunum^2       ;; Not the "default" solution
    endelse
  endfor

  ; make time binning of chus to evt data
  achu_comb = round(interpol(mask, chu.time, evta.time))
  bchu_comb = round(interpol(mask, chu.time, evtb.time))

  ; filter out bad CHUs and not single pixel hits (not GRADE 0)
  ida2=where(achu_comb eq chumask and evta.grade eq 0)
  evta=evta[ida2]
  a_engs=1.6+0.04*evta.pi

  idb2=where(bchu_comb eq chumask and evtb.grade eq 0)
  evtb=evtb[idb2]
  b_engs=1.6+0.04*evtb.pi

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Get the livetime info
  hka_file=maindir+subdir+'/'+ddname+'/hk/nu'+ddname+'A_fpm.hk'
  hka = mrdfits(hka_file, 1, hkahdr)
  hkatims=anytim(hka.time+anytim('01-Jan-2010'))

  t1a=anytim(min(evta.time)+anytim('01-Jan-2010'),/yoh,/trunc)
  t2a=anytim(max(evta.time)+anytim('01-Jan-2010'),/yoh,/trunc)

  lvida=where(hkatims ge anytim(t1a) and hkatims lt anytim(t2a))
  lvtcora=mean(hka[lvida].livetime)

  hkb_file=maindir+subdir+'/'+ddname+'/hk/nu'+ddname+'B_fpm.hk'
  hkb = mrdfits(hkb_file, 1, hkbhdr)
  hkbtims=anytim(hkb.time+anytim('01-Jan-2010'))

  t1b=anytim(min(evtb.time)+anytim('01-Jan-2010'),/yoh,/trunc)
  t2b=anytim(max(evtb.time)+anytim('01-Jan-2010'),/yoh,/trunc)

  lvidb=where(hkbtims ge anytim(t1b) and hkbtims lt anytim(t2b))
  lvtcorb=mean(hkb[lvidb].livetime)


  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Setup the pixel and binning sizes
  ; Get the same values if using evtah or evtbh
  ttype = where(stregex(evtah, "TTYPE", /boolean))
  xt = where(stregex(evtah[ttype], 'X', /boolean))
  xpos = (strsplit( (strsplit(evtah[ttype[max(xt)]], ' ', /extract))[0], 'E', /extract))[1]
  npix = sxpar(evtah, 'TLMAX'+xpos)
  pix_size = abs(sxpar(evtah,'TCDLT'+xpos))

  centerx = round(xc / pix_size) + npix * 0.5
  centery = round(yc / pix_size) + npix * 0.5
  im_size = 1037. / pix_size
  im_width = round(im_size * 2.)

  ; energy binning for output spectrum/datacube
  maxe=10.
  mine=1.6
  ebin=0.04
  eid=mine+ebin*indgen(1+(maxe-mine)/ebin)
  nengs=n_elements(eid)-1
  ims=intarr(nengs,2*im_width,2*im_width)

  ; Need to shift the pointing ? Maybe up to arcmin off AIA/RHESSI/XRT
  xshf=0. ; arcsec/pix_size
  yshf=0.
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Do FPMA

  engs=1.6+0.04*evta.pi
  evtx=evta.x-xshf
  evty=evta.y-yshf

  for i=0, nengs-1 do begin
    ide=where(engs ge eid[i] and engs lt eid[i+1],nid)
    print,eid[i]
    if (nid gt 1) then begin
      ; this data still  has the x opposite direction to standard solar coords

      ;      plot,evtx[ide],evty[ide],psym=1,xr=[1300,2000],yr=[1300,2000]

      ;    pixinds = (npix - evta[ide].x) + evta[ide].y * npix
      pixinds = evta[ide].x + evta[ide].y * npix
      im_hist = histogram(pixinds, min = 0, max = npix*npix-1, binsize = 1)
      im = reform(im_hist, npix, npix)
      im= im[(centerx-im_width):(centerx+im_width-1), (centery-im_width):(centery+im_width-1)]
      ims[i,*,*] = im
    endif
    if (nid eq 1) then begin
      pixinds = (npix - evtx[ide]) + evty[ide] * npix
      im_hist = lonarr(npix*npix)
      im_hist[pixinds]=1
      im = reform(im_hist, npix, npix)
      im= im[(centerx-im_width):(centerx+im_width-1), (centery-im_width):(centery+im_width-1)]
      ims[i,*,*] = im
    endif
  endfor

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Save all this out
  npp=n_elements(ims[0,0,*])
;  ; Save out whole image array
  subx=[0.0,npp-1]
  suby=[0.0,npp-1]
;  ; Save out just the MOS40 region !!! Wont work for other MOS tiles
  subx=[550,1150]
  suby=[700,1300]
  ims=ims[*,subx[0]:subx[1],suby[0]:suby[1]]

  pxs=pix_size
  x0=xc-npp*0.5*pxs+pxs*subx[0]
  y0=yc-npp*0.5*pxs+pxs*suby[0]

  newxc=x0+0.5*n_elements(ims[0,*,0])*pxs
  newyc=y0+0.5*n_elements(ims[0,0,*])*pxs

  dc_out={specs:ims,eng_edgs:eid,t1:t1a,t2:t2a,lvtfac:float(lvtcora),$
    x0:x0,y0:y0,dx:pxs,xc:newxc,yc:newyc,id:'MOS'+mnm+': FPMA'}

  save,file=dirout+'MOS'+mnm+'_EngCube_GO_FPMA.dat',dc_out

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Do FPMB

  ims=intarr(nengs,1688,1688)
  engs=1.6+0.04*evtb.pi

  evtx=evtb.x-xshf
  evty=evtb.y-yshf

  for i=0, nengs-1 do begin
    ide=where(engs ge eid[i] and engs lt eid[i+1],nid)
    print,eid[i]
    if (nid gt 1) then begin
      ;    pixinds = (npix - evta[ide].x) + evta[ide].y * npix
      pixinds = evta[ide].x + evta[ide].y * npix
      im_hist = histogram(pixinds, min = 0, max = npix*npix-1, binsize = 1)
      im = reform(im_hist, npix, npix)
      im= im[(centerx-im_width):(centerx+im_width-1), (centery-im_width):(centery+im_width-1)]
      ims[i,*,*] = im
    endif
    if (nid eq 1) then begin
      pixinds = (npix - evtx[ide]) + evty[ide] * npix
      im_hist = lonarr(npix*npix)
      im_hist[pixinds]=1
      im = reform(im_hist, npix, npix)
      im= im[(centerx-im_width):(centerx+im_width-1), (centery-im_width):(centery+im_width-1)]
      ims[i,*,*] = im
    endif
  endfor

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Save all this out
  npp=n_elements(ims[0,0,*])
  ims=ims[*,subx[0]:subx[1],suby[0]:suby[1]]

  pxs=pix_size
  x0=xc-npp*0.5*pxs+pxs*subx[0]
  y0=yc-npp*0.5*pxs+pxs*suby[0]

  newxc=x0+0.5*n_elements(ims[0,*,0])*pxs
  newyc=y0+0.5*n_elements(ims[0,0,*])*pxs

  dc_out={specs:ims,eng_edgs:eid,t1:t1b,t2:t2b,lvtfac:float(lvtcorb),$
    x0:x0,y0:y0,dx:pxs,xc:newxc,yc:newyc,id:'MOS'+mnm+': FPMB'}

  save,file=dirout+'MOS'+mnm+'_EngCube_GO_FPMB.dat',dc_out

end

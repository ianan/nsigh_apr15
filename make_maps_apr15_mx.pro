pro make_maps_apr15_mx,mosid=mosid,chuid=chuid,dirout=dirout

  ; For saving it out
  ; shifts, smoothing, combinations can be done in another file
  dirout='~/data/ns_data/obs4_maps/'
  if (n_elements(mosid) ne 1) then mosid=40;20
  mnm=string(mosid,format='(i2)')
  if (n_elements(chuid) ne 1) then chuid=4

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (n_elements(dirout) ne 1) then dirout='~/data/ns_data/obs4_maps/'
  subdir='201100'+mnm+'_Sol_15119_MOS0'+mnm
  ddname='201100'+mnm+'001'

  chm=[5,9,14,13]
  chmn='CHU'+['012','003','123','023']
  chumask=chm[chuid-1]
  chunam=chmn[chuid-1]

  xr=[-700,700]
  yr=[-300,1100]

  restore,file=dirout+'MOS'+mnm+'_EngCube_GO_FPMA.dat'
  ; For some reason my map2fits.pro complains about no L0,B0 or Rsun if not provided
  ang = pb0r(dc_out.t1,/arcsec,l0=l0)

  ima2=total(dc_out.specs[where(dc_out.eng_edgs ge 2),*,*],1)/dc_out.lvtfac
  ma2=make_map(ima2,dx=dc_out.dx,dy=dc_out.dx,xc=dc_out.xc,yc=dc_out.yc,$
    time=dc_out.t1,dur=anytim(dc_out.t2)-anytim(dc_out.t1),id='FPMA G0 >2 keV',l0=l0,b0=ang[1],rsun=ang[2])
  sub_map,ma2,sma2,xrange=xr,yrange=yr
  map2fits,sma2,'maps/MOS'+mnm+'_EngCube_GO_FPMA_E2.fits'

  ima24=total(dc_out.specs[where(dc_out.eng_edgs ge 2 and dc_out.eng_edgs le 4),*,*],1)/dc_out.lvtfac
  ma24=make_map(ima24,dx=dc_out.dx,dy=dc_out.dx,xc=dc_out.xc,yc=dc_out.yc,$
    time=dc_out.t1,dur=anytim(dc_out.t2)-anytim(dc_out.t1),id='FPMA G0 2-4 keV',l0=l0,b0=ang[1],rsun=ang[2])
  sub_map,ma24,sma24,xrange=xr,yrange=yr
  map2fits,sma24,'maps/MOS'+mnm+'_EngCube_GO_FPMA_E24.fits'

  ima46=total(dc_out.specs[where(dc_out.eng_edgs ge 4 and dc_out.eng_edgs le 6),*,*],1)/dc_out.lvtfac
  ma46=make_map(ima46,dx=dc_out.dx,dy=dc_out.dx,xc=dc_out.xc,yc=dc_out.yc,$
    time=dc_out.t1,dur=anytim(dc_out.t2)-anytim(dc_out.t1),id='FPMA G0 4-6 keV',l0=l0,b0=ang[1],rsun=ang[2])
  sub_map,ma46,sma46,xrange=xr,yrange=yr
  map2fits,sma46,'maps/MOS'+mnm+'_EngCube_GO_FPMA_E46.fits'


  restore,file=dirout+'MOS'+mnm+'_EngCube_GO_FPMB.dat'
  imb2=total(dc_out.specs[where(dc_out.eng_edgs ge 2),*,*],1)/dc_out.lvtfac
  mb2=make_map(imb2,dx=dc_out.dx,dy=dc_out.dx,xc=dc_out.xc,yc=dc_out.yc,$
    time=dc_out.t1,dur=anytim(dc_out.t2)-anytim(dc_out.t1),id='FPMB G0 >2 keV',l0=l0,b0=ang[1],rsun=ang[2])
  sub_map,mb2,smb2,xrange=xr,yrange=yr
  map2fits,smb2,'maps/MOS'+mnm+'_EngCube_GO_FPMB_E2.fits'

  imb24=total(dc_out.specs[where(dc_out.eng_edgs ge 2 and dc_out.eng_edgs le 4),*,*],1)/dc_out.lvtfac
  mb24=make_map(imb24,dx=dc_out.dx,dy=dc_out.dx,xc=dc_out.xc,yc=dc_out.yc,$
    time=dc_out.t1,dur=anytim(dc_out.t2)-anytim(dc_out.t1),id='FPMB G0 2-4 keV',l0=l0,b0=ang[1],rsun=ang[2])
  sub_map,mb24,smb24,xrange=xr,yrange=yr
  map2fits,smb24,'maps/MOS'+mnm+'_EngCube_GO_FPMB_E24.fits'

  imb46=total(dc_out.specs[where(dc_out.eng_edgs ge 4 and dc_out.eng_edgs le 6),*,*],1)/dc_out.lvtfac
  mb46=make_map(imb46,dx=dc_out.dx,dy=dc_out.dx,xc=dc_out.xc,yc=dc_out.yc,$
    time=dc_out.t1,dur=anytim(dc_out.t2)-anytim(dc_out.t1),id='FPMB G0 4-6 keV',l0=l0,b0=ang[1],rsun=ang[2])
  sub_map,mb46,smb46,xrange=xr,yrange=yr
  map2fits,smb46,'maps/MOS'+mnm+'_EngCube_GO_FPMB_E46.fits'

;  stop
end
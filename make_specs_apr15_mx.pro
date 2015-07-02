pro make_specs_apr15_mx,mosid=mosid,chuid=chuid,fpmid=fpmid,dirout=dirout

  ; For saving it out
  ; shifts, smoothing, combinations can be done in another file
  if (n_elements(dirout) ne 1) then dirout='~/data/ns_data/obs4_maps/'
  if (n_elements(mosid) ne 1) then mosid=40;20
  mnm=string(mosid,format='(i2)')
  if (n_elements(chuid) ne 1) then chuid=4

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subdir='201100'+mnm+'_Sol_15119_MOS0'+mnm
  ddname='201100'+mnm+'001'

  chm=[5,9,14,13]
  chmn='CHU'+['012','003','123','023']
  chumask=chm[chuid-1]
  chunam=chmn[chuid-1]

  @post_outset

  if (n_elements(fpmid) ne 1) then fpmid='FPMA'
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ; The new regions of interest
  ; in S/C arcsec
  xcsp=[360]
  ycsp=[390]

  rnm=['R1']
  nr=n_elements(rnm)
  wid=120.

  restore,file=dirout+'MOS'+mnm+'_EngCube_GO_'+fpmid+'.dat'

  lvt=dc_out.lvtfac
  xc=dc_out.xc
  yc=dc_out.yc
  pxs=dc_out.dx
  x0=dc_out.x0
  y0=dc_out.y0
  t1=dc_out.t1
  t2=dc_out.t2
  dur=anytim(t2)-anytim(t1)
  engs=dc_out.eng_edgs
  mide=get_edges(engs,/mean)


  ; Get the region coords in terms of indexing of dc_out.specs
  xcs=round((xcsp-x0)/pxs)
  ycs=round((ycsp-y0)/pxs)
  iwid=round(wid/pxs)

  specs=fltarr(n_elements(mide))

  xr=round(xcs[0]+0.5*[-iwid,iwid])
  yr=round(ycs[0]+0.5*[-iwid,iwid])
  specs=total(total(dc_out.specs[*,xr[0]:xr[1],yr[0]:yr[1]],2),2)
  specs[0]=0.0


  str_out={engs:engs,mide:mide,specs:specs,$
    t1:t1,t2:t2,dur:dur,lvt:dc_out.lvtfac,$
    pxs:pxs,xcsp:xcsp,ycsp:ycsp,wid:wid,ids:rnm}

  savegen,file='specs/MOS'+mnm+'_'+fpmid,str_out

  ; need IDL >8.0 for ct 74
  !p.multi=[0,2,1]
  loadct,74,/silent
  reverse_ct
  tvlct,r,g,b,/get
  r[0]=0
  g[0]=0
  b[0]=0
  r[1]=255
  g[1]=255
  b[1]=255
  tvlct,r,g,b

  ims=gauss_smooth(total(dc_out.specs,1),2)

  set_plot,'ps'
  device, /encapsulated, /color, /isolatin1,/inches, $
    bits=8, xsize=7, ysize=4,file='figs/Specsall_MOS'+mnm+'_'+fpmid+'.eps'
  !p.charsize=1.0

  plot_image,ims^0.5,tit='MOS'+mnm+' '+fpmid,xtickf='(a1)',ytickf='(a1)',bot=1
  box_igh,xcs[0],ycs[0],iwid,color=0,thick=2
  xyouts,xcs[0]+0.5*iwid+5,ycs[0]-0.5*iwid-5,rnm[0],color=0,chars=0.7

  plot,mide,specs,title=rnm[0],xtit='Energy [keV]',ytit='NuSTAR counts',$
    /ylog,yrange=[0.9,700],xrange=[1,7],psym=10,ytickf='exp1',thick=4

  device,/close
  set_plot, mydevice

  ;  stop
end

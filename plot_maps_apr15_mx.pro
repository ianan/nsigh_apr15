pro plot_maps_apr15_mx,eid=eid,mosid=mosid,chuid=chuid

  ; For saving it out
  ; shifts, smoothing, combinations can be done in another file
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

  ;
  ; Warning - need a subdirectory called "maps" in your current working directory where the maps are
  ;           need a subdirectory called "figs" in your current working directory to put the figs in
  ;           For colortable 74 need IDL >8.0

  ;
  ; 3-Jun-2015  IGH

  ; What energy range do you want?
  if (n_elements(eid) ne 1) then eid='E2' ;'E24','E46'

  if (eid eq 'E2') then ename='>2keV'
  if (eid eq 'E24') then ename='2-4keV'
  if (eid eq 'E46') then ename='4-6keV'

  fits2map,'maps/MOS'+mnm+'_EngCube_GO_FPMA_'+eid+'.fits',ma
  fits2map,'maps/MOS'+mnm+'_EngCube_GO_FPMB_'+eid+'.fits',mb

  mab=ma
  mab.data=ma.data+mb.data
  mab.id='FPMA+FPMB G0 '+ename

  @post_outset

  ; Need IDL >8.0 for this color table
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

  dnl=1.5
  dmx=3

  xrn=[-700,700]
  yrn=[-300,1100]

  set_plot,'ps'
  device, /encapsulated, /color, / isolatin1,/inches, $
    bits=8, xsize=5, ysize=5,file='figs/MOS'+mnm+'_EngCube_GO_FPMA_'+eid+'.eps'
  plot_map,ma,/log,dmin=10^dnl,dmax=10^dmx,title='FPMA G0 '+ename+' '+anytim(ma.time,/yoh,/trunc),grid_sp=30,$
    gcolor=0,bot=1,position=[0.15,0.1,0.95,0.9],chars=1,xrange=xrn,yrange=yrn
  plot_map_cb_igh,[dnl,dmx],position=[0.2,0.25,0.5,0.27],color=0,chars=1,$
    cb_title='NuSTAR [log!D10!N counts]',bottom=1,format='(f4.1)'

  device,/close
  set_plot, mydevice

  set_plot,'ps'
  device, /encapsulated, /color, / isolatin1,/inches, $
    bits=8, xsize=5, ysize=5,file='figs/MOS'+mnm+'_EngCube_GO_FPMB_'+eid+'.eps'
  plot_map,mb,/log,dmin=10^dnl,dmax=10^dmx,title='FPMB G0 '+ename+' '+anytim(ma.time,/yoh,/trunc),grid_sp=30,$
    gcolor=0,bot=1,position=[0.15,0.1,0.95,0.9],chars=1,xrange=xrn,yrange=yrn
  plot_map_cb_igh,[dnl,dmx],position=[0.2,0.25,0.5,0.27],color=0,chars=1,$
    cb_title='NuSTAR [log!D10!N counts]',bottom=1,format='(f4.1)'

  device,/close
  set_plot, mydevice
  
 
  ;
  ;  set_plot,'ps'
  ;  device, /encapsulated, /color, / isolatin1,/inches, $
  ;    bits=8, xsize=5, ysize=5,file='figs/MOS'+mnm+'_EngCube_GO_FPMAB_'+eid+'.eps'
  ;  plot_map,mab,/log,dmin=10^dnl,dmax=10^dmx,title='FPMA+B G0 '+ename+' '+anytim(ma.time,/yoh,/trunc),grid_sp=30,$
  ;    gcolor=0,bot=1,position=[0.15,0.1,0.95,0.9],chars=1,xrange=xrn,yrange=yrn
  ;  plot_map_cb_igh,[dnl,dmx],position=[0.2,0.25,0.5,0.27],color=0,chars=1,$
  ;    cb_title='NuSTAR [log!D10!N]',bottom=1,format='(f4.1)'
  ;
  ;  device,/close
  ;  set_plot, mydevice
  ;************************************************************

  sr=2
  mas=ma
  mas.data=gauss_smooth(ma.data,sr)
  mbs=mb
  mbs.data=gauss_smooth(mb.data,sr)
  mabs=mab
  mabs.data=gauss_smooth(mab.data,sr)


  set_plot,'ps'
  device, /encapsulated, /color, / isolatin1,/inches, $
    bits=8, xsize=5, ysize=5,file='figs/MOS'+mnm+'_EngCube_GO_FPMA_'+eid+'_S.eps'
  plot_map,mas,/log,dmin=10^dnl,dmax=10^dmx,title='FPMA G0 '+ename+' '+anytim(ma.time,/yoh,/trunc),grid_sp=30,$
    gcolor=0,bot=1,position=[0.15,0.1,0.95,0.9],chars=1,xrange=xrn,yrange=yrn
  plot_map_cb_igh,[dnl,dmx],position=[0.2,0.25,0.5,0.27],color=0,chars=1,$
    cb_title='NuSTAR [log!D10!N counts]',bottom=1,format='(f4.1)'

  device,/close
  set_plot, mydevice

  set_plot,'ps'
  device, /encapsulated, /color, / isolatin1,/inches, $
    bits=8, xsize=5, ysize=5,file='figs/MOS'+mnm+'_EngCube_GO_FPMB_'+eid+'_S.eps'
  plot_map,mbs,/log,dmin=10^dnl,dmax=10^dmx,title='FPMB G0 '+ename+' '+anytim(ma.time,/yoh,/trunc),grid_sp=30,$
    gcolor=0,bot=1,position=[0.15,0.1,0.95,0.9],chars=1,xrange=xrn,yrange=yrn
  plot_map_cb_igh,[dnl,dmx],position=[0.2,0.25,0.5,0.27],color=0,chars=1,$
    cb_title='NuSTAR [log!D10!N counts]',bottom=1,format='(f4.1)'

  device,/close
  set_plot, mydevice
  ;
  ;
  ;  set_plot,'ps'
  ;  device, /encapsulated, /color, / isolatin1,/inches, $
  ;    bits=8, xsize=5, ysize=5,file='figs/MOS'+mnm+'_EngCube_GO_FPMAB_'+eid+'_S.eps'
  ;  plot_map,mabs,/log,dmin=10^dnl,dmax=10^dmx,title='FPMA+B G0 '+ename+' '+anytim(ma.time,/yoh,/trunc),grid_sp=30,$
  ;    gcolor=0,bot=1,position=[0.15,0.1,0.95,0.9],chars=1,xrange=xrn,yrange=yrn
  ;  plot_map_cb_igh,[dnl,dmx],position=[0.2,0.25,0.5,0.27],color=0,chars=1,$
  ;    cb_title='NuSTAR [log!D10!N]',bottom=1,format='(f4.1)'
  ;
  ;  device,/close
  ;  set_plot, mydevice


;    stop
end
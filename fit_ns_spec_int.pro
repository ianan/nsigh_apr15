pro fit_ns_spec_int,mosid=mosid,fpmid=fpmid,chuid=chuid,regid=regid, $
  fiter=fiter,pint=pint,puf=puf,newe=newe,monte=monte,plter=plter,ab=ab

  ; Do a single thermal spectral fit to the April data
  ; By default will whole AR region (R1) for MOS40 and FPMA over 2.5-7 keV
  ; with no pileup correction and using the original 0.04 keV binning

  ; Optional inputs
  ; mosid   -  ID of mosiac tile/file (10-26 for Orbit 1, 30-46 for Orbit 2)
  ; fpmid   -  FPM to use ('FPMA' or 'FPMB', default is 'FPMA')
  ; chuid   -  Which CHU state to filter for (default here is CHU23 in MOS40)
  ; regid   -  Region to fit (just 0 here as only one region)
  ; fiter  -  Energy range to fit over (default is 2.5 to 7 keV]
  ; pint    -  Ininitial parameters for the fit (default is [42,0.5,1.0])
  ; puf     -  Pileup fraction to do the pileup correction (default is 0.)
  ;                   should be 0.25*(grade 21-24)/(all grades)
  ; newe    -  New energy binning of the array (deafult is 0, don't do rebinning)
  ; monte   -  Do the monte carlo analysis on the errors?
  ; plter   -  Energy range to plot
  ; ab      -  Keyword to do a combined spectral fit for A+B (default 0) - fpmid has no effect for /ab

  ; Warning -  the spectra to fit are expected to be in subdir "specs" (also output location for fit info)
  ;            the arf/rmf from nupipeline or the resulting rsp file should be in subdir "rsp"
  ;            a count spectrum (data, fit and residuals) will go into subdir "figs"

  ; 15-Jun-15 IGH - Begun as modification of Nov2014 analysis script
  ; 24-Jun-15 IGH - Correctly inputting spectrum, errors and livetime (latter is livetime frac*duration)
  ;                    DRM should be DRM=RMF*ARF/dE to get it in the units ospex wants
  ;                    Monte-carlo error analysis also included and works now for NuSTAR data
  ; 25-Jun-15 IGH - Include link to *.pro to make the rsp from rmf and arf (which you'd also need)
  ; 29-Jun-15 IGH - Added plter option to change the energy range of the output spectra
  ;               - Added keyword to do combined A+B spectral fit

  @post_outset
  ; Which pointing to use
  if (n_elements(mosid) ne 1) then mosid=40;20
  mnm=string(mosid,format='(i2)')
  if (n_elements(chuid) ne 1) then chuid=4
  ; Which FPM to use
  if (n_elements(fpmid) ne 1) then fpmid='FPMA'
  ; Which region to do
  if (n_elements(regid) ne 1) then regid=0
  ; Energy range to do the fit over
  if (n_elements(fiter) ne 2) then fiter=[2.5, 7.]
  ; Starting point of the fit parameters
  if (n_elements(pint) ne 3) then pin=[42, 0.5, 1]
  ; puf is 0.25*(grade21-24/all grades)
  if (n_elements(puf) eq 0) then puf=0.0;0.11*.25/100.
  ; New energy binning - this is the energy bin midpoints
  if (n_elements(newe) eq 0) then newe=0.;1.7+findgen(42)*0.2;1.8+findgen(21)*0.4
  ; What energy range to plot over ?
  if (n_elements(plter) ne 2) then plter=[1.6,7.2]

  if keyword_set(ab) then begin
    fpmid='FPMAB'
    restgen,file='specs/MOS'+mnm+'_FPMA',ss
    restgen,file='specs/MOS'+mnm+'_FPMB',ssb
    ; normalise the B spectrum to the A livetime (in secs)
    spectrum=reform(ss.specs)+(ss.lvt*ss.dur)*reform(ssb.specs)/(ssb.lvt*ssb.dur)
    nspec=n_elements(spectrum)
    sengs=ss.engs
    smide=ss.mide
    regnm=ss.ids[regid]
    dur=ss.dur
    lvt=ss.lvt
   endif else begin
    restgen,file='specs/MOS'+mnm+'_'+fpmid,ss
    spectrum=reform(ss.specs)
    nspec=n_elements(spectrum)
    sengs=ss.engs
    smide=ss.mide
    regnm=ss.ids[regid]
    dur=ss.dur
    lvt=ss.lvt
  endelse

  ; Do the interpolation to a different energy binning?
  ; Interpol the best way of doing this - lose stuff with the odd few counts ???
  intnam=''
  if (n_elements(newe) gt 1) then begin
    intnam='INT'+string(1000+100*(newe[1]-newe[0]),format='(i4)')
    spectrum=interpol(spectrum/(sengs[1]-sengs[0]),smide,newe)*(newe[1]-newe[0])
    smide=newe
    nspec=n_elements(smide)
    sengs=[newe-0.5*(newe[1]-newe[0]),newe[nspec-1]+0.5*(newe[1]-newe[0])]
  endif

  ; Do the rough pileup correction
  ;   subtract observation*pileup fraction shifted to double energy
  puftit=''
  pufnam=''
  if (puf gt 0.) then begin
    puftit=string(puf*1e2,format='(f6.2)')+'%'
    pufnam='_'+strmid(string(1e6+puf*1e6,format='(i7)'),1,6)
    pdu_spec=interpol(puf*spectrum,2*smide,smide)
    idg=where(smide ge 2*min(smide))
    spectrum[idg]=(spectrum[idg]-pdu_spec[idg]) >0.
  endif

  if keyword_set(ab) then begin
    fl_rsp=['rsp/'+regnm+'FPMA'+intnam+'_RSP.genx','rsp/'+regnm+'FPMB'+intnam+'_RSP.genx']
    fpmids=['FPMA','FPMB']
  endif else begin
    fl_rsp='rsp/'+regnm+fpmid+intnam+'_RSP.genx'
    fpmids=fpmid
  endelse

  ; Only make the DRM file if it doesn't already exist
  for pp=0,n_elements(fpmids)-1 do begin
    if (file_test(fl_rsp[pp]) eq 0) then begin
      nupdir='rsp/'+regnm+fpmids[pp]+'/'
      ; These are the outputs from nupipeline
      flrmf=file_search(nupdir,'*.rmf')
      flarf=file_search(nupdir,'*.arf')

      ; The idl procedures to make these are mostly at
      ; https://lost-contact.mit.edu/afs/physics.wisc.edu/home/craigm/lib/idl/spectral/
      ; vcol2arr.pro and pointer_value.pro are at
      ; https://lost-contact.mit.edu/afs/physics.wisc.edu/home/craigm/lib/idl/util/
      rmfread, flrmf, rmf_str, compressed='UNCOMPRESSED'
      rmf = transpose(*(rmf_str.data))
      edges = *(rmf_str.ebins)
      fxbopen, unit, flarf, 'SPECRESP', hh, errmsg=err
      fxbreadm, unit, ['ENERG_LO','ENERG_HI','SPECRESP'], elo, ehi, arf, errmsg=err
      fxbclose, unit
      rsp=rmf
      for rr = 0, n_elements(rmf[0,*])-1 do rsp[*,rr]=rmf[*,rr]*arf[rr]

      edges=edges[*,0:n_elements(ss.mide)-1]
      rsp=rsp[0:n_elements(ss.mide)-1,0:n_elements(ss.mide)-1]

      ; Need to change the energy binning for rsp?
      if (n_elements(newe) gt 1) then begin
        new_edges=fltarr(2,nspec)
        new_edges[0,*]=sengs[0:nspec-1]
        new_edges[1,*]=sengs[1:nspec]
        new_rsp_t=fltarr(nspec,n_elements(ss.mide))
        for yy=0,n_elements(ss.mide)-1 do new_rsp_t[*,yy]=interpol(rsp[*,yy],ss.mide,smide)
        new_rsp=fltarr(nspec,nspec)
        for xx=0, nspec-1 do new_rsp[xx,*]=interpol(new_rsp_t[xx,*],ss.mide,smide)
        edges=new_edges
        rsp=new_rsp
      endif
      ; Seems that we need to divide by dE in here - otherwise EM varies with dE !!
      rsp=rsp
      savegen,file=fl_rsp[pp],edges,rsp
    endif
  endfor
  ; load the DRM
  ; If doing A+B combined need to load both and add
  ; assumes already have run a and b separately
  if keyword_set(ab) then begin
    fl_rspa='rsp/'+regnm+'FPMA'+intnam+'_RSP.genx'
    fl_rspb='rsp/'+regnm+'FPMB'+intnam+'_RSP.genx'
    restgen,file=fl_rspa,edges,rsp
    restgen,file=fl_rspb,edgesb,rspb
    rsp=rsp+rspb
  endif else begin
    restgen,file=fl_rsp,edges,rsp
  endelse

  dE=edges[1,0]-edges[0,0]

  ; Assuming saved rsp is just RMF*ARF
  ; so need to /dE to get into (counts/photons)/keV can check via o->getdata(class=’spex_drm’)
  rsp=rsp/dE

  ; start up ospec in noninteractive mode
  set_logenv, 'OSPEX_NOINTERACTIVE', '1
  o = ospex()
  o->set,  spex_autoplot_enable=0, spex_fit_progbar=0, spex_fitcomp_plot_resid=0,spex_fit_manual=0
  o->set, spex_data_source='SPEX_USER_DATA'
  z = o->get(/spex_data_origunits)
  z.data_name = 'NuSTAR'
  o->set, spex_data_origunits=z
  o->set, spex_detectors=fpmid
  o->set, spex_ut_edges=anytim([ss.t1,ss.t2])
  o->set, spex_fit_time_interval=o->get(/spex_ut_edges)
  o->set, spex_ct_edges=edges

  ; ss.spec just counts so need to convert to rate
  spec_in=spectrum
  ; error estimate just Poisson on actual detected counts then convert to lvtcor rate
  err_in=sqrt(spectrum)
  ; Set the livetime
  ltime_in=fltarr(n_elements(spectrum))+(dur*lvt)
  o->set, spectrum=spec_in, errors=err_in, livetime=ltime_in
  ; for speed just include the response that covers the energy range of the spectrum
  o->set, spex_respinfo=rsp

  ; setup the fit and then do it
  o->set, spex_erange=fiter
  o->set, fit_function='vth'
  o->set, fit_comp_minima= [1e-10, 0.1, 0.5]
  o->set, fit_comp_maxima= [1e10,2.0, 1.5]
  o->set, fit_comp_param=pin
  o->set, fit_comp_free = [1, 1, 0]
  o->dofit, /all

  ; get the resulting fitted spectrum model and paramteres
  model = o -> calc_func_components(spex_units='flux',/counts,/all_func,this_interval=0)
  ;  obs=o->calc_summ(item='data_count_flux',errors=err,this_interval=0)
  ; Doing it this way below to check the model actually matches the data we input to ospex
  obs=spec_in/dE/lvt/dur
  err=err_in/dE/lvt/dur

  pmodel = o -> calc_func_components(spex_units='flux',/photons,/all_func,this_interval=0)
  pobs=o->calc_summ(item='data_photon_flux',errors=perr,this_interval=0)
  id0=where(obs eq 0,nid0)
  if (nid0 gt 0) then begin
    pobs[id0]=0.
    perr[id0]=0.
  endif

  engs=get_edges(model.ct_energy,/mean)
  np=n_elements(engs)

  parm=o->get(/spex_summ_params)
  parmerr=o->get(/spex_summ_sigmas)
  chisq=o->get(/spex_summ_chisq)
  tkev=0.08617

  eranfit=o->get(/spex_erange)

  ; Save all the spectra and fit info out
  res={eranfit:eranfit,parm:parm,parmerr:parmerr,tmk:parm[1]/tkev,$
    engs:engs,cnt_mod:model.yvals[*],cnt_flx:obs,ecnt_flx:err,$
    ph_mod:pmodel.yvals[*],pht_flx:pobs,epht_flx:perr,chisq:chisq}

  savegen,file='specs/MOS'+mnm+'_'+regnm+'_'+fpmid+pufnam+intnam+'_fit_vth',res

  ; Make a plot of the spectra and fit
  loadct,0,/silent
  tube_line_colors
  !p.thick=2
  @post_outset
  !p.multi=[0,1,2]
  set_plot,'ps'
  device, /encapsulated, /color, /isolatin1, $
    /inches, bits=8, xsize=4, ysize=5,file='figs/SCN_MOS'+mnm+'_'+regnm+'_'+fpmid+pufnam+intnam+'_VTH.eps'

  ylim=[2e0,2e4]
  xlim=plter
  gd=where(engs lt 15,ngd)

  plot,[1,1],[1,1],xrange=xlim,/ylog,$
    yrange=ylim,ytit='counts s!U-1!N cm!U-2!N keV!U-1!N',ytickf='exp1',$
    title='MOS'+mnm+' '+fpmid+' '+regnm+' '+puftit,position=[0.175,0.3,0.975,0.94],xtickf='(a1)';,xtit='Energy [keV]'
  oploterr,engs,obs,err,color=0,psym=10,thick=5,/nohat,errthick=3,errcolor=150
  oplot,engs,model.yvals[*],color=1,thick=5

  oplot,eranfit[0]*[1,1],ylim,line=2,thick=2,color=150
  oplot,eranfit[1]*[1,1],ylim,line=2,thick=2,color=150

  xyouts,9.2e3,10.8e3,string(parm[1]/tkev,format='(f4.1)')+string(177B)+string(parmerr[1]/tkev,format='(f4.2)')+' MK ('+$
    string(parm[1],format='(f4.2)')+' keV)',color=1,chars=1.,/device,align=1
  xyouts,9.2e3,10.3e3,string(parm[0]*1e2,format='(f4.1)')+string(177B)+string(parmerr[0]*1e2,format='(f4.2)')+$
    ' x10!U47!N cm!U-3!N',color=1,chars=1.,/device,align=1
  xyouts,9.2e3,9.8e3, '!Mc!3!U2!N: '+strcompress(string(chisq,format='(f5.1)'),/rem),color=1,chars=1.,/device,align=1

  plot,xlim,[-4.5,4.5],xr=xlim,xtit='Energy [keV]',/nodata,position=[0.175,0.1,0.975,0.3],ytit='(Obs-Mod)/Err'
  oplot,xlim,[0,0],line=1,color=150
  oplot,engs,(obs-model.yvals[*])/err,thick=4,psym=10

  oplot,eranfit[0]*[1,1],[-4.5,4.5],line=2,thick=2,color=150
  oplot,eranfit[1]*[1,1],[-4.5,4.5],line=2,thick=2,color=150

  device,/close
  set_plot, mydevice


  ;  stop

  ;  ; Make a plot of the spectra in photon space
  ;  set_plot,'ps'
  ;  device, /encapsulated, /color, /HELVETICA, $
  ;    /inches, bits=8, xsize=4, ysize=4,file='figs/SPH_MOS'+mnm+'_'+regnm+'_'+fpmid+'_VTH.eps'
  ;
  ;  ylim=[1e-1,1e8]
  ;  xlim=[1.6,7.2]
  ;  gd=where(engs lt 15,ngd)
  ;
  ;  plot,[1,1],[1,1],xtit='Energy [keV]',xrange=xlim,/ylog,$
  ;    yrange=ylim,ytit='photons s!U-1!N cm!U-2!N keV!U-1!N',ytickf='exp1',$
  ;    title='MOS'+mnm+' '+fpmid+' '+regnm
  ;  oploterr,engs,pobs,perr,color=0,psym=10,thick=5,/nohat,errthick=3,errcolor=150
  ;  oplot,engs,pmodel.yvals[*],color=1,thick=5
  ;
  ;  oplot,eranfit[0]*[1,1],ylim,line=2,thick=2,color=150
  ;  oplot,eranfit[1]*[1,1],ylim,line=2,thick=2,color=150
  ;
  ;  xyouts,9.0e3,8.5e3,string(parm[1]/tkev,format='(f4.1)')+' MK ('+$
  ;    string(parm[1],format='(f4.2)')+' keV)',color=1,chars=1.,/device,align=1
  ;  xyouts,9.0e3,8.0e3,string(parm[0]*1e2,format='(f4.1)')+' x10!U47!N cm!U-3!N',color=1,chars=1.,/device,align=1
  ;  xyouts,9.0e3,7.5e3, '!Mc!3!U2!N: '+strcompress(string(chisq,format='(f5.1)'),/rem),color=1,chars=1.,/device,align=1
  ;
  ;  device,/close
  ;  set_plot, mydevice

  ;  stop

  ; Do the monte carlo error analysis on the best fit parameters
  if keyword_set(monte) then begin
    ; by default assume want to do 1e4 iterations
    o -> monte_carlo, niter=1e4, interval=0,$
      savefile='specs/MC_MOS'+mnm+'_'+regnm+'_'+fpmid+pufnam+intnam+'_fit_vth.sav'
    restore,file='specs/MC_MOS'+mnm+'_'+regnm+'_'+fpmid+pufnam+intnam+'_fit_vth.sav'
    spex_monte_carlo_results, $
      savefile='specs/MC_MOS'+mnm+'_'+regnm+'_'+fpmid+pufnam+intnam+'_fit_vth.sav', out_struct=mcs

    nx=74.
    ny=74.
    p1 = reform(savep[0,*])
    p2 = reform(savep[1,*])
    dx = (max(p1) - min(p1)) / (nx-1)
    dy = (max(p2) - min(p2)) / (ny-1)
    x = min(p1) + findgen(nx)*dx
    y = min(p2) + findgen(ny)*dy
    xbin = value_locate(x, p1)
    ybin = value_locate(y, p2)
    map = fltarr(nx,ny)
    chi=savechi
    for m = 0L,n_elements(chi)-1 do map[xbin[m], ybin[m]]++

    xrmc=minmax(x)*1d2
    yrmc=minmax(y)/tkev

    loadct,39
    !p.multi=0
    set_plot,'ps'
    device, /encapsulated, /color, /isolatin1, $
      /inches, bits=8, xsize=4, ysize=4,file='figs/SCN_MOS'+mnm+'_'+regnm+'_'+fpmid+pufnam+intnam+'_VTH_MC.eps'
    tvim,map,/norm,xrange=xrmc,yrange=yrmc,pcharsize=1.,$
      xtit='EM x10!U47!N cm!U-3!N',ytit='T [MK]',title='Fit err (white), MC 67% (red)'

    em=parm[0]*1d2
    eme=parmerr[0]*1d2
    tmk=parm[1]/tkev
    tmke=parmerr[1]/tkev

    oplot,em*[1,1],yrmc,color=255
    oplot,xrmc,tmk*[1,1],color=255

    oplot,(em+eme)*[1,1],yrmc,color=255
    oplot,(em-eme)*[1,1],yrmc,color=255
    oplot,xrmc,(tmk+tmke)*[1,1],color=255
    oplot,xrmc,(tmk-tmke)*[1,1],color=255

    oplot,1e2*mcs[0].mode*[1,1],yrmc,lines=2,color=250
    oplot,xrmc,mcs[1].mode*[1,1]/tkev,lines=2,color=250
    oplot,1e2*mcs[0].v67[0]*[1,1],yrmc,lines=2,color=250
    oplot,1e2*mcs[0].v67[1]*[1,1],yrmc,lines=2,color=250
    oplot,xrmc,mcs[1].v67[0]*[1,1]/tkev,lines=2,color=250
    oplot,xrmc,mcs[1].v67[1]*[1,1]/tkev,lines=2,color=250

    device,/close
    set_plot, mydevice

  endif

end

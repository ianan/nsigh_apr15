pro batch_do_all

  ; A quick analysis of the 29-April-2015 data to get and fit the spectra of the AR which XRT also observed
  ; Just looking at one mosaic tile
  ;
  ; Produce the data cube (e,x,y)
  make_datacube_apr15_mx
  
  ; Produce a map per FPM and energy (>2, 2-4, 4-6)
  make_maps_apr15_mx
  
  ; Plot the maps
  plot_maps_apr15_mx,eid='E2'
  
  make_specs_apr15_mx,fpmid='FPMA'
  make_specs_apr15_mx,fpmid='FPMB'
  
  ; Got to ds9 to make the *.reg file - so it original *cl.evt not the sunpos one
  ; Then run nupipeline to get the *.arf and *.rmf - using make_arfrmf.txt and have nuproducts running
  
  fit_ns_spec_int,fpmid='FPMA',plter=[2.3,5.2],fiter=[2.5,5.0];,/monte
  fit_ns_spec_int,fpmid='FPMB',plter=[2.3,5.2],fiter=[2.5,5.0];,/monte
  ; Do them for A+B combined
  fit_ns_spec_int,/ab,plter=[2.3,5.2],fiter=[2.5,5.0]
  
  ; Do the fits again but using 0.2 keV binning
  fit_ns_spec_int,fpmid='FPMA',plter=[2.3,5.2],fiter=[2.5,5.0],newe=1.7+findgen(42)*0.2;,/monte
  fit_ns_spec_int,fpmid='FPMB',plter=[2.3,5.2],fiter=[2.5,5.0],newe=1.7+findgen(42)*0.2;,/monte
  ; Do them for A+B combined
  fit_ns_spec_int,/ab,plter=[2.3,5.2],fiter=[2.5,5.0],newe=1.7+findgen(42)*0.2

  stop
end
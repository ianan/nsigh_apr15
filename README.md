# nsigh_apr15
Scripts for a quick end to end analysis of a mosaic tile from the 29-April-2015 observation - this can take you from the processed NuSTAR events files to a fitted OSPEX spectrum.

For more help see https://heasarc.gsfc.nasa.gov/docs/nustar/analysis/nustar_swguide.pdf

Other than these scripts you will need the following (though some parts will work without all the bits)
  - The NuSTAR data from the April 2015 observations, in the original directory structure, post pipeline to stage 2, and co-ordinates converted to sun centered. Need this to create the data cubes (using make_data_cube_apr15_mx.pro).
  - Ds9 (http://ds9.si.edu/site/Home.html) to produce the *.reg region file - one is alredy provided here ds9_R1.reg
  - NuSTAR software installed to run nuproducts (c.f. make_arfrmf.txt) on the data to generate *.arf and *.rmf files (which are used by fit_ns_spec_int.pro to make the response matrix, used in the fitting)
  - Some extra IDL scripts used to make the response matrix from *.rmf and *.arf - available here https://lost-contact.mit.edu/afs/physics.wisc.edu/home/craigm/lib/idl/spectral/ and https://lost-contact.mit.edu/afs/physics.wisc.edu/home/craigm/lib/idl/util/
  - IDL (at least version 8 - not tested on earlier ones)
 
Running batch_do_all.pro will take you through all the steps and each code is documented individually, basically need to do
  - make_datacube_apr15_mx            
    - Turns the eventlist into a datacube (E,x,y)
  - make_maps_apr15_mx                
    - Turns the datacube into maps of >2kev, 2-4 kV and 4-6keV
  - plot_maps_apr15_mx,eid='E2'  
    - Plot one of the maps
  - make_specs_apr15_mx,fpmid='FPMA'
    - Makes the spectrum for a region from the datacube 
  - Use nuproducts (make_arfmf.txt) to make the *.arf and *.rmf
    - none supplied here 
  - fit_ns_spec_int,fpmid='FPMA',plter=[2.3,5.2],fiter=[2.5,5.0]
    - Generates the /rsp/*RSP.genx and does the single thermal OSPEX fit 
  



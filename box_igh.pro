pro box_igh,xc,yc,wid,color=color,thick=thick
  ; Draw a box on another plot - centered at xc,yc of width wid
  ; IGH
  if (n_elements(color) ne 1) then color=0
  if (n_elements(thick) ne 1) then thick=1
  xr=xc+0.5*[-wid,wid]
  yr=yc+0.5*[-wid,wid]
  oplot,xr,yr[0]*[1,1],color=color,thick=thick
  oplot,xr,yr[1]*[1,1],color=color,thick=thick
  oplot,xr[0]*[1,1],yr,color=color,thick=thick
  oplot,xr[1]*[1,1],yr,color=color,thick=thick

end
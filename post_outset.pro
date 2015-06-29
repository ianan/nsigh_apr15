; just setup some of the plotting and postscript options
; IGH

clearplot
!y.range = 0
!x.range = 0
!p.title = ''
!x.title = ''
!y.title = ''
!y.style = 17
!x.style = 17
current_plot = !d.name
set_plot,'x'
device,retain = 2, decomposed=0
mydevice = !d.name
;; Make IDL use device/hardware fonts
!p.font = 0
;; Other things to make graphs nicer
!p.color = 255  ;255 for white line
!p.background = 0   ;0 for balck background
!p.thick = 1
!p.charthick = 1
!p.charsize = 1
!p.symsize = 1
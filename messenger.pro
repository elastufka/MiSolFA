;;;;;;;;;;;;;;;;;;;;;;;;;;
;messenger 2014 September 1
;;;;;;;;;;;;;;;;;;;;;;;;;;
this_file='xrs2014244_data.fits'
;counts per s per cm2 per keV
ccc=mrdfits(this_file,1)
;energy bins in keV
eee=mrdfits(this_file,2)
;get center of each energy bin
this_e=(eee.e_max+eee.e_min)/2.

;select two energy range (below calibration is less certain)
elist0=where( (this_e ge 2.6) AND (this_e le 5.8) )
elist1=where( (this_e ge 6.1) AND (this_e le 7.3) )

;sum over energy bins
m0=total(ccc.rate(elist0,*),1)
m1=total(ccc.rate(elist1,*),1)
m0_sep=m0

;time in utplot format
mtime=anytim('2014-09-01')+ccc.time


;popen,xsize=20./2.54,ysize=10/2.54
clear_utplot
loadct2,5
th3=4
chs=1.4
utplot,mtime,m0/max(m0),xstyle=1,timerange=[mtime(45),mtime(150)],yrange=[0,1.05],ystyle=1,title='messenger XRS: count rate (solid) and derivative (dotted)'
outplot,mtime,m1/max(m1),color=1,thick=th3
outplot,mtime,m0/max(m0),color=6,thick=th3
outplot,mtime,deriv(m1)/max(deriv(m1)),color=1,linestyle=2,thick=th3
outplot,mtime,deriv(m0)/max(deriv(m0)),color=6,linestyle=2,thick=th3
xyouts,0.14,0.85,'2.6-5.8 keV',color=6,/normal,charsize=chs
xyouts,0.14,0.78,'6.1-7.3 keV',color=1,/normal,charsize=chs
;pclose

;here the two routines to calculate the GOES class given flare temperature and EM:

goes_fluxes,30,1,flong,fshort,sat=15
print,goes_value2class(flong)

;example is for T=30 MK and EM=1d49 cm-3 

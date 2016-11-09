;energy in keV
eedd=findgen(600-5)/2.+3.


;thermal emission for MIXI proposal flare
this_em=7d47/1d49          ;in units of 10d49
this_te=20./11.6           ;in keV   1 keV = 11.6 M
;to get thermal spectrum in photons/s/cm2/keV
thth=f_vth(eedd,[this_em,this_te])


;nonthermal spectrum in photons/s/cm2/keV

this_norm=3         ;norm (value of lower power law at 50 keV)
this_below=1.7      ;power law index below break
this_break=13.5     ;break energy in keV
this_gamma=3.7      ;power law index above break


ntnt=f_bpow(eedd,[this_norm,this_below,this_break,this_gamma])

loadct,5
plot_oo,eedd,thth,xtitle='energy [keV',ytitle='photons s!U-1!N cm!U-2!N keV!U-1!N',xstyle=1,yrange=[0.001,max(thth)]
oplot,eedd,ntnt,color=44
oplot,eedd,thth,color=122

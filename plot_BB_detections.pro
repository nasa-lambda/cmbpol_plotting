pro plot_BB_detections

; This makes a plot of observed B-mode power from CMB experiments with
; significant detections. Data from BICEP2+Keck October 2015, 
; BICEP2+Keck/Planck January 2015, POLARBEAR, and SPTpol are included.
; The plotted BICEP2+Keck data are the CMB component from a spectral 
; decomposition into CMB, dust, and synchrotron components. The plotted 
; BICEP2+Keck/Planck data have had the dust foreground subtracted based 
; on the measured cross-power between Planck and BICEP2+Keck. The other
; data plotted have not had any dust foreground subtraction.
; For comparison, a theoretical curve for a LCDM model with tensor-to-scalar
; ratio r=0.1 is also plotted as a solid curve. The inflationary and 
; gravitational lensing components are plotted separately as dashed
; and dotted curves, respectively.
; The data are plotted to a postscript file BB_detections.ps.

; The experimental data are read from a file bb_data_2015nov.txt, which
; should be copied to the user's local directory.
; Sources for the data and more information are given in 
; http://lambda.gsfc.nasa.gov/graphics/bb_upperlimits/

; The theoretical data are read from files made by the BICEP2 team,
; from their calculations using the March 2013 version of CAMB.
; They should be copied to the user's local directory from
; http://lambda.gsfc.nasa.gov/data/suborbital/BICEP2/B2_3yr_camb_planck_lensed_uK_20140314.txt
; http://lambda.gsfc.nasa.gov/data/suborbital/BICEP2/B2_3yr_camb_planck_withB_uK_20140314.txt

; This code calls readcol.pro and associated routines from the IDL 
; Astronomy User's Library, http://idlastro.gsfc.nasa.gov/


; read data from experiments with BB detections

readcol,'bb_data_2015nov.txt',format='A,I,F,I,F,F,F',experiment,$
 l_min,l_center,l_max,BB,sigma_BB_minus,sigma_bb_plus,skipline=3,numline=36

set_plot,'ps'

device,filename='BB_detections.ps',/color,bits=8

!p.thick=3
!p.charthick=3
!x.thick=3
!y.thick=3

loadct,12

; set up plot without data

plot_oo,[1.,1e+4],[1.e-2,1.e+3],ps=3,xtitle='Multipole !8l!3',$
 ytitle='!8l!3(!8l!3+1)C!d!8l!3!u!8BB!n/!32!4p  !3[!4l!3K!u2!n]',$
 yr=[1.e-3,0.6e+0],ys=1,xr=[1.,4000],xs=1,/nodata


; add BICEP2+Keck/Planck detections to plot

s=where(experiment eq 'BICEP2+Keck/Planck')

oplot,l_center(s(1:*)),bb(s(1:*)),ps=4,symsize=0.6,color=70

errplot,l_center(s(1:*)),bb(s(1:*))+sigma_bb_plus(s(1:*)),$
 bb(s(1:*))-sigma_bb_minus(s(1:*)),width=0,color=70

l_mid=(l_min+l_max)/2.

halfwidth=l_mid-l_min

xerr,l_mid(s(1:*)),bb(s(1:*)),halfwidth(s(1:*)),color=70

; oplot 95% confidence upper limit for first BICEP2+Keck/Planck bin
;  where foreground-cleaned BB < 0
 
; get upper limit from Gaussian likelihood truncated at zero,
;  take mean of 10 realizations

realiz=fltarr(10)

for i=0,9 do begin

  get_limit,bb(s(0)),sigma_bb_plus(s(0)),0.95,limit

  realiz(i) = limit

endfor

usersym,[-1,1,0,0,-1,0,1],[0,0,0,-4.5,-3.5,-4.5,-3.5],thick=3.

plots,l_center(s(0)),mean(realiz),psym=8,color=70

xerr,l_mid(s(0)),mean(realiz),halfwidth(s(0)),color=70


; add BICEP2+Keck Oct2015 detections to plot

s=where(experiment eq 'BICEP2+Keck')

oplot,l_center(s(1:*)),bb(s(1:*)),ps=4,symsize=0.6,color=100

;oplot,l_center(s),bb(s),ps=4,symsize=0.6,color=100

errplot,l_center(s(1:*)),bb(s(1:*))+sigma_bb_plus(s(1:*)),$
 bb(s(1:*))-sigma_bb_minus(s(1:*)),width=0,color=100

;errplot,l_center(s),bb(s)+sigma_bb_plus(s),$
; bb(s)-sigma_bb_minus(s),width=0,color=100

l_mid=(l_min+l_max)/2.

halfwidth=l_mid-l_min

xerr,l_mid(s(1:*)),bb(s(1:*)),halfwidth(s(1:*)),color=100

;xerr,l_mid(s),bb(s),halfwidth(s),color=100

; oplot 95% confidence upper limit for first BICEP2+Keck bin
;  where CMB BB = 0

usersym,[-1,1,0,0,-1,0,1],[0,0,0,-4.5,-3.5,-4.5,-3.5],thick=3.

plots,l_center(s(0)),9.63e-03,psym=8,color=100

xerr,l_mid(s(0)),9.63e-03,halfwidth(s(0)),color=100


; add Polarbear detections to plot

s = where(experiment eq 'POLARBEAR' and l_center ne 1500)

oplot,l_center(s),bb(s),ps=4,symsize=0.6,color=20

errplot,l_center(s),bb(s)+sigma_bb_plus(s),bb(s)-sigma_bb_minus(s),width=0,color=20

xerr,l_mid(s),bb(s),halfwidth(s),color=20

; oplot 95% confidence upper limit for Polarbear l_center=1500 bin

s = where(experiment eq 'POLARBEAR' and l_center eq 1500)

realiz=fltarr(10)

for i=0,9 do begin

  get_limit,bb(s(0)),sigma_bb_plus(s(0)),0.95,limit

  realiz(i) = limit

endfor

plots,l_center(s),mean(realiz),psym=8,color=20

xerr,l_mid(s),mean(realiz),halfwidth(s),color=20


; add SPTpol detections to plot

s = where(experiment eq 'SPTpol')

oplot,l_center(s(0:3)),bb(s(0:3)),ps=4,symsize=0.6,color=200

errplot,l_center(s(0:3)),bb(s(0:3))+sigma_bb_plus(s(0:3)),$
 bb(s(0:3))-sigma_bb_minus(s(0:3)),width=0,color=200

xerr,l_mid(s(0:3)),bb(s(0:3)),halfwidth(s(0:3)),color=200

; oplot 95% confidence upper limit for SPTpol l_center=2100 bin

s = where(experiment eq 'SPTpol' and l_center eq 2100)

realiz=fltarr(10)

for i=0,9 do begin

  get_limit,bb(s(0)),sigma_bb_plus(s(0)),0.95,limit

  realiz(i) = limit

endfor

plots,l_center(s),mean(realiz),psym=8,color=200

xerr,l_mid(s),mean(realiz),halfwidth(s),color=200


alt_xyouts,.05,.81,'SPTpol',charsize=0.8,color=200,/log
alt_xyouts,.05,.85,'Polarbear',charsize=0.8,color=20,/log
alt_xyouts,.05,.89,'BICEP2+Keck/Planck',charsize=0.8,color=70,/log
alt_xyouts,.05,.93,'BICEP2+Keck',charsize=0.8,color=100,/log


; oplot LCDM predictions

; read BICEP2 team results for r=0.1

readcol,'B2_3yr_camb_planck_lensed_uK_20140314.txt',$
 l_p1_lens,c_BB_p1_lens,format='I,X,X,X,F',skipline=14

readcol,'B2_3yr_camb_planck_withB_uK_20140314.txt',$
 l_p1_inflation,c_BB_p1_inflation,format='I,X,X,X,F',skipline=14

; get theoretical spectrum for r=0.1, call it c_BB_p1

; start with lensing component

c_BB_p1=c_BB_p1_lens

; add inflationary component

s=indgen(n_elements(l_p1_inflation))

c_BB_p1(s)=c_BB_p1(s)+c_BB_p1_inflation

loadct,0

oplot,l_p1_lens,c_BB_p1  ; theor. curve for r=0.1

oplot,l_p1_inflation,c_BB_p1_inflation,linestyle=2  ; inflationary component, r=0.1

oplot,l_p1_lens,c_BB_p1_lens,linestyle=1   ; lensing component

xyouts,3.,1.3*c_BB_p1_inflation(2),'r=0.1',charsize=0.8


device,/close

!p.thick=1
!p.charthick=1
!x.thick=1
!y.thick=1

loadct,0

end


pro alt_xyouts,xfrac,yfrac,str,charsize=charsize,color=color,log=log

; inputs
; xfrac - print string starting at this fraction of plot range in x
; yfrac -                                                         y
; str   - string to print

if (not keyword_set(log)) then begin

xyouts,!x.crange(0) + xfrac*(!x.crange(1)-!x.crange(0)),$
       !y.crange(0) + yfrac*(!y.crange(1)-!y.crange(0)),$
       str,charsize=charsize,color=color

endif else begin

xyouts,10^(!x.crange(0) + xfrac*(!x.crange(1)-!x.crange(0))),$
       10^(!y.crange(0) + yfrac*(!y.crange(1)-!y.crange(0))),$
       str,charsize=charsize,color=color

endelse

end


pro get_limit,mn,sigma,confidence,upper_limit

; gets upper limit from Gaussian likelihood truncated at zero

; input 
; mn, sigma - mean and sigma of Gaussian likelihood
; confidence - e.g., 0.95 for 95% confidence

like=randomn(seed,100000)

like=like*sigma

like=like+mn

sel=where(like ge 0.)

like=like(sel)

np=n_elements(sel)

frac=fltarr(np)

order=sort(like)

like=like(order)

for i=0L,np-1 do frac(i) = total(like(0:i),/double)/total(like,/double)

upper_limit=like(min(where(frac gt confidence)))

end



pro xerr,x,y,sigx,color=color

; oplots error bars in x


!psym=0

np=n_elements(y)

xmin=x-sigx

xmax=x+sigx

for i=0,np-1 do begin

  oplot,[xmin(i),xmax(i)],[y(i),y(i)],color=color

endfor

end

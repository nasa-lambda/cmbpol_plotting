pro plot_bb_bounds

; Makes a B-mode power spectrum plot showing 95% confidence upper limits
; for B-mode power from different experiments. Data from ACTPol, BICEP1, 
; BOOMERanG, CAPMAP, CBI, DASI, MAXIPOL, QUaD, QUIET-Q, QUIET-W, and WMAP 
; are included.
; For comparison, theoretical curves for a LCDM model with tensor-to-scalar 
; ratio r=0.1 and r=0.01 are also plotted as solid curves. The inflationary
; and gravitational lensing components are plotted separately as dashed
; and dotted curves, respectively.
; The data are plotted to a postscript file BB_bounds.ps.

; The experimental data are read from a file bb_data_2015apr.txt, which
; should be copied to the user's local directory.
; Sources for the data and more information are given in 
; http://lambda.gsfc.nasa.gov/graphics/bb_upperlimits/

; For the QUIET-W results, the larger of the upper limits from the two
; pipelines is plotted for each l-bin.

; The theoretical data are read from files made by the BICEP2 team,
; from their calculations using the March 2013 version of CAMB.
; They should be copied to the user's local directory from
; http://lambda.gsfc.nasa.gov/data/suborbital/BICEP2/B2_3yr_camb_planck_lensed_uK_20140314.txt
; http://lambda.gsfc.nasa.gov/data/suborbital/BICEP2/B2_3yr_camb_planck_withB_uK_20140314.txt

; This code calls readcol.pro and associated routines from the IDL 
; Astronomy User's Library, http://idlastro.gsfc.nasa.gov/


set_plot,'ps'

device,filename='BB_bounds.ps',/color,bits=8

!p.thick=4
!p.charthick=4
!x.thick=4
!y.thick=4

; read and plot 95% confidence upper limits

readcol,'bb_data_2015apr.txt',expt,l_min,l_max,bb_limit,format='A,I,I,F',skipline=35

unique_list,expt,expt_uniq,expt_index

n_uniq = n_elements(expt_uniq)   ; number of experiments with upper limit data

; set up arrays of IDL color tables and colors to
;  use for the different experiments

ctables=[27,12,22,12,4,12,4,12,10,12,12]

colors=[130,20,90,100,240,70,131,120,210,150,200]

delta_y = (.96-.5)/n_uniq   ; separation in y for xyouts

plot_oo,[1.,1.e+4],[1.e-2,1.e+3],ps=3,xtitle='Multipole !8l!3',$
 ytitle='!8l!3(!8l!3+1)C!d!8l!3!u!8BB!n/!32!4p  !3[!4l!3K!u2!n]',$
 yr=[1.e-4,1.e+3],ys=1,/nodata

for i=0,n_uniq-1 do begin

  sel = where(expt eq expt_uniq(i))

  loadct,ctables(i)

  oplot_bb,l_min(sel),l_max(sel),bb_limit(sel),color=colors(i)

  alt_xyouts,.05,.95-i*delta_y,expt_uniq(i),charsize=0.8,color=colors(i),/log

endfor

; oplot lambda CDM predictions

; read BICEP2 team results for r=0.1

readcol,'B2_3yr_camb_planck_lensed_uK_20140314.txt',l_p1_lens,c_BB_p1_lens,format='I,X,X,X,F',skipline=14

readcol,'B2_3yr_camb_planck_withB_uK_20140314.txt',l_p1_inflation,c_BB_p1_inflation,format='I,X,X,X,F',skipline=14

; get theoretical spectrum for r=0.1, call it c_bb_p1

; start with lensing component

c_bb_p1=c_BB_p1_lens

; add inflationary component

s=indgen(n_elements(l_p1_inflation))

c_bb_p1(s)=c_bb_p1(s)+c_BB_p1_inflation

; get theoretical spectrum for r=0.01

c_bb_p01_from_bicep2=c_BB_p1_lens

s=indgen(n_elements(l_p1_inflation))

c_bb_p01_from_bicep2(s)=c_bb_p01_from_bicep2(s)+0.1*c_BB_p1_inflation

loadct,0

oplot,l_p1_lens,c_bb_p1  ; theor. curve for r=0.1

oplot,l_p1_inflation,c_BB_p1_inflation,linestyle=2  ; inflationary component, r=0.1

oplot,l_p1_lens,c_bb_p01_from_bicep2  ; theor. curve for r=0.01

oplot,l_p1_inflation,0.1*c_BB_p1_inflation,linestyle=2   ; inflationary component, r=0.01

oplot,l_p1_lens,c_BB_p1_lens,linestyle=1   ; lensing component

xyouts,3.,1.5*c_BB_p1_inflation(2),'r=0.1',charsize=0.8

xyouts,3.,1.5*c_bb_p01_from_bicep2(2),'r=0.01',charsize=0.8

device,/close

set_plot,'x' 

!p.thick=1
!p.charthick=1
!x.thick=1
!y.thick=1

loadct,0

end


pro oplot_bb,lmin,lmax,bb,color=color

np=n_elements(lmin)

for i=0,np-1 do begin

  oplot,[lmin(i),lmax(i)],[bb(i),bb(i)],ps=0,color=color

endfor

end


;--------------------------------------------------------------------
; Procedure unique_list
;
; IDL procedure to convert a sorted array of elements into a sorted 
; array of unique elements and return the index range of each unique
; element within the original array.  This is useful when grouping
; pixel observations for averaging.
;
; Writen By: BA Franz
;
; Inputs:
;     inlist - input array of elements (must be sorted)
;
; Outputs:
;     outlist - output array of unique elements of inlist
;     index   - a 2xN_Unique array containing the index of the 1st and
;               last occurance of each element of outlist within inlist.
;
;--------------------------------------------------------------------
;
pro unique_list,inlist,outlist,index

nobs = n_elements(inlist)
index = lonarr(2,nobs)
outlist = inlist

i1 = 0L
count = 0L
for i = 1L,nobs-1 do begin
    if (inlist(i) ne inlist(i1)) then begin
        outlist(count) = inlist(i1)
        index(0:1,count) = [i1,i-1]
        count = count+1
        i1 = i
    endif
endfor
index(0:1,count) = [i1,nobs-1]
outlist(count) = inlist(i1)

outlist = outlist(0:count)
index = index(0:1,0:count)

return
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


pro alt_xyouts,xfrac,yfrac,str,charsize=charsize,color=color,log=log

; alternative to IDL's xyouts procedure

; inputs
; xfrac - print the string starting at this fraction of plot range in x
; yfrac -                                                             y
; str   - string to print
; log - set this keyword if plot to be annotated is log-log

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

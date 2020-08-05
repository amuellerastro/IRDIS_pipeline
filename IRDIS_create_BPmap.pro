@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@headfits.pro
@fxposit.pro
@mrd_hread.pro
@get_eso_keyword.pro
@mrdfits.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@robust_mean.pro
@avg.pro
@array_indices.pro
@mwrfits.pro
@fxaddpar.pro
@strnumber.pro
@mpfitfun.pro
@mpfit.pro
@cgcolor.pro
@cggetcolorstate.pro
@cgsnapshot.pro
@cgcolor24.pro
@writefits.pro
@mkhdr.pro
@sxaddpar.pro
@get_date.pro
@daycnv.pro
@sxdelpar.pro
@sxpar.pro
@cgimage.pro
@image_dimensions.pro
@cgdefcharsize.pro
@setdefaultvalue.pro
@cgdefaultcolor.pro
@cgerase.pro
@cgsetcolorstate.pro
@clipscl.pro
@cgresizeimage.pro

function one_gauss, x, p

  z1 = (x-p[1])/p[2]
  g1 = p[0]*exp(-z1^2./2.d0)
  fit = g1

  return, fit

end


function two_gauss, x, p

  z1 = (x-p[1])/p[2]
  z2 = (x-p[4])/p[5]
  g1 = p[0]*exp(-z1^2./2.d0)
  g2 = p[3]*exp(-z2^2./2.d0)
  fit = g1+g2

  return, fit

end

pro IRDIS_create_BPmap

;===============================================================================================

;manually defined detector mask

readcol, '/home/amueller/work/IDLlibs/AO/SPHERE/IRDIS_mask.txt', dx, dy, format='d,d', /silent
dim = 800

;===============================================================================================

;data path and 2 FF with different exposures, e.g. 4s and 10s

readcol, '/home/amueller/work/IDLlibs/AO/SPHERE/datapaths.txt', tmp, format='a', /silent
print, ''
for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]

print, ''
read, 'Select Path: ', selp
path = tmp[selp-1]+'/RAW/IRDIS/'

;===============================================================================================

;file selection

file = file_search(path+'*.fits', count=nfiles)

out = strarr(nfiles)
object = strarr(nfiles)
type = strarr(nfiles)
catg = strarr(nfiles)
opti2 = strarr(nfiles)
dit = strarr(nfiles)
icor = strarr(nfiles)
combind = strarr(nfiles)
filt = strarr(nfiles)

for i=0,nfiles-1 do begin

  pos1 = strpos(file[i], '/', /reverse_search)
  pos2 = strpos(file[i], '.fits', /reverse_search)
  out[i] = strmid(file[i], pos1+1, pos2-pos1-1)

  hdr = headfits(file[i], exten=0, /silent)

  object[i] = strcompress(get_eso_keyword(hdr,'OBJECT'),/rem)
  type[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DPR TYPE'),/rem)
  catg[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DPR CATG'),/rem)
  opti2[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS1 OPTI2 NAME'),/rem)
  dit[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DET SEQ1 DIT'),/rem)
  dit[i] = strmid(dit[i], 0, strpos(dit[i], '.'))
  icor[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS COMB ICOR'),/rem)
  combind[i] = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO INS4 FILT2 NAME'),/rem)	;Assembly for infrared neutral density
  filt[i] = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO INS1 FILT NAME'),/rem)	;IRDIS filter unique ID

endfor

;selection of FF frames

idx = where(catg eq 'CALIB' and type eq 'FLAT,LAMP')
tfile = file[idx]
tout = out[idx]
ttype = type[idx]
tcatg = catg[idx]
tdit = dit[idx]
topti2 = opti2[idx]
ticor = icor[idx]
tcombind = combind[idx]
tfilt = filt[idx]

print, ''
for i=0,n_elements(idx)-1 do begin

  print, strcompress(i+1,/rem), tout[i], tcatg[i], ttype[i], tdit[i], tfilt[i], topti2[i], tcombind[i], ticor[i], format='(a3, a32, a15, a17, a5, a8, a8, a10, a13)'
  print, '----------------------------------------------------------------------------------------------------------------'

endfor

fff = strarr(2)
print, ''
read, 'Select low exptime: ', sel
fff[0] = tfile[sel-1]
read, 'Select high exptime: ', sel
fff[1] = tfile[sel-1]

sel = 0

;selection of DARK frame

idx = where(catg eq 'CALIB' or catg eq '')
tfile = file[idx]
tout = out[idx]
ttype = type[idx]
tcatg = catg[idx]
tdit = dit[idx]
topti2 = opti2[idx]
ticor = icor[idx]
tcombind = combind[idx]
tfilt = filt[idx]

print, ''
for i=0,n_elements(idx)-1 do begin

  print, strcompress(i+1,/rem), tout[i], tcatg[i], ttype[i], tdit[i], tfilt[i], topti2[i], tcombind[i], ticor[i], format='(a3, a32, a15, a17, a5, a8, a8, a10, a13)'
  print, '----------------------------------------------------------------------------------------------------------------'

endfor

; print, ''
; read, 'Master Dark: ', sel
; bgf = tfile[sel-1]


;===============================================================================================

;read in of data
ff1 = mrdfits(fff[0], 0, /silent)
ff2 = mrdfits(fff[1], 0, /silent)

; bg = (mrdfits(bgf, 0, hdrbg, /silent))[*,*,0]

; ditbg = double(get_eso_keyword(hdrbg, 'HIERARCH ESO DET SEQ1 DIT'))

;===============================================================================================
;===============================================================================================
; 
; ;create a bad pixel mask form the dark by fitting 2 Gaussians to the histogram of counts
; 
; bpbg = intarr(2048,1024)
; 
; idxs = sort(bg)
; bgs = bg[idxs]
; bgx = dindgen(n_elements(bgs))
; 
; idx = where(bgs gt 0.d0)
; if (idx[0] ne -1) then begin
; 
;   bgs = bgs[idx]
;   bgx = dindgen(n_elements(bgs))
; 
; endif
; 
; idx = where(bgs lt 5.*median(bgs))
; if (idx[0] ne -1) then begin
; 
;   bgs = bgs[idx]
;   bgx = dindgen(n_elements(bgs))
; 
; endif
; 
; ; window, 1
; 
; hd = histogram(bgs, binsize=0.1, locations=hx)
; 
; ; window, 2
; nbins = floor((max(bgs)-min(bgs))/0.1+1.)
; ; hx = dindgen(nbins)*0.1d0+0.05d0
; 
; if (n_elements(hx) ne nbins) then begin
; 
;   print, ''
;   print, 'Number of bins not computed correclty. Stop.'
;   stop
; 
; endif
; 
; dumerr = n_elements(hx)
; dumerr[*] = 1.
; peak = max(hd, idxmax)
; 
; if (ditbg gt 16.d0) then begin
; 
;   start_val = [peak, hx[idxmax], 3.5, 0.5*peak, 2.*hx[idxmax], 3.5]
; 
;   pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},6)
;   ;force amplitude to be positive
;   pi[0].limited(0) = 1
;   pi[0].limited(1) = 1
;   pi[0].limits(0) = 1.d-15
;   pi[0].limits(1) = 3.*peak
;   ;force fit to find peak at the estimated position
;   pi[1].limited(0) = 1
;   pi[1].limited(1) = 1
;   pi[1].limits(0) = 1.d-15
;   pi[1].limits(1) = 1.5*hx[idxmax]
;   ;do not allow unreasonable large fwhms
;   pi[2].limited(0) = 1
;   pi[2].limited(1) = 1
;   pi[2].limits(0) = 1.d-15
;   pi[2].limits(1) = 20.
;   ;force amplitude to be positive
;   pi[3].limited(0) = 1
;   pi[3].limited(1) = 1
;   pi[3].limits(0) = 1.d-15
;   pi[3].limits(1) = peak
;   ;force fit to find peak at the estimated position
;   pi[4].limited(0) = 1
;   pi[4].limited(1) = 1
;   pi[4].limits(0) = 1.2*hx[idxmax]
;   pi[4].limits(1) = 3.5*hx[idxmax]
;   ;do not allow unreasonable large fwhms
;   pi[5].limited(0) = 1
;   pi[5].limited(1) = 1
;   pi[5].limits(0) = 1.d-15
;   pi[5].limits(1) = 20.
; 
;   result = mpfitfun('two_gauss', hx, hd, dumerr, start_val, weights=sqrt(hd), yfit=yfit, /quiet, parinfo=pi)
; 
;   fwhm1 = 2.d0*sqrt(2.d0*alog(2.d0))*result[2]
;   fwhm2 = 2.d0*sqrt(2.d0*alog(2.d0))*result[5]
; 
;   x1 = abs(result[1]-fwhm1)
;   x2 = abs(result[1]+fwhm1)
;   x3 = abs(result[4]-fwhm2)
;   x4 = abs(result[4]+fwhm2)
;   xr = [x1,x2,x3,x4]
;   ; x0 = floor(min(xr))
;   ; x1 = ceil(max(xr))
;   x0 = min(xr)
;   x1 = max(xr)
; 
; endif else begin
; 
;   start_val = [peak, hx[idxmax], 3.5]
; 
;   pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},3)
;   ;force amplitude to be positive
;   pi[0].limited(0) = 1
;   pi[0].limited(1) = 1
;   pi[0].limits(0) = 1.d-15
;   pi[0].limits(1) = 3.*peak
;   ;force fit to find peak at the estimated position
;   pi[1].limited(0) = 1
;   pi[1].limited(1) = 1
;   pi[1].limits(0) = -10.;1.d-15
;   pi[1].limits(1) = 1.5*hx[idxmax]
;   ;do not allow unreasonable large fwhms
;   pi[2].limited(0) = 1
;   pi[2].limited(1) = 1
;   pi[2].limits(0) = 1.d-15
;   pi[2].limits(1) = 20.
; 
;   result = mpfitfun('one_gauss', hx, hd, dumerr, start_val, weights=sqrt(hd), yfit=yfit, /quiet, parinfo=pi)
; 
;   fwhm1 = 2.d0*sqrt(2.d0*alog(2.d0))*result[2]
; 
; ;   x1 = abs(result[1]-fwhm1)
;   x1 = (result[1]-fwhm1)	;did this for J-band observations. Gaussian around 0
;   x2 = abs(result[1]+fwhm1)
;   xr = [x1,x2]
;   x0 = min(xr)
;   x1 = max(xr)
; 
; endelse
; 
; window, 0, xs=1000, ys=600
; !p.multi = [0,1,2]
; plot, bgx, bgs, psym=3, xst=1, yr=[-5.*median(bgs),5.*median(bgs)], charsize=2, yst=1, thick=3
;   plots, !x.crange, x0*[1,1], color=cgcolor('green')
;   plots, !x.crange, x1*[1,1], color=cgcolor('green')
; ; cghistoplot, bgs, histdata=hd, reverse_indices=r, bin=0.1, charsize=2, axiscolorname=cgcolor('white'), ytitle=''
; plot, hx, hd, xst=1, charsize=2, thick=3
;   oplot, hx, yfit, color=cgcolor('red')
;   g1 = one_gauss(hx, result[0:2])
;   if (ditbg gt 32.d0) then g2 = one_gauss(hx, result[3:5])
;   oplot, hx, g1, color=cgcolor('yellow')
;   if (ditbg gt 32.d0) then oplot, hx, g2, color=cgcolor('yellow')
;   plots, x0*[1,1], !y.crange, color=cgcolor('green')
;   plots, x1*[1,1], !y.crange, color=cgcolor('green')
; 
; !p.multi = [0,1,0]
; 
; print, ''
; print, 'Counts considered: ', x0, x1
; print, ''
; 
; quest = ''
; read, 'Limits OK? y/n: ', quest
; if (quest ne 'y') then begin
; 
;   read, 'Lower Limit: ', x0
;   read, 'Upper Limit: ', x1
; 
;   plots, x0*[1,1], !y.crange, color=cgcolor('red')
;   plots, x1*[1,1], !y.crange, color=cgcolor('red')
; 
; endif
; 
; idx = where(bg lt x0 or bg gt x1)
; 
; idxbad = array_indices(bpbg, idx)
; for i=0L,n_elements(idxbad[0,*])-1 do bpbg[idxbad[0,i], idxbad[1,i]] = 1

;===============================================================================================
;===============================================================================================

;only consider 800x800pixel of left and right detector respectively

ff1l = ff1[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]
ff2l = ff2[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]
ff1r = ff1[1024:2047,*]
ff2r = ff2[1024:2047,*]
ff1r = ff1r[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]
ff2r = ff2r[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]

;===============================================================================================

;apply detector mask
; 
; mask = dblarr(2048,1024)
; mask[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1] = 1.
; mask[dx[1]-dim/2+1024:dx[1]+dim/2-1+1024, dy[1]-dim/2:dy[1]+dim/2-1] = 1.
; 
; ff1 = ff1*mask
; ff2 = ff2*mask

;===============================================================================================

;define some arrays
; bp1 = intarr(2048,1024)
; bp1[*,*] = 1
; bp2 = intarr(2048,1024)
; bp2[*,*] = 1
; bp = intarr(2048,1024)
bp1l = intarr(800,800)
bp1l[*,*] = 1
bp1r = bp1l
bp2l = intarr(800,800)
bp2l[*,*] = 1
bp2r = bp1l
bp = intarr(2048,1024)
bp[*,*] = 1

;compute ratio of both FFs and reject outliers
; ratio = ff2/ff1
ratiol = ff2l/ff1l
ratior = ff2r/ff1r

rml = robust_mean(ratiol, 4., sigma, numrej, goodind=goodl)
rmr = robust_mean(ratior, 4., sigma, numrej, goodind=goodr)

; window, 0, xs=1500, ys=800
; !p.multi=[0,1,2]
; 
;   plot, ratio, xst=1, yst=1, yr=[-20,20], psym=3
;   plot, ratio[good], psym=3, xst=1, yst=1, yr=[-20,20]
; 
; !p.multi=[0,1,0]

; idxgoodl = array_indices(bp1l, goodl)
; for i=0L,n_elements(idxgoodl[0,*])-1 do bp1l[idxgoodl[0,i], idxgoodl[1,i]] = 0
bp1l[*,*] = 0

; idxgoodr = array_indices(bp1r, goodr)
; for i=0L,n_elements(idxgoodr[0,*])-1 do bp1r[idxgoodr[0,i], idxgoodr[1,i]] = 0
bp1r[*,*] = 0

;make a selection of bad pixels (in this case pixels which do not tget illuminated) based on flux
rml = robust_mean(ff2l, 4., sigma, numrej, goodind=goodl)

idxgoodl = array_indices(bp2l, goodl)
for i=0L,n_elements(idxgoodl[0,*])-1 do bp2l[idxgoodl[0,i], idxgoodl[1,i]] = 0

rmr = robust_mean(ff2r, 4., sigma, numrej, goodind=goodr)
idxgoodr = array_indices(bp2r, goodr)
for i=0L,n_elements(idxgoodr[0,*])-1 do bp2r[idxgoodr[0,i], idxgoodr[1,i]] = 0

;combine both BP maps

bpl = bp1l+bp2l
idx0 = where(bpl gt 0.)
idx = array_indices(bpl, idx0)
for i=0L,n_elements(idx[0,*])-1 do bpl[idx[0,i], idx[1,i]] = 1

bpr = bp1r+bp2r
idx0 = where(bpr gt 0.)
idx = array_indices(bpr, idx0)
for i=0L,n_elements(idx[0,*])-1 do bpr[idx[0,i], idx[1,i]] = 1

;===============================================================================================

;put sub frames back together

; bp = intarr(2048,1024)
; bp[*,*] = 1
bp[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1] = bpl
bp[1024+dx[1]-dim/2:1024+dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1] = bpr

;put BP map of FF and Dark together

bp_final = bp;bpbg+
idx = where(bp_final gt 0)
idxbad = array_indices(bp_final, idx)
for i=0L,n_elements(idxbad[0,*])-1 do bp_final[idxbad[0,i], idxbad[1,i]] = 1

; mwrfits, bp_final, path+'static_badpixels.fits'
writefits, path+'static_badpixels.fits', bp_final

;===============================================================================================

window, 2, xs=1024, ys=512
cgimage, bp_final, stretch=2

print, ''
print, 'Percantage of bad Pixel: ', 100.*double(n_elements(where(bp_final eq 1.)))/double(n_elements(bp_final))
print, ''

stop
end
@linspace.pro
@get_eso_keyword.pro
@strc.pro
@fixpix_mod.pro
@showsym.pro
@mrdfits.pro
@fxposit.pro
@mrd_hread.pro
@fxpar.pro
@gettok.pro
@valid_num.pro
@mrd_skip.pro
@dist_circle.pro
@proceeding_text.pro
@mwrfits.pro
@fxaddpar.pro
@detabify.pro
@fxparpos.pro
@cgimage.pro
@image_dimensions.pro
@cgdefcharsize.pro
@setdefaultvalue.pro
@cgdefaultcolor.pro
@cggetcolorstate.pro
@cgerase.pro
@cgsetcolorstate.pro
@cgcolor.pro
@strsplit.pro
@cgsnapshot.pro
@cgcolor24.pro
@sdevscl.pro
@cgscalevector.pro
@fpufix.pro
@cgresizeimage.pro
@congrid.pro
@cgplot.pro
@cgbitget.pro
@convert_to_type.pro
@cgcheckforsymbols.pro
@colorsareidentical.pro
@array_indices.pro
@tvcircle.pro
@cgplots.pro
@mpfit2dpeak.pro
@mpfit.pro
@mpfit2dfun.pro
@mpfitfun.pro
@sixlin.pro
@fxmove.pro
@match.pro
@mrd_struct.pro
@sxaddhist.pro
@get_screen_size.pro
@headfits.pro
@la_cosmic_MOD.pro
@reverse.pro
@fxhmake.pro
@fxwrite.pro
@check_fits.pro
@get_date.pro
@daycnv.pro
@host_to_ieee.pro
@readcol.pro
@remchar.pro
@strnumber.pro
@writefits.pro
@mkhdr.pro
@sxaddpar.pro
@get_date.pro
@daycnv.pro
@sxdelpar.pro
@sxpar.pro
@scale_image_am.pro

function moffat_slope, x, y, p

;   u = ( (x-p[3])/p[1] )^2. + ( (y-p[4])/p[2] )^2.
;   moffat = p[0]/(u + 1.d0)^p[5]
;   moffat = mpfit2dpeak_moffat_mod(x, y, p)

  widx  = abs(p[1]) > 1e-20 & widy  = abs(p[2]) > 1e-20 
  xp    = x-p[3]            & yp    = y-p[4]
  theta = p[5]
  c  = cos(theta) & s  = sin(theta)
  u = ( (xp * (c/widx) - yp * (s/widx))^2 + $
                (xp * (s/widy) + yp * (c/widy))^2 )

  moffat = p[0] / (u + 1)^p[6]

  slope = p[7]*x+p[8]*y

  fit = moffat+slope

  return, fit

end

; ; Moffat Function
; function mpfit2dpeak_moffat_mod, x, y, p
; 
;   u = mpfit2dpeak_u_mod(x, y, p)
;   return, p[0] / (u + 1)^p[6]
; 
; end
; 
; function mpfit2dpeak_u_mod, x, y, p
; 
;   widx  = abs(p[1]) > 1e-20 & widy  = abs(p[2]) > 1e-20 
;   xp    = x-p[3]            & yp    = y-p[4]
;   theta = p[5]
;   ;if theta NE 0 then begin
;       c  = cos(theta) & s  = sin(theta)
;       return, ( (xp * (c/widx) - yp * (s/widx))^2 + $
;                 (xp * (s/widy) + yp * (c/widy))^2 )
; ;   endif else begin
; ;       return, (xp/widx)^2 + (yp/widy)^2
; ;   endelse
; 
; end


pro IRDIS_star_center

;===============================================================================================

readcol, '/home/amueller/work/IDLlibs/AO/SPHERE/IRDIS_mask.txt', dx, dy, format='d,d', /silent

;====================================================================================================

readcol, '/home/amueller/work/IDLlibs/AO/SPHERE/datapaths.txt', tmp, format='a', /silent
print, ''
for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]

print, ''
read, 'Select Path: ', selp
path = tmp[selp-1]+'/RAW/IRDIS/'

;====================================================================================================

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

;selection of STAR CENTER frame

idx = where(catg eq 'SCIENCE' and type eq 'OBJECT,CENTER')
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

print, ''
sel = intarr(2)
read, 'Star Center (first and last number): ', sel
scf = tfile[indgen(sel[1]-sel[0]+1)+sel[0]-1]
sel = 0

;selection of CALIBRATION frame

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

print, ''
; read, 'Master Dark: ', sel
; bgf = tfile[sel-1]
; read, 'Bad Pixel Mask: ', sel
; bpf = tfile[sel-1]
bpf = path+'static_badpixels.fits'
; read, 'Master Flat: ', sel
; fff = tfile[sel-1]
fff = path+'instrument_flat.fits'
; read, 'Sky BG: ', sel
; skf = tfile[sel-1]
match = strmatch(tfile, '*sky_background*')
idx = where(match eq 1)
skf = tfile[idx]

distquest = ''
; read, 'With distortion correction? y/n ', distquest

;====================================================================================================

;cut images down to 800x800px
;using bad pixel map of FF too and all the edges are full with bad pixels -> takes a very long time for fixpix_mod
; xc0 = 512
; yc0 = 512
dim = 800

;====================================================================================================

screendim = get_screen_size()

; if (distquest eq 'y') then begin
; 
;   restore, path+'distortion_coeffs_left.sav', /v
;   kxl = coeffl.kx
;   kyl = coeffl.ky
; 
;   restore, path+'distortion_coeffs_right.sav',/v
;   kxr = coeffr.kx
;   kyr = coeffr.ky
; 
; endif


; bg = mrdfits(bgf, 0, /silent)
ff = mrdfits(fff, 0, /silent)
bp = mrdfits(bpf, 0, /silent)
sk = mrdfits(skf, 0, hdrsk, /silent)

; bgl = bg[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1,*]
ffl = ff[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1,*]
bpl = bp[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1,*]
skl = sk[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]

; bgr = bg[1024:2047,*]
ffr = ff[1024:2047,*]
bpr = bp[1024:2047,*]
skr = sk[1024:2047,*]

; bgr = bgr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1,*]
ffr = ffr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1,*]
bpr = bpr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1,*]
skr = skr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]

;====================================================================================================

nfiles = n_elements(scf)	;number of files/cubes

;dummy read in to get dimensions
tmp = mrdfits(scf[0], 0, hdrsc, /silent)
sz = size(tmp)

if (n_elements(sz) eq 6) then begin

  sc = dblarr(2048,1024,sz[3],nfiles) 
  nim = sz[3]	;number of images per cube

endif else begin

  sc = dblarr(2048,1024,1,nfiles)
  nim = 1

endelse

pacx = dblarr(nfiles)
pacy = dblarr(nfiles)
date = strarr(nfiles)
jd = dblarr(nfiles)

for i=0,nfiles-1 do begin

  tmp = mrdfits(scf[i], 0, hdrsc, /silent)
  sc[*,*,*,i] = tmp

  pacx[i] = double(get_eso_keyword(hdrsc, 'HIERARCH ESO INS1 PAC X'))/18.d0
  pacy[i] = double(get_eso_keyword(hdrsc, 'HIERARCH ESO INS1 PAC Y'))/18.d0
  date[i] = get_eso_keyword(hdrsc, 'DATE-OBS')
  jd[i] = double(get_eso_keyword(hdrsc, 'MJD-OBS'))+2400000.5d0
  if (i eq 0) then dit = double(get_eso_keyword(hdrsc,'HIERARCH ESO DET SEQ1 DIT'))

  exptime = double(get_eso_keyword(hdrsc, 'HIERARCH ESO DET SEQ1 EXPTIME'))
  ndit = double(get_eso_keyword(hdrsc, 'HIERARCH ESO DET NDIT'))

endfor

scl = sc[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1,*,*]

scr = sc[1024:2047,*,*,*]
scr = scr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1,*,*]

;====================================================================================================

;bad pixel correction of images and division by FF

sccl = scl
sccr = scr

newjd = dblarr(nim,nfiles)

if (dit gt 10.) then readn = 10. else readn = 4.	;DIT of 10s chosen abritrarily

for i=0,nfiles-1 do begin

  print, ''
  print, 'Reducing file '+strcompress(i+1,/rem)+' / '+strcompress(nfiles,/rem)
  print, ''

  for j=0,nim-1 do begin

    print, ''
    print, 'Reducing frame '+strcompress(j+1,/rem)+' / '+strcompress(nim,/rem)
    print, ''

    newjd[j,i] = jd[i] + ((exptime/ndit)*j + (exptime/ndit)/2.)/86400.d0/2.d0

    tmp = (sccl[*,*,j,i]-skl)/ffl
    fixpix_mod, tmp, bpl, outim, npix=24, /weight;, /silent
    sccl[*,*,j,i] = outim
;     la_cosmic_MOD, sccl[*,*,j,i], outim, path, gain=1.75, readn=readn, skyval=median(skl,/even)
;     sccl[*,*,j,i] = outim

    tmp = (sccr[*,*,j,i]-skr)/ffr
    fixpix_mod, tmp, bpr, outim, npix=24, /weight;, /silent
    sccr[*,*,j,i] = outim
;     la_cosmic_MOD, sccr[*,*,j,i], outim, path, gain=1.75, readn=readn, skyval=median(skr,/even)
; ;     sccr[*,*,j,i] = outim

;     if (distquest eq 'y') then begin
; 
;       sccl[*,*,j,i] = poly_2d(sccl[*,*,j,i], kxl, kyl, 2, cubic=-0.5)
;       sccr[*,*,j,i] = poly_2d(sccr[*,*,j,i], kxr, kyr, 2, cubic=-0.5)
; 
;     endif

  endfor

endfor

;====================================================================================================
;====================================================================================================

;find center LEFT

window, 0, xs=dim+200, ys=dim+200
device, cursor_standard=2
cgimage, scale_image_am(sccl[*,*,0,0]), /axis;, minvalue=0., maxvalue=3000.;, stretch=10

!mouse.button = 0
print, ''
print, 'Select 4 Stars (clockwise)'

radius = 7.
i = 0
xc = dblarr(4)
yc = dblarr(4)
mx = dblarr(4)
cutim = dblarr(2.*radius+1, 2.*radius+1, 4, nim, nfiles)	;selected sub image to fit position of object

;manually select planet candidate
;while !mouse.button ne 4 do begin
repeat begin

  cursor, x, y, 3, /data

  if (!mouse.button eq 1) then begin

    !mouse.button = 0

    x = ceil(x) & y=ceil(y)

    ;refine centering by identifying max intensity
    tmp_cutim = sccl[x-radius:x+radius, y-radius:y+radius,0,0]

    mx[i] = max(tmp_cutim, location)
    idxmax = array_indices(tmp_cutim, location)

    xc[i] = x+(idxmax[0]-radius)
    yc[i] = y+(idxmax[1]-radius)

    ;final cutted region
    for xx=0,nfiles-1 do begin

      for yy=0,nim-1 do begin

	cutim[*,*,i,yy,xx] = sccl[xc[i]-radius:xc[i]+radius, yc[i]-radius:yc[i]+radius,yy,xx]

      endfor

    endfor

;cutim[*,*,i] = im[xc[i]-radius:xc[i]+radius, yc[i]-radius:yc[i]+radius]

    tvcircle, /data, 5, xc[i], yc[i], color=cgcolor('green')
    ;draw a box corresponding of the slected radius
    x0 = xc[i]-radius
    x1 = xc[i]+radius
    y0 = yc[i]-radius
    y1 = yc[i]+radius
    plots, [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], /data, color=cgcolor('green'), thick=3
    xyouts, x0, y1, strcompress(i+1), /data, alignment=1, color=cgcolor('green')

    ;window, 1, xs=10.*csz, ys=10.*csz
    ;cgimage, cutim[*,*,i], minvalue=stretch[0], maxvalue=stretch[1]

    i = i+1

  endif

endrep until (i eq 4)


print, ''
print, 'Fitting of Star Positions'

xposl = dblarr(4, nim, nfiles) & xposlerr = xposl
yposl = dblarr(4, nim, nfiles) & yposlerr = yposl

for i=0,3 do begin

  for j=0,nfiles-1 do begin

    for k=0,nim-1 do begin

;       estimates = [median(cutim[*,*,i,k,j]), mx[i]-median(cutim[*,*,i,k,j]), 4., 4., radius, radius, 0., 1.]
;       xa = dindgen(2*radius+1)+1.d0 & ya = xa

      ;weights as a function of distance
;       dist_circle, weight, 2.*radius+1
;       weights = 1./sqrt(weight+1.)


;       weights = cutim[*,*,i,k,j]
;       idx1 = where(cutim[*,*,i,k,j] eq 0.)	;e.g. beta Pic close to the center which is masked out
;       if (idx1[0] ne -1) then begin
; 	idx2 = array_indices(weights, idx1)
; 	for l=0,n_elements(idx2[0,*])-1 do weights[idx2[0,l],idx2[1,l]] = 1.d5
;       endif
;       sign = signum(cutim[*,*,i,k,j])
;       idx1 = where(sign eq -1.)
;       if (idx1[0] ne -1) then begin
; 	idx2 = array_indices(weights, idx1)
; 	for l=0,n_elements(idx2[0,*])-1 do weights[idx2[0,l],idx2[1,l]] = 1.d5
;       endif
;       weights = 1./sqrt(weights)


      xa = (dindgen(2.*radius+1.)+1.d0)#(dblarr(2.*radius+1.)+1.)
      ya = (dblarr(2.*radius+1.)+1.)#(dindgen(2.*radius+1.)+1.d0)

      dumerr = xa
      dumerr[*,*] = 1.

      ;		amplitude		fwhm	fwhm	xc	yc	rot	power	slope x,y
      estimates = [mx[i]-median(cutim[*,*,i,k,j]), 2., 2., radius, radius, 0.1, 1., 1., 1.]

      pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},9)
      pi[1].limited(0) = 1
      pi[1].limits(0) = 1.0d0
      pi[1].limited(1) = 1
      pi[1].limits(1) = 4.0d0
      pi[2].limited(0) = 1
      pi[2].limits(0) = 1.0d0
      pi[2].limited(1) = 1
      pi[2].limits(1) = 4.0d0
      pi[3].limited(0) = 1
      pi[3].limits(0) = radius-3.d0
      pi[3].limited(1) = 1
      pi[3].limits(1) = radius+3.d0
      pi[4].limited(0) = 1
      pi[4].limits(0) = radius-3.d0
      pi[4].limited(1) = 1
      pi[4].limits(1) = radius+3.d0
      pi[5].limited(0) = 1
      pi[5].limits(0) = -2.*!DPI
      pi[5].limited(1) = 1
      pi[5].limits(1) = 2.*!DPI

      A = mpfit2dfun('moffat_slope', xa, ya, cutim[*,*,i,k,j], dumerr, estimates, dof=dof, chisq=chisq, perror=perror, yfit=yfit, parinfo=pi, /quiet);, weights=weights)

      xposl[i,k,j] = A[3]-radius+xc[i];+1
      yposl[i,k,j] = A[4]-radius+yc[i];-1
      xposlerr[i,k,j] = perror[3]
      yposlerr[i,k,j] = perror[4]

;       pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},8)
;       pi[0].limited(0) = 1
;       pi[0].limits(0) = 0.d0
;       pi[0].limited(1) = 1
;       pi[0].limits(1) = 5.*estimates[0]
;       pi[2].limited(0) = 1
;       pi[2].limits(0) = 2.0d0
;       pi[2].limited(1) = 1
;       pi[2].limits(1) = 6.0d0
;       pi[3].limited(0) = 1
;       pi[3].limits(0) = 2.0d0
;       pi[3].limited(1) = 1
;       pi[3].limits(1) = 6.0d0
;       pi[4].limited(0) = 1
;       pi[4].limits(0) = radius-3.d0
;       pi[4].limited(1) = 1
;       pi[4].limits(1) = radius+3.d0
;       pi[5].limited(0) = 1
;       pi[5].limits(0) = radius-3.d0
;       pi[5].limited(1) = 1
;       pi[5].limits(1) = radius+3.d0
; 
;       yfit = mpfit2dpeak(cutim[*,*,i,k,j], A, xa, ya, /moffat, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /tilt, /quiet, weights=weight, parinfo=pi)
; 
;       xposl[i,k,j] = A[4]-radius+xc[i];+1
;       yposl[i,k,j] = A[5]-radius+yc[i];-1

;       window, 4, xs=400, ys=400, title='Original'
;       cgimage, cutim[*,*,i,k,j], minvalue=min(cutim[*,*,i,k,j]), maxvalue=max(cutim[*,*,i,k,j])
;       window, 6, xs=400, ys=400, title='Fit'
;       cgimage, yfit, minvalue=min(yfit), maxvalue=max(yfit)
;       window, 5, xs=400, ys=400, title='RMS'
;       cgimage, cutim[*,*,i,k,j]-yfit, minvalue=min(cutim[*,*,i,k,j]-yfit), maxvalue=max(cutim[*,*,i,k,j]-yfit)
;       print, A, max(cutim[*,*,i,k,j])-median(cutim[*,*,i,k,j])
;       hak

    endfor

  endfor

  proceeding_text,loop=(4), i=i, prompt='> Fitting Star Position   '+string(i+1,form='(I4)')

endfor

window, 1, xs=dim/2, ys=dim/2, title='Left Detector', xpos=screendim[0]-800, ypos=0
plot, xposl[0,*], yposl[0,*], psym=1, /yn, xr=[floor(min(xposl))-10,ceil(max(xposl))+10], yr=[floor(min(yposl))-10,ceil(max(yposl))+10], xst=1, yst=1
oplot, xposl[1,*], yposl[1,*], psym=1
oplot, xposl[2,*], yposl[2,*], psym=1
oplot, xposl[3,*], yposl[3,*], psym=1

;linear fit through opposite images

posl = dblarr(2,nim,nfiles)

for i=0,nfiles-1 do begin

  for j=0,nim-1 do begin

    sixlin, [xposl[0,j,i], xposl[2,j,i]], [yposl[0,j,i], yposl[2,j,i]], n1, dn1, m1, dm1
    sixlin, [xposl[1,j,i], xposl[3,j,i]], [yposl[1,j,i], yposl[3,j,i]], n2, dn2, m2, dm2

    ;compute intersection = star center

    posl[0,j,i] = (n2[0]-n1[0])/(m1[0]-m2[0])
    posl[1,j,i] = m1[0]*posl[0,j,i]+n1[0]

    xal1 = linspace(xposl[0,j,i], xposl[2,j,i], 10)
    yal1 = linspace(yposl[0,j,i], yposl[2,j,i], 10)
    linl1 = m1[0]*xal1+n1[0]
    oplot, xal1, linl1
    xal2 = linspace(xposl[1,j,i], xposl[3,j,i], 10)
    yal2 = linspace(yposl[1,j,i], yposl[3,j,i], 10)
    linl2 = m2[0]*xal2+n2[0]
    oplot, xal2, linl2

    ;bring coordinates back to original detector size
;     posl[0,j,i] = posl[0,j,i]+112.d0-pacx[i]
;     posl[1,j,i] = posl[1,j,i]+112.d0-pacy[i]
    posl[0,j,i] = posl[0,j,i]+dx[0]-dim/2.-pacx[i]
    posl[1,j,i] = posl[1,j,i]+dy[0]-dim/2.-pacy[i]

  endfor

endfor



;====================================================================================================
;====================================================================================================

;find center RIGHT

window, 0, xs=dim, ys=dim
device, cursor_standard=2
cgimage, scale_image_am(sccr[*,*,0,0]), /axis;, minvalue=0., maxvalue=3000.;, stretch=10

!mouse.button = 0
print, ''
print, 'Select 4 Stars (clockwise)'

;radius = 5.
i = 0
xc = dblarr(4)
yc = dblarr(4)
mx = dblarr(4)
cutim = dblarr(2.*radius+1, 2.*radius+1, 4, nim, nfiles)	;selected sub image to fit position of object

;manually select planet candidate
;while !mouse.button ne 4 do begin
repeat begin

  cursor, x, y, 3, /data

  if (!mouse.button eq 1) then begin

    !mouse.button = 0

    x = ceil(x) & y=ceil(y)

    ;refine centering by identifying max intensity
    tmp_cutim = sccr[x-radius:x+radius, y-radius:y+radius,0,0]

    mx[i] = max(tmp_cutim, location)
    idxmax = array_indices(tmp_cutim, location)

    xc[i] = x+(idxmax[0]-radius)
    yc[i] = y+(idxmax[1]-radius)

    ;final cutted region
    for xx=0,nfiles-1 do begin

      for yy=0,nim-1 do begin

	cutim[*,*,i,yy,xx] = sccr[xc[i]-radius:xc[i]+radius, yc[i]-radius:yc[i]+radius,yy,xx]

      endfor

    endfor

;cutim[*,*,i] = im[xc[i]-radius:xc[i]+radius, yc[i]-radius:yc[i]+radius]

    tvcircle, /data, 5, xc[i], yc[i], color=cgcolor('green')
    ;draw a box corresponding of the slected radius
    x0 = xc[i]-radius
    x1 = xc[i]+radius
    y0 = yc[i]-radius
    y1 = yc[i]+radius
    plots, [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], /data, color=cgcolor('green'), thick=3
    xyouts, x0, y1, strcompress(i+1), /data, alignment=1, color=cgcolor('green')

    ;window, 1, xs=10.*csz, ys=10.*csz
    ;cgimage, cutim[*,*,i], minvalue=stretch[0], maxvalue=stretch[1]

    i = i+1

  endif

endrep until (i eq 4)


print, ''
print, 'Fitting of Star Positions'

xposr = dblarr(4, nim, nfiles) & xposrerr = xposr
yposr = dblarr(4, nim, nfiles) & yposrerr = yposr

for i=0,3 do begin

  for j=0,nfiles-1 do begin

    for k=0,nim-1 do begin

;       estimates = [median(cutim[*,*,i,k,j]), mx[i], 4., 4., radius, radius, 0., 1.]
;       xa = dindgen(2*radius+1)+1.d0 & ya = xa

      ;weights as a function of distance
;       dist_circle, weight, 2.*radius+1
;       weight = 1./sqrt(weight+1.)

      ;weights of unity if not NAN
  ;weight = cutim[*,*,i,j]
  ;tmp = cutim[*,*,i,j]	;remove weights
  ;for uu=0,n_elements(cutim[0,*,i,j])-1 do begin
  ; 
  ;for vv=0,n_elements(cutim[0,*,i,j])-1 do begin
  ; 
  ; 	if (finite(cutim[uu,vv,i,j]) ne 1) then begin
  ; 
  ; 	  weight[uu,vv] = 0.
  ; 	  tmp[uu,vv] = 0.d0
  ; 
  ; 	endif else begin
  ; 
  ; 	  tmp[uu,vv] = cutim[uu,vv,i,j]
  ; 	  weight[uu,vv] = 1.
  ; 
  ; 	endelse
  ; 
  ;endfor
  ; 
  ;endfor


      xa = (dindgen(2*radius+1)+1.d0)#(dblarr(2.*radius+1)+1)
      ya = (dblarr(2.*radius+1)+1)#(dindgen(2*radius+1)+1.d0)

      dumerr = xa
      dumerr[*,*] = 1.

      ;		amplitude		fwhm	fwhm	xc	yc	rot	power	slope x,y
      estimates = [mx[i]-median(cutim[*,*,i,k,j]), 2., 2., radius, radius, 0.1, 1., 1., 1.]

      pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},9)
      pi[1].limited(0) = 1
      pi[1].limits(0) = 1.0d0
      pi[1].limited(1) = 1
      pi[1].limits(1) = 4.0d0
      pi[2].limited(0) = 1
      pi[2].limits(0) = 1.0d0
      pi[2].limited(1) = 1
      pi[2].limits(1) = 4.0d0
      pi[3].limited(0) = 1
      pi[3].limits(0) = radius-3.d0
      pi[3].limited(1) = 1
      pi[3].limits(1) = radius+3.d0
      pi[4].limited(0) = 1
      pi[4].limits(0) = radius-3.d0
      pi[4].limited(1) = 1
      pi[4].limits(1) = radius+3.d0
      pi[5].limited(0) = 1
      pi[5].limits(0) = -2.*!DPI
      pi[5].limited(1) = 1
      pi[5].limits(1) = 2.*!DPI

      A = mpfit2dfun('moffat_slope', xa, ya, cutim[*,*,i,k,j], dumerr, estimates, dof=dof, chisq=chisq, perror=perror, yfit=yfit, parinfo=pi, weights=weight, /quiet)

      xposr[i,k,j] = A[3]-radius+xc[i];+1
      yposr[i,k,j] = A[4]-radius+yc[i];-1
      xposrerr[i,k,j] = perror[3]
      yposrerr[i,k,j] = perror[4]

;       yfit = mpfit2dpeak(cutim[*,*,i,k,j], A, xa, ya, /moffat, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /tilt, /quiet, weights=weight)
; 
;       xposr[i,k,j] = A[4]-radius+xc[i];+1
;       yposr[i,k,j] = A[5]-radius+yc[i];-1

;       window, 4, xs=400, ys=400, title='Original'
;       cgimage, cutim[*,*,i,k,j], minvalue=min(cutim[*,*,i,k,j]), maxvalue=max(cutim[*,*,i,k,j])
;       window, 6, xs=400, ys=400, title='Fit'
;       cgimage, yfit, minvalue=min(yfit), maxvalue=max(yfit)
;       window, 5, xs=400, ys=400, title='RMS'
;       cgimage, cutim[*,*,i,k,j]-yfit, minvalue=min(cutim[*,*,i,k,j]-yfit), maxvalue=max(cutim[*,*,i,k,j]-yfit)
;       print, A, max(cutim[*,*,i,k,j])-median(cutim[*,*,i,k,j])
;       hak

    endfor

  endfor

  proceeding_text,loop=(4), i=i, prompt='> Fitting Star Position   '+string(i+1,form='(I4)')

endfor

window, 2, xs=dim/2, ys=dim/2, title='Right Detector', xpos=screendim[0]-dim/2, ypos=0
plot, xposr[0,*], yposr[0,*], psym=1, /yn, xr=[floor(min(xposr))-10,ceil(max(xposr))+10], yr=[floor(min(yposr))-10,ceil(max(yposr))+10], xst=1, yst=1
oplot, xposr[1,*], yposr[1,*], psym=1
oplot, xposr[2,*], yposr[2,*], psym=1
oplot, xposr[3,*], yposr[3,*], psym=1

;linear fit through opposite images

posr = dblarr(2,nim,nfiles)

for i=0,nfiles-1 do begin

  for j=0,nim-1 do begin

    sixlin, [xposr[0,j,i], xposr[2,j,i]], [yposr[0,j,i], yposr[2,j,i]], n1, dn1, m1, dm1
    sixlin, [xposr[1,j,i], xposr[3,j,i]], [yposr[1,j,i], yposr[3,j,i]], n2, dn2, m2, dm2

    ;compute intersection = star center

    posr[0,j,i] = (n2[0]-n1[0])/(m1[0]-m2[0])
    posr[1,j,i] = m1[0]*posr[0,j,i]+n1[0]

    xar1 = linspace(xposr[0,j,i], xposr[2,j,i], 10)
    yar1 = linspace(yposr[0,j,i], yposr[2,j,i], 10)
    linr1 = m1[0]*xar1+n1[0]
    oplot, xar1, linr1
    xar2 = linspace(xposr[1,j,i], xposr[3,j,i], 10)
    yar2 = linspace(yposr[1,j,i], yposr[3,j,i], 10)
    linr2 = m2[0]*xar2+n2[0]
    oplot, xar2, linr2

    ;bring coordinates back to original detector size
;    posr[0,j,i] = posr[0,j,i]+112.d0-pacx[i]
;    posr[1,j,i] = posr[1,j,i]+112.d0-pacy[i]
    posr[0,j,i] = posr[0,j,i]+dx[1]-dim/2.-pacx[i]
    posr[1,j,i] = posr[1,j,i]+dy[1]-dim/2.-pacy[i]

  endfor

endfor


;====================================================================================================

;average values per cube, othwerwise a large trend could be artifically caused when extrapolating later resulting in a shift of the target

if (nim gt 1) then begin

  finposl = dblarr(2,nfiles)
  finposr = dblarr(2,nfiles)
  finnewjd = dblarr(nfiles)

  for i=0,nfiles-1 do begin

    finposl[0,i] = mean(posl[0,*,i])
    finposl[1,i] = mean(posl[1,*,i])

    finposr[0,i] = mean(posr[0,*,i])
    finposr[1,i] = mean(posr[1,*,i])

    finnewjd[i] = mean(newjd[*,i])

  endfor

endif else begin

  finposl = reform(posl)
  finposr = reform(posr)
  finnewjd = reform(newjd)

endelse



;====================================================================================================

;save result

st = replicate({jd:0.d0, center_left_x:0.d0, center_left_y:0.d0, center_right_x:0.d0, center_right_y:0.d0},nfiles)
for i=0,nfiles-1 do begin

  st[i].jd = finnewjd[i]
  st[i].center_left_x = reform(finposl[0,i])-1.d0
  st[i].center_left_y = reform(finposl[1,i])-1.d0
  st[i].center_right_x = reform(finposr[0,i])-1.d0
  st[i].center_right_y = reform(finposr[1,i])-1.d0

endfor

; st = replicate({jd:0.d0, center_left_x:0.d0, center_left_y:0.d0, center_right_x:0.d0, center_right_y:0.d0},nfiles*nim)
; for i=0,nfiles-1 do begin
; 
;   st[i*nim:i*nim+nim-1].jd = finnewjd[i]
;   st[i*nim:i*nim+nim-1].center_left_x = reform(finposl[0,i])-1.d0
;   st[i*nim:i*nim+nim-1].center_left_y = reform(finposl[1,i])-1.d0
;   st[i*nim:i*nim+nim-1].center_right_x = reform(finposr[0,i])-1.d0
;   st[i*nim:i*nim+nim-1].center_right_y = reform(finposr[1,i])-1.d0
; 
; endfor


mwrfits, st, path+'star_center.fits'
;writefits, path+'star_center.fits', st

;====================================================================================================

window, 0, xs=600, ys=600
!p.multi=[0,1,2]
plot, finposl[0,*], finposl[1,*], psym=2, /yn, charsize=2;, color=cgcolor('green')
plot, finposr[0,*], finposr[1,*], psym=2, /yn, charsize=2;, color=cgcolor('blue')
!p.multi=[0,1,0]

stop	
end

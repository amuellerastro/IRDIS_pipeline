@showsym.pro
@la_cosmic_MOD.pro
@ts_diff.pro
@headfits.pro
@fxposit.pro
@mrd_hread.pro
@get_eso_keyword.pro
@readcol.pro
@remchar.pro
@gettok.pro
@strnumber.pro
@readfits.pro
@sxpar.pro
@valid_num.pro
@mrdfits.pro
@fxpar.pro
@mrd_skip.pro
@fixpix_mod.pro
@dist_circle.pro
@cgloadct.pro
@cgimage.pro
@image_dimensions.pro
@cgdefcharsize.pro
@setdefaultvalue.pro
@cgdefaultcolor.pro
@cggetcolorstate.pro
@cgerase.pro
@cgsetcolorstate.pro
@cgcolor.pro
@cgsnapshot.pro
@cgcolor24.pro
@clipscl.pro
@cgscalevector.pro
@fpufix.pro
@cgresizeimage.pro
@cgcheckforsymbols.pro
@cgplot.pro
@cgbitget.pro
@convert_to_type.pro
@colorsareidentical.pro
@cntrd.pro
@tvcircle.pro
@cgplots.pro
@sixlin.pro
@cv_coord.pro
@where_xyz.pro
@gcntrd.pro
@polywarp.pro

pro IRDIS_distortion

;assumptions:
;holes in pinhole mask are evenly spaced
;field distortions are smallest in the center of the image

;=============================================================================================

readcol, '/home/amueller/work/IDLlibs/AO/SPHERE/IRDIS_mask.txt', dx, dy, format='d,d', /silent

;=============================================================================================

readcol, '/home/amueller/work/IDLlibs/AO/SPHERE/datapaths.txt', tmp, format='a', /silent
print, ''
for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]

print, ''
read, 'Select Path: ', selp
path = tmp[selp-1]+'/RAW/IRDIS/'

;=============================================================================================

dim = 800

;=============================================================================================

reffile = '/home/amueller/work/IDLlibs/AO/SPHERE/irdis_distortion_points.dat'

readcol, reffile, xr, yr, format='d,d', /silent

diff = ts_diff(xr,1)
idx = where(diff lt 0.)
distance = abs(median(diff[idx]))	;reference distance between the points

;=============================================================================================

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

;selection of distortion frame

idx = where(catg eq 'CALIB' and type eq 'LAMP,DISTORT')
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
read, 'Distortion Frame: ', sel
dgf = tfile[sel-1]

;selection of CALIBRATION frames

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
read, 'Bad Pixel Mask: ', sel
bpf = tfile[sel-1]
read, 'Master Flat: ', sel
fff = tfile[sel-1]
read, 'Master Dark: ', sel
bgf = tfile[sel-1]


;=============================================================================================

;read in files

bg = readfits(bgf, hdrbg, exten_no=0, /silent)
ff = mrdfits(fff, 0, hdrff, /silent)
bp = mrdfits(bpf, 0, hdrbg, /silent)
dg = mrdfits(dgf, 0, hdrdg, /silent)

bgl = bg[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]
ffl = ff[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]
bpl = bp[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]
dgl = dg[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]

bgr = bg[1024:2047,*]
ffr = ff[1024:2047,*]
bpr = bp[1024:2047,*]
dgr = dg[1024:2047,*]

bgr = bgr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]
ffr = ffr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]
bpr = bpr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]
dgr = dgr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]

;=============================================================================================

;image cosmetics

rawbgffl = (dgl-bgl)/ffl
fixpix_mod, rawbgffl, bpl, outim, npix=24, /weight, /silent
dgl = outim

rawbgffr = (dgr-bgr)/ffr
fixpix_mod, rawbgffr, bpr, outim, npix=24, /weight, /silent
dgr = outim

;=============================================================================================

;left detector

;---------------------------------------------------------------------------------------------

;identify central spots

cgloadct, 3
window, 0, xs=dim+200, ys=dim+200
device, cursor_standard=2
  cgimage, dgl, /axis, stretch=2, title='Left detector', axkeywords={ticklen:-0.02}; minvalue=stretch[0], maxvalue=stretch[1]
  oplot, !x.crange, dim/2.*[1,1], linestyle=2
  oplot, dim/2.*[1,1], !y.crange, linestyle=2

!mouse.button = 0
print, ''
print, 'Select a central spot'

cursor, x, y, 3, /data

if (!mouse.button eq 1) then begin

  !mouse.button = 0
  cntrd, dgl, x, y, xc0, yc0, 3.5, extendbox=2

endif

;find the 4 surrounding spots

cntrd, dgl, xc0-distance, yc0, xc1, yc1, 3.5
cntrd, dgl, xc0+distance, yc0, xc2, yc2, 3.5
cntrd, dgl, xc0, yc0-distance, xc3, yc3, 3.5
cntrd, dgl, xc0, yc0+distance, xc4, yc4, 3.5

cntrd, dgl, xc0-distance, yc0+distance, xc5, yc5, 3.5
cntrd, dgl, xc0+distance, yc0+distance, xc6, yc6, 3.5
cntrd, dgl, xc0-distance, yc0-distance, xc7, yc7, 3.5
cntrd, dgl, xc0+distance, yc0-distance, xc8, yc8, 3.5


tvcircle, /data, 5, xc0, yc0, color=cgcolor('red')
tvcircle, /data, 5, xc1, yc1, color=cgcolor('green')
tvcircle, /data, 5, xc2, yc2, color=cgcolor('green')
tvcircle, /data, 5, xc3, yc3, color=cgcolor('green')
tvcircle, /data, 5, xc4, yc4, color=cgcolor('green')
tvcircle, /data, 5, xc5, yc5, color=cgcolor('green')
tvcircle, /data, 5, xc6, yc6, color=cgcolor('green')
tvcircle, /data, 5, xc7, yc7, color=cgcolor('green')
tvcircle, /data, 5, xc8, yc8, color=cgcolor('green')

;---------------------------------------------------------------------------------------------

;get abritary rotation by fitting the spots

; sixlin, [xc1,xc0,xc2], [yc1,yc0,yc2], ah, sigah, bh, sigbh	;horizontal
; sixlin, [xc3,xc0,xc4], [yc3,yc0,yc4], av, sigav, bv, sigbv	;vertical
; sixlin, [xc1,xc2], [yc1,yc2], ah, sigah, bh, sigbh	;horizontal
; sixlin, [xc5,xc6], [yc5,yc6], ah2, sigah2, bh2, sigbh2	;horizontal
; sixlin, [xc7,xc8], [yc7,yc8], ah3, sigah3, bh3, sigbh3	;horizontal

sixlin, [xc5,xc4], [yc5,yc4], ah, sigah, bh, sigbh	;horizontal
sixlin, [xc4,xc6], [yc4,yc6], ah2, sigah2, bh2, sigbh2	;horizontal
sixlin, [xc1,xc0], [yc1,yc0], ah3, sigah3, bh3, sigbh3	;horizontal
sixlin, [xc0,xc2], [yc0,yc2], ah4, sigah4, bh4, sigbh4	;horizontal
sixlin, [xc7,xc3], [yc7,yc3], ah5, sigah5, bh5, sigbh5	;horizontal
sixlin, [xc3,xc8], [yc3,yc8], ah6, sigah6, bh6, sigbh6	;horizontal

; sixlin, [xc3,xc4], [yc3,yc4], av, sigav, bv, sigbv	;vertical
; sixlin, [xc5,xc7], [yc5,yc7], av2, sigav2, bv2, sigbv2	;vertical
; sixlin, [xc6,xc8], [yc6,yc8], av3, sigav3, bv3, sigbv3	;vertical

angh = atan(bh[0])/!dtor
angh2 = atan(bh2[0])/!dtor
angh3 = atan(bh3[0])/!dtor
angh4 = atan(bh4[0])/!dtor
angh5 = atan(bh5[0])/!dtor
angh6 = atan(bh6[0])/!dtor

ang = median([angh,angh2,angh3,angh4,angh5,angh6])
print, ang
; angh = atan(bh[0])/!dtor
; ; angv = atan(bv[0])/!dtor;-90.
; angh2 = atan(bh2[0])/!dtor
; ; angv2 = atan(bv2[0])/!dtor;-90.
; angh3 = atan(bh3[0])/!dtor
; ; angv3 = atan(bv3[0])/!dtor;-90.
; ; ang = mean([angh,angv])
; ang = mean([angh,angh2,angh3])

window, 4
plot, [angh, angh2, angh3, angh4, angh5, angh6], xr=[-1,6], xst=1, psym=2, symsize=2, /yn
oplot, !x.crange, ang*[1,1]

;=============================================================================================

;generate a synthetic grid as estimate positions for cntrd and as reference to be compared with the measurements

maxspots = floor(dim/distance)
xguess = dblarr(maxspots)
yguess = dblarr(maxspots)

x = (dindgen(maxspots)*distance+1.) # (dblarr(maxspots) + 1)
y = (dblarr(maxspots) + 1) # (dindgen(maxspots)*distance+1.)

tx = x[*]
ty = y[*]

polar = dblarr(2,n_elements(tx))
for i=0,n_elements(tx)-1 do polar[*,i] = cv_coord(/degrees, /double, from_rect=[tx[i]-dim/2.,ty[i]-dim/2.], /to_polar)

polar[0,*] = polar[0,*]+ang

rect = dblarr(2,n_elements(tx))
for i=0,n_elements(tx)-1 do rect[*,i] = cv_coord(/degrees, /double, from_polar=polar[*,i], /to_rect)

x = reform(rect[0,*])+dim/2 & y = reform(rect[1,*])+dim/2

if (min(x) lt 0.) then x = x+distance/2.
if (min(y) lt 0.) then y = y+distance/2.

;=============================================================================================

;bring synthetic grid over xc0, yc0

tx = x[*]
ty = y[*]
td = sqrt((tx-dim/2.)^2.+(ty-dim/2.)^2.)

;most likely point of the synthetic grid which we now have to shift over xc0, yc0
dum = min(td, idx)
dx = x[idx]-xc0
dy = y[idx]-yc0

xref = x-dx
yref = y-dy

;cut out regions with clusters of bad pixels

ind = where_xyz(xref lt 305. or xref gt 334. or yref gt 241.5 or yref lt 213., XIND=xind);, YIND=yind, ZIND=zind)
nrej1 = n_elements(xref)-n_elements(ind)
; print, 'Rejected '+strcompress(n_elements(xref)-n_elements(ind),/rem)+' points'
xref = xref[ind]
yref = yref[ind]

;do not consider the edges because of bad pixel map (usually vignetted)
;reject every point from with a central distance greater than 495px
d = sqrt((xref-dim/2.)^2.+(yref-dim/2.)^2.)
idxd = where(d lt 500.)
if (idxd[0] ne -1) then begin

  nrej2 = n_elements(xref)-n_elements(idxd)
;   print, 'Rejected '+strcompress(n_elements(xref)-n_elements(idxd),/rem)+' points'

  xref = xref[idxd]
  yref = yref[idxd]

endif

oplot, xref, yref, psym=sym(1), symsize=0.5, color=cgcolor('green')

;=============================================================================================

;measure the postiions of all spots

mx = dblarr(n_elements(xref)) & my = mx

for i=0,n_elements(xref)-1 do begin

  cntrd, dgl, xref[i]+(0.01*xref[i]-(0.01*xc0/2.)), yref[i], t1, t2, 3.5, extendbox=2, /silent
  ;gcntrd, dgl, xref[i]+(0.01*xref[i]-(0.01*xc0/2.)), yref[i], t1, t2, 3.5, /silent
  mx[i] = t1
  my[i] = t2

endfor

;do not consider failed detections
idxx = where(mx ne -1.)
ntotx = n_elements(xref)
if (idxx[0] ne -1) then begin

  mx = mx[idxx]
  my = my[idxx]
  xref = xref[idxx]
  yref = yref[idxx]

endif

idxy = where(my ne -1.)
ntoty = n_elements(yref)
if (idxy[0] ne -1) then begin

  mx = mx[idxy]
  my = my[idxy]
  xref = xref[idxy]
  yref = yref[idxy]

endif

nrej3 = (ntotx-n_elements(idxx))+(ntoty-n_elements(idxy))
;print, 'Rejected '+strcompress((ntotx-n_elements(idxx))+(ntoty-n_elements(idxy)),/rem)+' points'

print, ''
print, 'Rejected '+strcompress(nrej1+nrej2+nrej3,/rem)+' points'

oplot, mx, my, psym=sym(1), symsize=1.5, color=cgcolor('blue')

;=============================================================================================

;compute distortion

polywarp, mx, my, xref, yref, 3, kx, ky, status=st
corrim = poly_2d(dgl, kx, ky, 2, cubic=-0.5)

cgloadct, 3
window, 2, xs=dim+200, ys=dim+200
cgimage, corrim, /axis, stretch=2, title='Left detector'; minvalue=stretch[0], maxvalue=stretch[1]
oplot, xref, yref, psym=sym(1), symsize=0.5, color=cgcolor('green')
oplot, !x.crange, dim/2.*[1,1], linestyle=2
oplot, dim/2.*[1,1], !y.crange, linestyle=2

coeffl = {kx:kx, ky:ky}
save, coeffl, filename=path+'distortion_coeffs_left.sav'

;=============================================================================================
;=============================================================================================


;right detector

;---------------------------------------------------------------------------------------------

;identify central spots

cgloadct, 3
window, 1, xs=dim+200, ys=dim+200
device, cursor_standard=2
  cgimage, dgr, /axis, stretch=2, title='Right detector', axkeywords={ticklen:-0.02}; minvalue=stretch[0], maxvalue=stretch[1]
  oplot, !x.crange, dim/2.*[1,1], linestyle=2
  oplot, dim/2.*[1,1], !y.crange, linestyle=2

!mouse.button = 0
print, ''
print, 'Select a central spot'

cursor, x, y, 3, /data

if (!mouse.button eq 1) then begin

  !mouse.button = 0
  cntrd, dgr, x, y, xc0, yc0, 3.5, extendbox=2

endif

;find the 4 surrounding spots

cntrd, dgr, xc0-distance, yc0, xc1, yc1, 3.5;, extendbox=2
cntrd, dgr, xc0+distance, yc0, xc2, yc2, 3.5;, extendbox=2
cntrd, dgr, xc0, yc0-distance, xc3, yc3, 3.5;, extendbox=2
cntrd, dgr, xc0, yc0+distance, xc4, yc4, 3.5;, extendbox=2

cntrd, dgr, xc0-distance, yc0+distance, xc5, yc5, 3.5;, extendbox=2
cntrd, dgr, xc0+distance, yc0+distance, xc6, yc6, 3.5;, extendbox=2
cntrd, dgr, xc0-distance, yc0-distance, xc7, yc7, 3.5;, extendbox=2
cntrd, dgr, xc0+distance, yc0-distance, xc8, yc8, 3.5;, extendbox=2


tvcircle, /data, 5, xc0, yc0, color=cgcolor('red')
tvcircle, /data, 5, xc1, yc1, color=cgcolor('green')
tvcircle, /data, 5, xc2, yc2, color=cgcolor('green')
tvcircle, /data, 5, xc3, yc3, color=cgcolor('green')
tvcircle, /data, 5, xc4, yc4, color=cgcolor('green')
tvcircle, /data, 5, xc5, yc5, color=cgcolor('green')
tvcircle, /data, 5, xc6, yc6, color=cgcolor('green')
tvcircle, /data, 5, xc7, yc7, color=cgcolor('green')
tvcircle, /data, 5, xc8, yc8, color=cgcolor('green')

;---------------------------------------------------------------------------------------------

;get abritary rotation by fitting the spots

; sixlin, [xc1,xc0,xc2], [yc1,yc0,yc2], ah, sigah, bh, sigbh	;horizontal
; sixlin, [xc3,xc0,xc4], [yc3,yc0,yc4], av, sigav, bv, sigbv	;vertical
; sixlin, [xc1,xc2], [yc1,yc2], ah, sigah, bh, sigbh	;horizontal
; sixlin, [xc5,xc6], [yc5,yc6], ah2, sigah2, bh2, sigbh2	;horizontal
; sixlin, [xc7,xc8], [yc7,yc8], ah3, sigah3, bh3, sigbh3	;horizontal

sixlin, [xc5,xc4], [yc5,yc4], ah, sigah, bh, sigbh	;horizontal
sixlin, [xc4,xc6], [yc4,yc6], ah2, sigah2, bh2, sigbh2	;horizontal
sixlin, [xc1,xc0], [yc1,yc0], ah3, sigah3, bh3, sigbh3	;horizontal
sixlin, [xc0,xc2], [yc0,yc2], ah4, sigah4, bh4, sigbh4	;horizontal
sixlin, [xc7,xc3], [yc7,yc3], ah5, sigah5, bh5, sigbh5	;horizontal
sixlin, [xc3,xc8], [yc3,yc8], ah6, sigah6, bh6, sigbh6	;horizontal

; sixlin, [xc3,xc4], [yc3,yc4], av, sigav, bv, sigbv	;vertical
; sixlin, [xc5,xc7], [yc5,yc7], av2, sigav2, bv2, sigbv2	;vertical
; sixlin, [xc6,xc8], [yc6,yc8], av3, sigav3, bv3, sigbv3	;vertical

angh = atan(bh[0])/!dtor
angh2 = atan(bh2[0])/!dtor
angh3 = atan(bh3[0])/!dtor
angh4 = atan(bh4[0])/!dtor
angh5 = atan(bh5[0])/!dtor
angh6 = atan(bh6[0])/!dtor

ang = median([angh,angh2,angh3,angh4,angh5,angh6])
print, ang
; angh = atan(bh[0])/!dtor
; ; angv = atan(bv[0])/!dtor;-90.
; angh2 = atan(bh2[0])/!dtor
; ; angv2 = atan(bv2[0])/!dtor;-90.
; angh3 = atan(bh3[0])/!dtor
; ; angv3 = atan(bv3[0])/!dtor;-90.
; ; ang = mean([angh,angv])
; ang = mean([angh,angh2,angh3])

window, 5
plot, [angh, angh2, angh3, angh4, angh5, angh6], xr=[-1,6], xst=1, psym=2, symsize=2, /yn
oplot, !x.crange, ang*[1,1]

;=============================================================================================

;generate a synthetic grid as estimate positions for cntrd and as reference to be compared with the measurements

maxspots = floor(dim/distance)
xguess = dblarr(maxspots)
yguess = dblarr(maxspots)

x = (dindgen(maxspots)*distance+1.) # (dblarr(maxspots) + 1)
y = (dblarr(maxspots) + 1) # (dindgen(maxspots)*distance+1.)

tx = x[*]
ty = y[*]

polar = dblarr(2,n_elements(tx))
for i=0,n_elements(tx)-1 do polar[*,i] = cv_coord(/degrees, /double, from_rect=[tx[i]-dim/2.,ty[i]-dim/2.], /to_polar)

polar[0,*] = polar[0,*]+ang

rect = dblarr(2,n_elements(tx))
for i=0,n_elements(tx)-1 do rect[*,i] = cv_coord(/degrees, /double, from_polar=polar[*,i], /to_rect)

x = reform(rect[0,*])+dim/2 & y = reform(rect[1,*])+dim/2

if (min(x) lt 0.) then x = x+distance/2.
if (min(y) lt 0.) then y = y+distance/2.

;=============================================================================================

;bring synthetic grid over xc0, yc0

tx = x[*]
ty = y[*]
td = sqrt((tx-dim/2.)^2.+(ty-dim/2.)^2.)

;most likely point of the synthetic grid which we now have to shift over xc0, yc0
dum = min(td, idx)
dx = x[idx]-xc0
dy = y[idx]-yc0

xref = x-dx
yref = y-dy

;cut out regions with clusters of bad pixels

ind = where_xyz(xref lt 228.7 or xref gt 243.8 or yref gt 195.1 or yref lt 182.1, XIND=xind);, YIND=yind, ZIND=zind)
nrej1 = n_elements(xref)-n_elements(ind)
; print, 'Rejected '+strcompress(n_elements(xref)-n_elements(ind),/rem)+' points'
xref = xref[ind]
yref = yref[ind]

ind = where_xyz(xref lt 310.8 or xref gt 328.2 or yref gt 242. or yref lt 215., XIND=xind);, YIND=yind, ZIND=zind)
nrej1a = n_elements(xref)-n_elements(ind)
; print, 'Rejected '+strcompress(n_elements(xref)-n_elements(ind),/rem)+' points'
xref = xref[ind]
yref = yref[ind]

;do not consider the edges because of bad pixel map (usually vignetted)
;reject every point from with a central distance greater than 495px
d = sqrt((xref-dim/2.)^2.+(yref-dim/2.)^2.)
idxd = where(d lt 500.)
if (idxd[0] ne -1) then begin

  nrej2 = n_elements(xref)-n_elements(idxd)
;   print, 'Rejected '+strcompress(n_elements(xref)-n_elements(idxd),/rem)+' points'

  xref = xref[idxd]
  yref = yref[idxd]

endif

oplot, xref, yref, psym=sym(1), symsize=0.5, color=cgcolor('green')

;=============================================================================================

;measure the postiions of all spots

mx = dblarr(n_elements(xref)) & my = mx

for i=0,n_elements(xref)-1 do begin

  cntrd, dgr, xref[i]+(0.01*xref[i]-(0.01*xc0/2.)), yref[i], t1, t2, 3.5, extendbox=2, /silent
  ;gcntrd, dgr, xref[i]+(0.01*xref[i]-(0.01*xc0/2.)), yref[i], t1, t2, 3.5, /silent
  mx[i] = t1
  my[i] = t2

endfor

;do not consider failed detections
idxx = where(mx ne -1.)
ntotx = n_elements(xref)
if (idxx[0] ne -1) then begin

  mx = mx[idxx]
  my = my[idxx]
  xref = xref[idxx]
  yref = yref[idxx]

endif

idxy = where(my ne -1.)
ntoty = n_elements(yref)
if (idxy[0] ne -1) then begin

  mx = mx[idxy]
  my = my[idxy]
  xref = xref[idxy]
  yref = yref[idxy]

endif

nrej3 = (ntotx-n_elements(idxx))+(ntoty-n_elements(idxy))
;print, 'Rejected '+strcompress((ntotx-n_elements(idxx))+(ntoty-n_elements(idxy)),/rem)+' points'

print, ''
print, 'Rejected '+strcompress(nrej1+nrej1a+nrej2+nrej3,/rem)+' points'

oplot, mx, my, psym=sym(1), symsize=1.5, color=cgcolor('blue')

;=============================================================================================

;compute distortion

polywarp, mx, my, xref, yref, 3, kx, ky, status=st
corrim = poly_2d(dgr, kx, ky, 2, cubic=-0.5)

cgloadct, 3
window, 3, xs=dim+200, ys=dim+200
cgimage, corrim, /axis, stretch=2, title='Right detector'; minvalue=stretch[0], maxvalue=stretch[1]
oplot, xref, yref, psym=sym(1), symsize=0.5, color=cgcolor('green')
oplot, !x.crange, dim/2.*[1,1], linestyle=2
oplot, dim/2.*[1,1], !y.crange, linestyle=2

coeffr = {kx:kx, ky:ky}
save, coeffr, filename=path+'distortion_coeffs_right.sav'


stop
end
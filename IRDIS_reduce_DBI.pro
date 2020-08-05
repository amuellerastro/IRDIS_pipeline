@array_indices.pro
@fixpix_mod.pro
@strc.pro
@ten.pro
@arrdelete.pro
@caldat.pro
@calcpos.pro
@addpm.pro
@adstring.pro
@sixty.pro
@imgclean.pro
@strn.pro
@extrac.pro
@skyline.pro
@stdev.pro
@avg.pro
@starchck.pro
@vect.pro
; @fshift.pro
@shift_sub.pro
@fftshift.pro
@eq2hor_MOD.pro
@radec.pro
@precess.pro
@premat.pro
@co_nutate.pro
@nutate.pro
@poly.pro
@cirrange.pro
@isarray.pro
@co_aberration.pro
@sunpos.pro
@ct2lst.pro
@hadec2altaz.pro
@co_refract.pro
@parangle.pro
@mwrfits.pro
@fxaddpar.pro
@mrdfits.pro
@fxposit.pro
@mrd_hread.pro
@fxpar.pro
@gettok.pro
@valid_num.pro
@mrd_skip.pro
@fxmove.pro
@match.pro
@mrd_struct.pro
@get_eso_keyword.pro
@strsplit.pro
@setdefaultvalue.pro
@dist_circle.pro
@detabify.pro
@fxparpos.pro
@headfits.pro
@fxhmake.pro
@fxwrite.pro
@check_fits.pro
@host_to_ieee.pro
@fxhclean.pro
@sxdelpar.pro
@la_cosmic_MOD.pro
@fxread.pro
@fxhread.pro
@ieee_to_host.pro
@get_date.pro
@daycnv.pro
@reverse.pro
@readcol.pro
@remchar.pro
@strnumber.pro
@linspace.pro
@sixlin.pro
@circint_MOD.pro
@mean.pro
@moment.pro
@mpfit2dpeak.pro
@mpfit.pro
@mpfit2dfun.pro
@readfits.pro
@sxpar.pro
@writefits.pro
@mkhdr.pro
@sxaddpar.pro
@fits_add_checksum.pro
@checksum32.pro
@n_bytes.pro
@is_ieee_big.pro
@fits_ascii_encode.pro


; function func_surface, x, y, p
; 
; ;   c = poly(x,p)
; ;   c = poly(y,c)
; ; 
; ;   return, c
; 
;   end
; 
; function func_surface_p1, x, y, p
; 
;   fit = p[0] + p[1]*x + p[2]*y
;   return, fit
; 
; end
; 
; function func_surface_p2, x, y, p
; 
;   fit = p[0] + p[1]*x + p[2]*y + p[3]*x*y + p[4]*x^2. + p[5]*y^3.
;   return, fit
; 
; end

; function func_surface_p3, x, y, p
; 
;   fit = p[0] + p[1]*x + p[2]*y + p[3]*x*y + p[4]*x^2. + p[5]*y^3. + p[6]*x^3. + p[7]*x^2.*y + p[8]*x*y^2. + p[9]*y^3.
;   return, fit
; 
; end


pro IRDIS_reduce_DBI

;=============================================================================================

; star = ['HD95086', 'TWHya', 'betaPic', 'V471Tau', 'GJ1046', 'Sirius', 'HD159073', 'TCha']

;select star and get its needed properties

filestar = file_search('/home/amueller/work/IDLlibs/AO/TargetProperties/Targets/*.sav', count=nstars)
stars = strarr(nstars)
for i=0,nstars-1 do begin

  stars[i] = strmid(filestar[i], 56, strlen(filestar[i])-56)
  stars[i] = strmid(stars[i], 0, strlen(stars[i])-4)

endfor

nstars = n_elements(stars)
print, ''
for i=0,nstars-1 do print, strcompress(i+1, /rem), ' ', stars[i]
selstar = ''
read, 'Select Star?: ', selstar
star = stars[selstar-1]
;selstar = '1'
restore, filestar[selstar-1], /verbose

ra_st = ra
dec_st = dec


;=============================================================================================

readcol, '/home/amueller/work/IDLlibs/AO/SPHERE/IRDIS_mask.txt', dx, dy, format='d,d', /silent

;=============================================================================================

readcol, '/home/amueller/work/IDLlibs/AO/SPHERE/datapaths.txt', tmp, format='a', /silent
print, ''

idx = where(strmatch(tmp, '*'+star+'*') eq 1)
tmp = tmp[idx]

for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]

print, ''
read, 'Select Path: ', selp
path = tmp[selp-1]+'/RAW/IRDIS/'

if (strmatch(path, '*'+star+'*') ne 1) then begin

  print, ''
  print, 'Wrong star selected! Stop.'
  stop

endif

;=============================================================================================

;IRDIS always has 2048x1024px for both chips, 1 chip is 1024x1024
; xc = 512.
; yc = 512.

;cut images down to 800x800px
;using bad pixel map of FF too and all the edges are full with bad pixels -> takes a very long time for fixpix_mod
; xc = 512
; yc = 512
dim = 800

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

;selection of SCIENCE frames

idx = where(catg eq 'SCIENCE' and type eq 'OBJECT')
if (n_elements(idx) lt 2) then idx = where(catg eq 'SCIENCE' and type eq 'OBJECT,CENTER')
; idx = where(catg eq 'CALIB' and type eq 'OBJECT,ASTROMETRY')
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
read, 'Science Frames (first and last number): ', sel
rawf = tfile[indgen(sel[1]-sel[0]+1)+sel[0]-1]

nfiles = n_elements(rawf)

sel = 0

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
; read, 'Bad Pixel Mask: ', sel
; bpf = tfile[sel-1]
bpf = path+'static_badpixels.fits'
; read, 'Master Flat: ', sel
; fff = tfile[sel-1]
fff = path+'instrument_flat.fits'
; read, 'Star Center (0 if none): ', sel
; if (sel gt 0) then begin
;   scf = tfile[sel-1]
;   scflag = 'y'
; endif else begin
;   scflag = 'n'
; endelse
match = strmatch(tfile, '*star_center*')
idx = where(match eq 1)
if (idx[0] ne -1) then begin
  scf = tfile[idx]
  scflag = 'y'
endif else begin
  scflag = 'n'
endelse
; read, 'Sky BG (0 if none): ', sel
; if (sel gt 0) then begin
;   skf = tfile[sel-1]
;   skyflag = 'y'
;   darkflag = 'n'
; endif else begin
;   skyflag = 'n'
; endelse
match = strmatch(tfile, '*sky_background*')
idx = where(match eq 1)
if (idx[0] ne -1) then begin
  skf = tfile[idx]
  skyflag = 'y'
  darkflag = 'n'
endif else begin
  skyflag = 'n'
endelse
if (skyflag eq 'n') then begin
  read, 'Master Dark: ', sel
  if (sel gt 0) then begin
    bgf = tfile[sel-1]
    darkflag = 'y'
  endif else begin
    darkflag = 'n'
  endelse
endif
distquest = ''
; read, 'With distortion correction? y/n: ', distquest
read, 'Anamorphic Correction? y/n: ', distquest

;=============================================================================================

;read in files

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

if (darkflag eq 'y') then bg = readfits(bgf, hdrbg, exten_no=0, /silent)
ff = mrdfits(fff, 0, hdrff, /silent)
bp = mrdfits(bpf, 0, hdrbg, /silent)
if (skyflag eq 'y') then sk = mrdfits(skf, 0, hdrsk, /silent)
if (scflag eq 'y') then sc = mrdfits(scf, 1, hdrsc, /silent)



; ;correct for anamorphism
; does not work here because BP map gets destroyed
; if (distquest eq 'y') then begin
; 
;   fac = 1.006
; 
;   if (darkflag eq 'y') then begin
; 
;     bg = congrid(bg, 2048, 1024.*fac, /center, cubic=-0.5)
;     bg = bg[*, 3:1026]
; 
;   endif
; 
;   ff = congrid(ff, 2048, 1024.*fac, /center, cubic=-0.5)
;   ff = ff[*, 3:1026]
;   bp = congrid(bp, 2048, 1024.*fac, /center, cubic=-0.5)
;   bp = bp[*, 3:1026]
; 
;   if (skyflag eq 'y') then begin
; 
;     sk = congrid(sk, 2048, 1024.*fac, /center, cubic=-0.5)
;     sk = sk[*, 3:1026]
; 
;   endif
; 
; endif


if (darkflag eq 'y') then bgl = bg[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]
ffl = ff[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]
bpl = bp[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]
if (skyflag eq 'y') then skl = sk[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]

if (darkflag eq 'y') then bgr = bg[1024:2047,*]
ffr = ff[1024:2047,*]
bpr = bp[1024:2047,*]
if (skyflag eq 'y') then skr = sk[1024:2047,*]

if (darkflag eq 'y') then bgr = bgr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]
ffr = ffr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]
bpr = bpr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]
if (skyflag eq 'y') then skr = skr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]

;=============================================================================================

;values based on http://wiki.oamp.fr/sphere/AstrometricCalibration
;PUPOFFSET = 135.87d0;+/-0.03 deg 
PUPOFFSET = 135.99d0;+/-0.11 deg 
;True_North= -1.764d0;+/-0.010 deg
True_North= 0.0d0 ;dummy, doinf TN in Astrometry_Photometry

;=============================================================================================

;Observatory parameters

;from http://www.eso.org/sci/facilities/paranal/astroclimate/site.html
lat = ten(-24.d0, 37.d0, 30.300d0)	;for UT3
lon = ten(-70.d0, 24.d0, 9.896d0)	;for UT3
altitude = 2635.43d0

;=============================================================================================

dum = mrdfits(rawf[0], 0, hdrraw, /silent)
sz = size(dum)
if (n_elements(sz) eq 5) then sz3 = 1 else sz3 = sz[3]
parang = dblarr(nfiles, sz3)
hang = dblarr(nfiles, sz3)
PA_onsky = parang
newjd = parang
date = dblarr(nfiles, sz3, 6)
left = dblarr(dim,dim,sz3)
right = dblarr(dim,dim,sz3)

;---------------------------------------------------------------------------------------------

for i=0,nfiles-1 do begin

  print, ''
  print, 'Reducing file '+strcompress(i+1,/rem)+' / '+strcompress(nfiles,/rem)
  print, ''

  ;-------------------------------------------------------------------------------------------

  im = mrdfits(rawf[i], 0, hdrraw, /silent)

  ;output name
  pos1 = strpos(rawf[i], '/', /reverse_search)
  pos2 = strpos(rawf[i], '.fits', /reverse_search)
  out = strmid(rawf[i], pos1+1, pos2-pos1-1)
  exptime = double(get_eso_keyword(hdrraw, 'HIERARCH ESO DET SEQ1 EXPTIME'))
  ndit = double(get_eso_keyword(hdrraw, 'HIERARCH ESO DET NDIT'))
  dit = double(get_eso_keyword(hdrraw,'HIERARCH ESO DET SEQ1 DIT'))
  jd = double(get_eso_keyword(hdrraw, 'MJD-OBS'))+2400000.5d0

  ;jitter
  pacx = double(get_eso_keyword(hdrraw, 'HIERARCH ESO INS1 PAC X'))/18.d0	;in px
  pacy = double(get_eso_keyword(hdrraw, 'HIERARCH ESO INS1 PAC Y'))/18.d0	;in px

  pres = double(get_eso_keyword(hdrraw, 'HIERARCH ESO TEL AMBI PRES START'))
  temp = double(get_eso_keyword(hdrraw, 'HIERARCH ESO TEL AMBI TEMP'))+273.15d0

  ra = [double(strmid(ra_st,0,2)),double(strmid(ra_st,2,2)),double(strmid(ra_st,4,6))]
  sign = strmid(dec_st,0,1)
  if (sign eq '-') then sign = -1.d0 else sign = 1.d0
  dec = [sign*double(strmid(dec_st,1,2)), double(strmid(dec_st,3,2)), double(strmid(dec_st,5,6))]
  ra = ten(ra)*15.
  dec = ten(dec)

  epoch = double(get_eso_keyword(hdrraw, 'HIERARCH ESO TEL TARG EPOCH'))


  ;have to remove naxis keyword, othwerwise mrdfits complains
;   for j=0,n_elements(hdrraw)-1 do begin
; 
;     if (strmatch(hdrraw[j], '*NAXIS3*') eq 1) then newhdr = arrdelete(hdr, at=j, length=1)
; 
;   endfor

  ;-------------------------------------------------------------------------------------------

  for j=0,sz3-1 do begin
;   for j=4,4 do begin

    print, ''
    print, 'Reducing frame '+strcompress(j+1,/rem)+' / '+strcompress(sz3,/rem)
    print, ''

    newjd[i,j] = jd + ((exptime/ndit)*j + (exptime/ndit)/2.)/86400.d0
    caldat, newjd[i,j], mo, day, yr, hh, mm, ss

    ;get date of the form yyyy-mm-ddThh:mm:ss
    date[i,j,0] = yr
    date[i,j,1] = mo
    date[i,j,2] = day
    date[i,j,3] = hh
    date[i,j,4] = mm
    date[i,j,5] = ss

    ;get year of observations for new RA, DEC, i.e. epoch
    hour = hh+mm/60.d0+ss/3600.d0
    if ((double(yr) mod 4.) eq 0.) then date2 = yr+mo/12.d0+day/366.d0+hour/8766.0d0 $
	else date2 = yr+mo/12.d0+day/365.d0+hour/8766.0d0	;takes leap year into account

;     calcpos, ra, dec, pma, pmd, radvel, plx, epoch, date2, epoch, date2, epoch, 'Y', 'Y', outra, outdec

    if (finite(radvel[0] ne 1)) then radvel = 0.
    ;if (finite(plx[0] ne 1)) then plx = 1./100.
    if (finite(plx) ne 1) then begin
      ;plx = 1./100.
      read, 'prallax in arcsec: ', plx
      ;read, 'error in arcsec: ', eplx
    endif

    ;compute current cordinates and parallactic angle at time of observation
    dt = date2-epoch

    eq2hor_MOD, ra, dec, pma, pmd, radvel, plx, dt, newjd[i,j], alt, az, ha, lat=lat, lon=lon, altitude=altitude, pres=pres, temp=temp, outra=outra, outdec=outdec;, /verbose

    hang[i,j] = ha/15.
    parang[i,j] = parangle(ha/15., outdec, lat)
;     print, parang

    if (hang[i,j] gt 12.) then hang[i,j] = hang[i,j]-24.
    if (hang[i,j] lt -12.) then hang[i,j] = hang[i,j]+24.


    if (outdec gt lat and parang[i,j] lt 0.) then parang[i,j] = parang[i,j]+360.d0

    PA_onsky[i,j] = parang[i,j] + True_North + PUPOFFSET

    ;-------------------------------------------------------------------------------------------

    ;split image
    iml = im[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1,j]
    imr = im[1024:2047,*,j]
    imr = imr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]

    ;-------------------------------------------------------------------------------------------

;do image cosmetics

    if (dit gt 10.) then readn = 10. else readn = 4.	;DIT of 10s chosen abritrarily

;left detector

    rawl = iml
    rawbgl = rawl;-bgl
    if (skyflag eq 'y') then rawbgffl = (rawbgl-skl)/ffl else rawbgffl = (rawbgl-bgl)/ffl
    fixpix_mod, rawbgffl, bpl, outim, npix=24, /weight;, /silent
    ;outim = maskinterp(rawbgff, (-1.)*(bpl-1.), 1, 6, )
    corrl = outim

    ;check for 'cosmic rays' or left overs
    ;imgclean, outim, star_psf_sens=0.25, sky_val_samp_sz=16, hi_scan_height=5, /fine_clean, phase1_iter=4
;     if (skyflag eq 'y') then la_cosmic_MOD, corrl, outim, path, gain=1.75, readn=readn, skyval=median(skl,/even) else la_cosmic_MOD, corrl, outim, path, gain=1.75, readn=readn
;     crl = outim
    crl = corrl

    idxnan = where(finite(crl) ne 1)
    if (idxnan[0] ne -1) then begin

      print, 'NaNs present. Skipping CR correction for this frame.'
      print, ''
      crl = corrl

    endif

;right detector

    rawr = imr
    rawbgr = rawr;-bgr
    if (skyflag eq 'y') then rawbgffr = (rawbgr-skr)/ffr else rawbgffr = (rawbgr-bgr)/ffr
    fixpix_mod, rawbgffr, bpr, outim, npix=24, /weight;, /silent
    corrr = outim

    ;check for 'cosmic rays' or left overs
    ;imgclean, outim, star_psf_sens=0.25, sky_val_samp_sz=16, hi_scan_height=5, /fine_clean, phase1_iter=4
;     if (skyflag eq 'y') then la_cosmic_MOD, corrr, outim, path, gain=1.75, readn=readn, skyval=median(skr,/even) else la_cosmic_MOD, corrr, outim, path, gain=1.75, readn=readn
;     crr = outim
    crr = corrr

    idxnan = where(finite(crr) ne 1)
    if (idxnan[0] ne -1) then begin

      print, 'NaNs present. Skipping CR correction for this frame.'
      print, ''
      crr = corrr

    endif


    imleft = crl
    imright = crr

    ;-------------------------------------------------------------------------------------------

    ;subtract sky

;     imleft = imleft - skl
;     imright = imright - skr

    ;-------------------------------------------------------------------------------------------

;     ;surface fitting
; 
;     tmp1 = imleft
;     tmp2 = tmp1
;     weights = dblarr(dim, dim)
;     dumerr = weights
;     dumerr[*,*] = 1.
; ;     xa = dindgen(dim) & ya = xa
;     xa = (dindgen(dim)) # (dblarr(dim)+1)
;     ya = (dblarr(dim)) # (dindgen(dim)+1)
; 
; 
;     ;remove Speckles
;     robmean = robust_mean(tmp1, 1., sigma, numrej, goodind=good)
;     idxgood = array_indices(tmp2, good)
;     for k=0L,n_elements(idxgood[0,*])-1 do weights[idxgood[0,k], idxgood[1,k]] = 1.
; 
;     start_params = [median(imleft),0.,0.]	;p1
; ;     start_params = [median(imleft),0.1,0.1,0,0,0]	;p2
;     params = mpfit2dfun('func_surface', xa, ya, imleft, dumerr, start_params, weights=weights, yfit=fit)
; 
; ;     res = sfit(tmp1, 4)
; ;     x=poly2d(xa, ya, 2, 2, [27.3, 0.009, 0.4])


    ;-------------------------------------------------------------------------------------------

    ;removal of cross talk using medamp technique
    ;http://www.cfht.hawaii.edu/Instruments/Imaging/WIRCam/WIRCamCrosstalks.html

;     amppos = linspace(0,1024,33);-abs(dx-dim/2)
;     ampl = dblarr(32,1024,32) & ampr = ampl
; 
;     tmpl = dblarr(1024,1024)
;     tmpl[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1,*] = imleft
;     tmpr = dblarr(1024,1024)
;     tmpr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1,*] = imright
; 
;     for k=0,31 do begin
; 
;       ampl[*,*,k] = tmpl[amppos[k]:amppos[k+1]-1,*]
;       ampr[*,*,k] = tmpr[amppos[k]:amppos[k+1]-1,*]
; 
;     endfor
; 
;     mampl = median(ampl, /even, dimension=3)
;     mampr = median(ampr, /even, dimension=3)
; 
;     for k=0,31 do begin
; 
;       tmpl[k*32:k*32+31,*] = tmpl[k*32:k*32+31,*]-mampl
;       tmpr[k*32:k*32+31,*] = tmpr[k*32:k*32+31,*]-mampr
; 
;     endfor
; 
;     imleft = tmpl[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]
;     imright = tmpr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]

    ;-------------------------------------------------------------------------------------------

    ;removal of fixed noise pattern / odd-even effect

    ;left detector upper and lower box median values
    ;right detector upper and lower box median values
;     ll = median(imleft[*,89:139], dimension=2)
;     ul = median(imleft[*,664:714], dimension=2)
;     lr = median(imright[*,89:139], dimension=2)
;     ur = median(imright[*,664:714], dimension=2)
; 
;     for k=0,799 do begin
; 
;       imleft[k,*] = imleft[k,*]-mean([ll[k],ul[k]])
;       imright[k,*] = imright[k,*]-mean([lr[k],ur[k]])
; 
;     endfor

    ;-------------------------------------------------------------------------------------------

;     if (distquest eq 'y') then begin
; 
;       imleft = poly_2d(imleft, kxl, kyl, 2, cubic=-0.5)
;       imright = poly_2d(imright, kxr, kyr, 2, cubic=-0.5)
; 
;     endif

    ;-------------------------------------------------------------------------------------------

    if (scflag eq 'y') then begin

      ;shift images to physical center

      ;if there were more than one star center frames taken interpolate position for current JD
      if (n_elements(sc.jd) gt 1) then begin

	sixlin, sc.jd, sc.center_left_x, a, ae, b, be
	center_left_x = b[0]*newjd[i,j]+a[0]
	sixlin, sc.jd, sc.center_left_y, a, ae, b, be
	center_left_y = b[0]*newjd[i,j]+a[0]

	sixlin, sc.jd, sc.center_right_x, a, ae, b, be
	center_right_x = b[0]*newjd[i,j]+a[0]
	sixlin, sc.jd, sc.center_right_y, a, ae, b, be
	center_right_y = b[0]*newjd[i,j]+a[0]

      endif else begin

	center_left_x = sc.center_left_x
	center_left_y = sc.center_left_y

	center_right_x = sc.center_right_x
	center_right_y = sc.center_right_y

      endelse

;       simleft = fshift(imleft, dx[0]-center_left_x-pacx, dy[0]-center_left_y-pacy)
;       simright = fshift(imright, dx[1]-center_right_x-pacx, dy[1]-center_right_y-pacy)

;       simleft = shift_sub(imleft, dx[0]-center_left_x-pacx, dy[0]-center_left_y-pacy)
;       simright = shift_sub(imright, dx[1]-center_right_x-pacx, dy[1]-center_right_y-pacy)

      wframe = 200.	;width of frame
      tmp = dblarr(dim+2.*wframe, dim+2.*wframe)
      tmp[wframe:wframe+dim-1, wframe:wframe+dim-1] = imleft
      stmp = fftshift(tmp, dx[0]-center_left_x-pacx, dy[0]-center_left_y-pacy)
      simleft = stmp[(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1,(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1]
      tmp = dblarr(dim+2.*wframe, dim+2.*wframe)
      tmp[wframe:wframe+dim-1, wframe:wframe+dim-1] = imright
      stmp = fftshift(tmp, dx[1]-center_right_x-pacx, dy[1]-center_right_y-pacy)
      simright = stmp[(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1,(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1]


      ;center at dim/2,dim/2
      ;simleft = fshift(imleft, 512.d0-center_left_x-pacx, 512.d0-center_left_y-pacy)
      ;simright = fshift(imright, 512.d0-center_right_x-pacx, 512.d0-center_right_y-pacy)

      ;-------------------------------------------------------------------------------------------

      left[*,*,j] = simleft
      right[*,*,j] = simright

    endif else begin	;no star center file was recorded, i.e. we can assume no coronograph and star is visible

      radius = 15.

      ;mximl = max(imleft, idxl)
      ;idxmaxl = array_indices(imleft, idxl)
      ;cutiml = imleft[idxmaxl[0]-radius:idxmaxl[0]+radius, idxmaxl[1]-radius:idxmaxl[1]+radius]
      
      find, imleft, fdx, fdy, flux, sharp, roundness, 2*stddev(imleft), 8., [-1,1], [0.2,1.0]
      if (n_elements(flux) gt 1) then begin
      
        tmp = max(flux,idxf)
        fdx = fdx[idxf]
        fdy = fdy[idxf]
      
      endif
      idxmaxl = [fdx,fdy]
      cutiml = imleft[idxmaxl[0]-radius:idxmaxl[0]+radius, idxmaxl[1]-radius:idxmaxl[1]+radius]
      estimates = [median(cutiml), max(flux), 4., 4., radius, radius, 0., 1.]
      xa = dindgen(n_elements(cutiml[0,*]))+1.d0 & ya = xa

;       weights = cutiml
;       idx1 = where(cutiml eq 0.)	;e.g. beta Pic close to the center which is masked out
;       if (idx1[0] ne -1) then begin
; 	idx2 = array_indices(weights, idx1)
; 	for k=0,n_elements(idx2[0,*])-1 do weights[idx2[0,k],idx2[1,k]] = 1.d5
;       endif
;       sign = signum(cutiml)
;       idx1 = where(sign eq -1.)
;       if (idx1[0] ne -1) then begin
; 	idx2 = array_indices(weights, idx1)
; 	for k=0,n_elements(idx2[0,*])-1 do weights[idx2[0,k],idx2[1,k]] = 1.d5
;       endif
;       weights = 1./sqrt(weights)

      yfit = mpfit2dpeak(cutiml, A, xa, ya, /moffat, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /tilt, /quiet);, weights=weights)

      center_left_x = A[4]-radius+idxmaxl[0]-1
      center_left_y = A[5]-radius+idxmaxl[1]-1

      ;mximr = max(imright, idxr)
      ;idxmaxr = array_indices(imright, idxr)
      ;cutimr = imright[idxmaxr[0]-radius:idxmaxr[0]+radius, idxmaxr[1]-radius:idxmaxr[1]+radius]
      ;estimates = [median(cutimr), mximr, 4., 4., radius, radius, 0., 1.]
      
      find, imright, fdx, fdy, flux, sharp, roundness, 2*stddev(imright), 8., [-1,1], [0.2,1.0]
      if (n_elements(flux) gt 1) then begin
      
        tmp = max(flux,idxf)
        fdx = fdx[idxf]
        fdy = fdy[idxf]
      
      endif
      idxmaxr = [fdx,fdy]
      cutimr = imright[idxmaxr[0]-radius:idxmaxr[0]+radius, idxmaxr[1]-radius:idxmaxr[1]+radius]
      estimates = [median(cutimr), max(flux), 4., 4., radius, radius, 0., 1.]
      
      
      xa = dindgen(n_elements(cutimr[0,*]))+1.d0 & ya = xa

;       weights = cutimr
;       idx1 = where(cutimr eq 0.)	;e.g. beta Pic close to the center which is masked out
;       if (idx1[0] ne -1) then begin
; 	idx2 = array_indices(weights, idx1)
; 	for k=0,n_elements(idx2[0,*])-1 do weights[idx2[0,k],idx2[1,k]] = 1.d5
;       endif
;       sign = signum(cutimr)
;       idx1 = where(sign eq -1.)
;       if (idx1[0] ne -1) then begin
; 	idx2 = array_indices(weights, idx1)
; 	for k=0,n_elements(idx2[0,*])-1 do weights[idx2[0,k],idx2[1,k]] = 1.d5
;       endif
;       weights = 1./sqrt(weights)

      yfit = mpfit2dpeak(cutimr, A, xa, ya, /moffat, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /tilt, /quiet);, weights=weights)

      center_right_x = A[4]-radius+idxmaxr[0]-1
      center_right_y = A[5]-radius+idxmaxr[1]-1

;       simleft = fshift(imleft, dim/2.d0-center_left_x, dim/2.d0-center_left_y)
;       simright = fshift(imright, dim/2.d0-center_right_x, dim/2.d0-center_right_y)
;       simleft = shift_sub(imleft, dim/2.d0-center_left_x, dim/2.d0-center_left_y)
;       simright = shift_sub(imright, dim/2.d0-center_right_x, dim/2.d0-center_right_y)

      wframe = 200.	;width of frame
      tmp = dblarr(dim+2.*wframe, dim+2.*wframe)
      tmp[wframe:wframe+dim-1, wframe:wframe+dim-1] = imleft
      stmp = fftshift(tmp, dim/2.d0-center_left_x, dim/2.d0-center_left_y)
      simleft = stmp[(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1,(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1]
      tmp = dblarr(dim+2.*wframe, dim+2.*wframe)
      tmp[wframe:wframe+dim-1, wframe:wframe+dim-1] = imright
      stmp = fftshift(tmp, dim/2.d0-center_right_x, dim/2.d0-center_right_y)
      simright = stmp[(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1,(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1]

      left[*,*,j] = simleft
      right[*,*,j] = simright

    endelse

    ;-------------------------------------------------------------------------------------------

  endfor

  ;-------------------------------------------------------------------------------------------

  idxnanl = where(finite(left) ne 1)
  idxnanr = where(finite(right) ne 1)

  if (idxnanl[0] ne -1) then begin

    print, 'NaNs present in left image. Stop.'
    stop

  endif

  if (idxnanr[0] ne -1) then begin

    print, 'NaNs present in right image. Stop.'
    stop

  endif

  ;-------------------------------------------------------------------------------------------

  ;Anamorphic correction
  if (distquest eq 'y') then begin

    fac = 1.006d0

    for j=0,sz3-1 do begin

      tmp = congrid(reform(left[*,*,j]), dim, dim*fac, /center, cubic=-0.5)
      diff = n_elements(tmp[0,*])-dim
      if (diff mod 2. ne 0.) then stop
      left[*,*,j] = tmp[*,diff/2:dim+diff/2-1]
      tmp = congrid(reform(right[*,*,j]), dim, dim*fac, /center, cubic=-0.5)
      right[*,*,j] = tmp[*,diff/2:dim+diff/2-1]

    endfor

  endif

  ;-------------------------------------------------------------------------------------------

;   mwrfits, left, path+out+'_left.fits', hdrraw, /silent
;   mwrfits, right, path+out+'_right.fits', hdrraw, /silent
  writefits, path+out+'_left.fits', left, hdrraw
  writefits, path+out+'_right.fits', right, hdrraw

  ;-------------------------------------------------------------------------------------------

  writefits, path+out+'_paral.fits', reform(PA_onsky[i,*])
  writefits, path+out+'_paral.fits', reform(hang[i,*]), /append

; 
;   if (scflag eq 'y') then begin
; 
;     openw, lun, path+out+'_fctable.txt', width=2000, /get_lun
; 
;       printf, lun, '# Field center and angles table'
;       printf, lun, '# '
;       printf, lun, '# This table gives the field centers and rotation angles'
;       printf, lun, '# for a set of frames.'
;       printf, lun, '# Format is: centers in pixel coordinates (0,0) at lower left'
;       printf, lun, '#            angles in degrees ccw wrt. x-axis.'
;       printf, lun, '# '
;       printf, lun, '#         JD   Centre Left x    Left y   Right x   Right y    Angle (deg)'
; 
;       for j=0,sz3-1 do begin
; 
; 	printf, lun, newjd[i,j], center_left_x-pacx, center_left_y-pacy, center_right_x-pacx, center_right_y-pacy, PA_onsky[i,j], format='(f18.10, 4f10.2, f15.8)'
; 
;       endfor
; 
;     close, lun
;     free_lun, lun
; 
;   endif else begin
; 
;     openw, lun, path+out+'_fctable.txt', width=2000, /get_lun
; 
;       printf, lun, '# Field center and angles table'
;       printf, lun, '# '
;       printf, lun, '# This table gives the field centers and rotation angles'
;       printf, lun, '# for a set of frames.'
;       printf, lun, '# Format is: centers in pixel coordinates (0,0) at lower left'
;       printf, lun, '#            angles in degrees ccw wrt. x-axis.'
;       printf, lun, '# '
;       printf, lun, '#         JD   Centre Left x    Left y   Right x   Right y    Angle (deg)'
; 
;       for j=0,sz3-1 do begin
; 
; 	printf, lun, -99., -99., -99., -99., -99., PA_onsky[i,j], format='(f18.10, 4f10.2, f15.8)'
; 
;       endfor
; 
;     close, lun
;     free_lun, lun
; 
;   endelse

  ;-------------------------------------------------------------------------------------------

endfor


stop
end

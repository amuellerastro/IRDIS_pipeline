@resistant_mean.pro
@mrdfits.pro
@fxposit.pro
@mrd_hread.pro
@fxpar.pro
@gettok.pro
@valid_num.pro
@mrd_skip.pro
@get_eso_keyword.pro
@strsplit.pro
@readcol.pro
@remchar.pro
@strnumber.pro
@mwrfits.pro
@fxaddpar.pro
@headfits.pro
@fixpix_mod.pro
@strc.pro
@dist_circle.pro
@la_cosmic_MOD.pro
@reverse.pro
@array_indices.pro
@mpfit2dpeak.pro
@mpfit.pro
@mpfit2dfun.pro
;@shift_sub.pro
@fftshift.pro
@fxhmake.pro
@fxwrite.pro
@check_fits.pro
@daycnv.pro
@host_to_ieee.pro
@writefits.pro
@mkhdr.pro
@sxaddpar.pro
@get_date.pro
@sxdelpar.pro
@sxpar.pro

pro IRDIS_extract_reference_PSF

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

;selection of OBJECT,FLUX frames

idx = where(catg eq 'SCIENCE' and type eq 'OBJECT,FLUX')
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

  print, strcompress(i+1,/rem)+' ', tout[i], tcatg[i], ttype[i], tdit[i], tfilt[i], topti2[i], tcombind[i], ticor[i], format='(a3, a32, a15, a17, a5, a8, a8, a10, a13)'
  print, '----------------------------------------------------------------------------------------------------------------'

endfor

print, ''
sel = intarr(2)
read, 'Object,Flux (first and last number): ', sel
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
read, 'Master Dark: ', sel
bgf = tfile[sel-1]
; read, 'Bad Pixel Mask: ', sel
; bpf = tfile[sel-1]
bpf = path+'static_badpixels.fits'
; read, 'Master Flat: ', sel
; fff = tfile[sel-1]
fff = path+'instrument_flat.fits'
; read, 'Star Center: ', sel
; scf = tfile[sel-1]
; read, 'Sky BG: ', sel
; skf = tfile[sel-1]

; distquest = ''
; read, 'With distortion correction? y/n: ', distquest

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

bg = mrdfits(bgf, 0, hdrbg, /silent)
ff = mrdfits(fff, 0, hdrff, /silent)
bp = mrdfits(bpf, 0, hdrbg, /silent)
; sk = mrdfits(skf, 0, hdrsk, /silent)
; sc = mrdfits(scf, 1, hdrsc, /silent)

bgl = bg[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]
ffl = ff[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]
bpl = bp[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]
; skl = sk[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]

bgr = bg[1024:2047,*]
ffr = ff[1024:2047,*]
bpr = bp[1024:2047,*]
; skr = sk[1024:2047,*]

bgr = bgr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]
ffr = ffr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]
bpr = bpr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]
; skr = skr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]

;=============================================================================================

tmp = mrdfits(rawf[0], 0, hdrraw, /silent)
sz = size(tmp)
tmp = 0

if (n_elements(sz) eq 5) then sz3 = 1 else sz3 = sz[3]
; sz3 = 20

left = dblarr(dim,dim,sz3,nfiles)
right = dblarr(dim,dim,sz3,nfiles)

;-------------------------------------------------------------------------------------------
;get filter

opti2 = strcompress(get_eso_keyword(hdrraw,'HIERARCH ESO INS1 OPTI2 NAME'),/rem)
dit = strcompress(get_eso_keyword(hdrraw,'HIERARCH ESO DET SEQ1 DIT'),/rem)

filter_left = strmid(opti2,2,2)
filter_right = strmid(opti2,2,1)+strmid(opti2,4,1)

if (opti2 eq 'CLEAR') then begin

  bb = strcompress(get_eso_keyword(hdrraw, 'HIERARCH ESO INS COMB IFLT'),/rem)
  if (strmatch(bb, '*Y*') eq 1) then bb = 'Y'
  if (strmatch(bb, '*J*') eq 1) then bb = 'J'
  if (strmatch(bb, '*H*') eq 1) then bb = 'H'
  if (strmatch(bb, '*Ks*') eq 1) then bb = 'K'

  filter_left = bb+'L'
  filter_right = bb+'R'

endif

;-------------------------------------------------------------------------------------------

for i=0,nfiles-1 do begin

  im = mrdfits(rawf[i], 0, hdr, /silent)
  if (i eq 0) then hdrref = hdr

  for j=0,sz3-1 do begin

    print, ''
    print, 'Reducing frame '+strcompress(j+1,/rem)+' / '+strcompress(sz3,/rem)
    print, ''

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
    rawbgl = rawl-bgl
    rawbgffl = rawbgl/ffl
    fixpix_mod, rawbgffl, bpl, outim, npix=24, /weight;, /silent
    corrl = outim

    ;check for 'cosmic rays' or left overs
;     la_cosmic_MOD, corrl, outim, path, gain=1.75, readn=readn
;     crl = outim
    crl = corrl

  ;   idxnan = where(finite(crl) ne 1)
  ;   if (idxnan[0] ne -1) then begin
  ; 
  ;     print, 'NaNs present. Skipping CR correction for this frame.'
  ;     print, ''
  ;     crl = corrl
  ; 
  ;   endif


  ;right detector

    rawr = imr
    rawbgr = rawr-bgr
    rawbgffr = rawbgr/ffr
    fixpix_mod, rawbgffr, bpr, outim, npix=24, /weight;, /silent
    corrr = outim

    ;check for 'cosmic rays' or left overs
;     la_cosmic_MOD, corrr, outim, path, gain=1.75, readn=readn
;     crr = outim
    crr = corrr

  ;   idxnan = where(finite(crr) ne 1)
  ;   if (idxnan[0] ne -1) then begin
  ; 
  ;     print, 'NaNs present. Skipping CR correction for this frame.'
  ;     print, ''
  ;     crr = corrr
  ; 
  ;   endif

    imleft = crl
    imright = crr

;     if (distquest eq 'y') then begin
; 
;       imleft = poly_2d(imleft, kxl, kyl, 2, cubic=-0.5)
;       imright = poly_2d(imright, kxr, kyr, 2, cubic=-0.5)
; 
;     endif

    ;-------------------------------------------------------------------------------------------

    ;subtract sky

    ;imleft = imleft - skl
    ;imright = imright - skr

    ;-------------------------------------------------------------------------------------------

    ;removal of cross talk using medamp technique
    ;http://www.cfht.hawaii.edu/Instruments/Imaging/WIRCam/WIRCamCrosstalks.html

    ;amppos = linspace(0,1024,33);-abs(dx-dim/2)
    ;ampl = dblarr(32,1024,32) & ampr = ampl
    ;
    ;tmpl = dblarr(1024,1024)
    ;tmpl[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1,*] = imleft
    ;tmpr = dblarr(1024,1024)
    ;tmpr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1,*] = imright
    ;
    ;for k=0,31 do begin
    ;
    ;  ampl[*,*,k] = tmpl[amppos[k]:amppos[k+1]-1,*]
    ;  ampr[*,*,k] = tmpr[amppos[k]:amppos[k+1]-1,*]
    ;
    ;endfor
    ;
    ;mampl = median(ampl, /even, dimension=3)
    ;mampr = median(ampr, /even, dimension=3)
    ;
    ;for k=0,31 do begin
    ;
    ;  tmpl[k*32:k*32+31,*] = tmpl[k*32:k*32+31,*]-mampl
    ;  tmpr[k*32:k*32+31,*] = tmpr[k*32:k*32+31,*]-mampr
    ;
    ;endfor
    ;
    ;imleft = tmpl[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1]
    ;imright = tmpr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]

    ;-------------------------------------------------------------------------------------------

  ;   ;removal of fixed noise pattern / odd-even effect
  ; 
  ;   ;left detector upper and lower box median values
  ;   ;right detector upper and lower box median values
  ;   ll = median(imleft[*,89:139], dimension=2)
  ;   ul = median(imleft[*,664:714], dimension=2)
  ;   lr = median(imright[*,89:139], dimension=2)
  ;   ur = median(imright[*,664:714], dimension=2)
  ; 
  ;   for k=0,799 do begin
  ; 
  ;     imleft[k,*] = imleft[k,*]-mean([ll[k],ul[k]])
  ;     imright[k,*] = imright[k,*]-mean([lr[k],ur[k]])
  ; 
  ;   endfor

    ;-------------------------------------------------------------------------------------------

    left[*,*,j,i] = imleft
    right[*,*,j,i] = imright

    ;-------------------------------------------------------------------------------------------

  endfor

endfor

;===========================================================================================

;shift PSF of each frame to the center by fitting a Moffat profile

sleft = left
sright = right

rad = 25

for i=0,nfiles -1 do begin

  ;median combine images to identify max value in a robust way (averaging out bad pixels)

  iml = median(left[*,*,*,i], dimension=3, /even)
  imr = median(right[*,*,*,i], dimension=3, /even)

;   tmp = max(iml, locl)
;   idxl = array_indices(iml, locl)
;   tmp = max(imr, locr)
;   idxr = array_indices(imr, locr)
    find, iml, fdx, fdy, flux, sharp, roundness, 2*stddev(iml), 8., [-1,1], [0.2,1.0]
    if (n_elements(flux) gt 1) then begin
      
        tmp = max(flux,idxf)
        fdx = fdx[idxf]
        fdy = fdy[idxf]
      
    endif
    idxl = round([fdx,fdy])
    
    find, imr, fdx, fdy, flux, sharp, roundness, 2*stddev(imr), 8., [-1,1], [0.2,1.0]
    if (n_elements(flux) gt 1) then begin
      
        tmp = max(flux,idxf)
        fdx = fdx[idxf]
        fdy = fdy[idxf]
      
    endif
    idxr = round([fdx,fdy])
    
  
;     print, idxl, idxr
;   	writefits, '~/Downloads/iml.fits', iml
; 	writefits, '~/Downloads/imr.fits', imr
; 
; 	idxl = [439,363]
; 	idxr = [376,359]
	
; stop
; if i eq 0 then begin
; idxl = [422,422]
; idxr = [358,422]
; endif
; if i eq 1 then begin
; idxl = [437,368]
; idxr = [373,369]
; endif

for j=0,sz3-1 do begin

    cut_left = left[idxl[0]-rad:idxl[0]+rad, idxl[1]-rad:idxl[1]+rad,j,i]
    cut_right = right[idxr[0]-rad:idxr[0]+rad, idxr[1]-rad:idxr[1]+rad,j,i]

    xa = dindgen(2*rad+1)+1.d0 & ya = xa

    ;left
    estimates = [median(cut_left), max(cut_left), 2., 2., rad, rad, 0., 1.]

;     weights = cut_left
;     idx1 = where(cut_left eq 0.)	;e.g. beta Pic close to the center which is masked out
;     if (idx1[0] ne -1) then begin
;       idx2 = array_indices(weights, idx1)
;       for k=0,n_elements(idx2[0,*])-1 do weights[idx2[0,k],idx2[1,k]] = 1.d5
;     endif
;     sign = signum(cut_left)
;     idx1 = where(sign eq -1.)
;     if (idx1[0] ne -1) then begin
;       idx2 = array_indices(weights, idx1)
;       for k=0,n_elements(idx2[0,*])-1 do weights[idx2[0,k],idx2[1,k]] = 1.d5
;     endif
;     weights = 1./sqrt(weights)

    yfit = mpfit2dpeak(cut_left, A, xa, ya, /moffat, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /tilt, /quiet);, weights=weights)

    xpos = A[4]-rad+idxl[0]
    ypos = A[5]-rad+idxl[1]

    ;sleft[*,*,j,i] = shift_sub(left[*,*,j,i], dim/2-xpos+1., dim/2-ypos+1.)
    wframe = 100.	;width of frame
    tmp = dblarr(dim+2.*wframe, dim+2.*wframe)
    tmp[wframe:wframe+dim-1, wframe:wframe+dim-1] = left[*,*,j,i]
    stmp = fftshift(tmp, dim/2-xpos+1., dim/2-ypos+1.)
    sleft[*,*,j,i] = stmp[(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1,(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1]

    ;stmp = shift_sub(tmp, dim/2-xpos+1., dim/2-ypos+1.)
    ;sleft[*,*,j,i] = stmp[(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1,(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1]


    ;right
    estimates = [median(cut_right), max(cut_right), 2., 2., rad, rad, 0., 1.]

;     weights = cut_right
;     idx1 = where(cut_right eq 0.)	;e.g. beta Pic close to the center which is masked out
;     if (idx1[0] ne -1) then begin
;       idx2 = array_indices(weights, idx1)
;       for k=0,n_elements(idx2[0,*])-1 do weights[idx2[0,k],idx2[1,k]] = 1.d5
;     endif
;     sign = signum(cut_right)
;     idx1 = where(sign eq -1.)
;     if (idx1[0] ne -1) then begin
;       idx2 = array_indices(weights, idx1)
;       for k=0,n_elements(idx2[0,*])-1 do weights[idx2[0,k],idx2[1,k]] = 1.d5
;     endif
;     weights = 1./sqrt(weights)

    yfit = mpfit2dpeak(cut_right, A, xa, ya, /moffat, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /tilt, /quiet);, weights=weights)

    xpos = A[4]-rad+idxr[0]
    ypos = A[5]-rad+idxr[1]

    ;sright[*,*,j,i] = shift_sub(right[*,*,j,i], dim/2-xpos+1., dim/2-ypos+1.)
    wframe = 100.	;width of frame
    tmp = dblarr(dim+2.*wframe, dim+2.*wframe)
    tmp[wframe:wframe+dim-1, wframe:wframe+dim-1] = right[*,*,j,i]
    stmp = fftshift(tmp, dim/2-xpos+1., dim/2-ypos+1.)
    sright[*,*,j,i] = stmp[(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1,(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1]

    ;stmp = shift_sub(tmp, dim/2-xpos+1., dim/2-ypos+1.)
    ;sright[*,*,j,i] = stmp[(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1,(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1]

  endfor

endfor

window, 2, xs=1000, ys=800
!p.multi=[0,1,2]
plot, sleft[dim/2,dim/2,*,0], psym=2, yr=[0,ceil(max(sleft))]
if (nfiles gt 1) then oplot, sleft[dim/2,dim/2,*,1], psym=2, color=cgcolor('red')
plot, sright[dim/2,dim/2,*,0], psym=2, yr=[0,ceil(max(sright))]
if (nfiles gt 1) then oplot, sright[dim/2,dim/2,*,1], psym=2, color=cgcolor('red')
!p.multi=[0,1,0]

;===========================================================================================

;reject bad frames based on statistics of left detector, right one should be correlated
;compute mean PSF

sdev = dblarr(nfiles)
mflux = dblarr(nfiles)
psf_left = dblarr(dim, dim, nfiles)
psf_right = dblarr(dim, dim, nfiles)

for i=0,nfiles-1 do begin

;   sdev[i] = stddev(sleft[dim/2,dim/2,*,i])
;   mflux[i] = median(sleft[dim/2,dim/2,*,i], /even)
;   idxgood = where(sleft[dim/2,dim/2,*,i] gt mflux[i]-sdev[i])

  medim = median(sleft[*,*,*,i], dim=3)
  sim = dblarr(sz3)
  for j=0,sz3-1 do sim[j] = stddev(reform(sleft[*,*,j,i])-medim)
  resistant_mean, sim, 3, t1, t2, nbad, goodvec=idxgood, badvec=idxbad

  psf_left[*,*,i] = mean(sleft[*,*,*,i], dim=3)
  psf_right[*,*,i] = mean(sright[*,*,*,i], dim=3)
  
;  psf_left[*,*,i] = mean(sleft[*,*,idxgood,i], dim=3)
;  psf_right[*,*,i] = mean(sright[*,*,idxgood,i], dim=3)

  rad = 13

  writefits, path+'PSF_'+filter_left+'_'+strcompress(i+1,/rem)+'.fits', sleft[dim/2-rad:dim/2+rad, dim/2-rad:dim/2+rad,idxgood,i]
  writefits, path+'PSF_'+filter_right+'_'+strcompress(i+1,/rem)+'.fits', sright[dim/2-rad:dim/2+rad, dim/2-rad:dim/2+rad,idxgood,i]


;   print, ''
;   print, 'Rejected '+strcompress(sz3-n_elements(idxgood),/rem)+' out of '+strcompress(sz3,/rem)

;   window, i
;   xa = dindgen(sz3)
;   plot, xa, sim, /yn, title=star
;   oplot, xa[idxgood], sim[idxgood], psym=2, color=cgcolor('green')
;   if (idxbad[0] ne -1) then oplot, xa[idxbad], sim[idxbad], psym=2, color=cgcolor('red')


endfor

;===========================================================================================

;average frames and cut out PSF

; rad = 13	;should be always a odd value
; 
; ;write out intermediate PSFs
; for i=0,nfiles-1 do begin
; 
;   writefits, path+'PSF_'+filter_left+'_'+strcompress(i+1,/rem)+'.fits', sleft[dim/2-rad:dim/2+rad, dim/2-rad:dim/2+rad,idxgood,i]
;   writefits, path+'PSF_'+filter_right+'_'+strcompress(i+1,/rem)+'.fits', sright[dim/2-rad:dim/2+rad, dim/2-rad:dim/2+rad,idxgood,i]
; 
; ;   mwrfits, psf_left[dim/2-rad:dim/2+rad, dim/2-rad:dim/2+rad,i], path+'PSF_'+filter_left+'_'+strcompress(i+1,/rem)+'.fits'
; ;   mwrfits, psf_right[dim/2-rad:dim/2+rad, dim/2-rad:dim/2+rad,i], path+'PSF_'+filter_right+'_'+strcompress(i+1,/rem)+'.fits'
; ;   writefits, path+'PSF_'+filter_left+'_'+strcompress(i+1,/rem)+'.fits', psf_left[dim/2-rad:dim/2+rad, dim/2-rad:dim/2+rad,i]
; ;   writefits, path+'PSF_'+filter_right+'_'+strcompress(i+1,/rem)+'.fits', psf_right[dim/2-rad:dim/2+rad, dim/2-rad:dim/2+rad,i]
; 
; endfor

if (nfiles gt 1) then begin

  fin_psf_left = mean(psf_left, dim=3)
  fin_psf_right = mean(psf_right, dim=3)

endif else begin

  fin_psf_left = psf_left
  fin_psf_right = psf_right

endelse


fin_psf_left = fin_psf_left[dim/2-rad:dim/2+rad, dim/2-rad:dim/2+rad]
fin_psf_right = fin_psf_right[dim/2-rad:dim/2+rad, dim/2-rad:dim/2+rad]


; mwrfits, fin_psf_left, path+'PSF_'+filter_left+'.fits'
; mwrfits, fin_psf_right, path+'PSF_'+filter_right+'.fits'
writefits, path+'PSF_'+filter_left+'.fits', fin_psf_left, hdrref
writefits, path+'PSF_'+filter_right+'.fits', fin_psf_right, hdrref

stop

;===========================================================================================

; 
; ;median combine frames
; 
; if (sz3 gt 1) then begin
; 
;   iml = median(left, dimension=3, /even)
;   imr = median(right, dimension=3, /even)
; 
; endif else begin
; 
;   iml = left
;   imr = right
; 
; endelse
; 
; ; iml = median(left, dimension=3, /even)
; ; imr = median(right, dimension=3, /even)
; 
; ;-------------------------------------------------------------------------------------------
; 
; ;cut out PSF and write file
; 
; rad = 25	;should be always a odd value
; 
; ;left detector
; 
; dum = max(iml, loc)
; idxl = array_indices(iml, loc)
; psf_left = iml[idxl[0]-2.*rad:idxl[0]+2.*rad, idxl[1]-2.*rad:idxl[1]+2.*rad]
; 
; xa = dindgen(n_elements(psf_left[0,*])) & ya=xa
; estimates = [median(psf_left,/even), max(psf_left), 4., 4., 2.*rad, 2.*rad, 0., 1.]
; pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},8)
; pi[1].limited(0) = 1
; pi[1].limits(0) = 1.d-15	;force amplitude to be positive
; pi[6].limited(0) = 1
; pi[6].limits(0) = 0.d0
; pi[6].limited(1) = 1
; pi[6].limits(1) = 2.d0*!DPI
; ; pi[2].limited(0) = 1
; ; pi[2].limits(0) = 2.d0
; ; pi[2].limited(1) = 1
; ; pi[2].limits(1) = 10.d0
; ; pi[3].limited(0) = 1
; ; pi[3].limits(0) = 2.d0
; ; pi[3].limited(1) = 1
; ; pi[3].limits(1) = 10.d0
; yfit = mpfit2dpeak(psf_left, A, xa, ya, /moffat, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /tilt, parinfo=pi, /quiet)
; 
; spsfl = shift_sub(psf_left, 2.*rad-A[4], 2.*rad-A[5])
; psf_left = spsfl[2.*rad-rad:2.*rad+rad, 2.*rad-rad:2.*rad+rad]
; 
; ;right detector
; 
; dum = max(imr, loc)
; idxr = array_indices(imr, loc)
; psf_right = imr[idxr[0]-2.*rad:idxr[0]+2.*rad, idxr[1]-2.*rad:idxr[1]+2.*rad]
; 
; xa = dindgen(n_elements(psf_right[0,*])) & ya=xa
; estimates = [median(psf_right,/even), max(psf_right), 4., 4., 2.*rad, 2.*rad, 0., 1.]
; pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},8)
; pi[1].limited(0) = 1
; pi[1].limits(0) = 1.d-15	;force amplitude to be positive
; pi[6].limited(0) = 1
; pi[6].limits(0) = 0.d0
; pi[6].limited(1) = 1
; pi[6].limits(1) = 2.d0*!DPI
; ; pi[2].limited(0) = 1
; ; pi[2].limits(0) = 2.d0
; ; pi[2].limited(1) = 1
; ; pi[2].limits(1) = 10.d0
; ; pi[3].limited(0) = 1
; ; pi[3].limits(0) = 2.d0
; ; pi[3].limited(1) = 1
; ; pi[3].limits(1) = 10.d0
; yfit = mpfit2dpeak(psf_right, A, xa, ya, /moffat, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /tilt, parinfo=pi, /quiet)
; 
; spsfr = shift_sub(psf_right, 2.*rad-A[4], 2.*rad-A[5])
; psf_right = spsfr[2.*rad-rad:2.*rad+rad, 2.*rad-rad:2.*rad+rad]
; 
; 
; ; dum = max(iml, loc)
; ; idxl = array_indices(iml, loc)
; ; psf_left = iml[idxl[0]-rad:idxl[0]+rad, idxl[1]-rad:idxl[1]+rad]
; ; 
; ; dum = max(imr, loc)
; ; idxl = array_indices(imr, loc)
; ; psf_right = imr[idxl[0]-rad:idxl[0]+rad, idxl[1]-rad:idxl[1]+rad]
; 
; 
; mwrfits, psf_left, path+'PSF_'+filter_left+'.fits'
; mwrfits, psf_right, path+'PSF_'+filter_right+'.fits'



stop
end

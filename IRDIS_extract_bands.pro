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
@detabify.pro
@fxparpos.pro
@dist.pro
@circint_MOD.pro
@cgcolor.pro
@cggetcolorstate.pro
@cgsnapshot.pro
@cgcolor24.pro
@robust_mean.pro
@avg.pro
@writefits.pro
@mkhdr.pro
@sxaddpar.pro
@get_date.pro
@daycnv.pro
@sxdelpar.pro
@sxpar.pro
@fits_add_checksum.pro
@checksum32.pro
@n_bytes.pro
@is_ieee_big.pro
@host_to_ieee.pro
@fits_ascii_encode.pro

pro IRDIS_extract_bands, qc=qc


;================================================================================

readcol, '/home/amueller/work/IDLlibs/AO/SPHERE/datapaths.txt', tmp, format='a', /silent
print, ''
for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]

print, ''
read, 'Select Path: ', selp
path = tmp[selp-1]+'/RAW/IRDIS/'

;================================================================================

file_left = file_search(path+'SPHER*_left.fits', count=nleft)
file_right = file_search(path+'SPHER*_right.fits', count=nright)

if (nleft ne nright) then stop

;file_para = file_search(path+'SPHER*_fctable.txt', count=ntable)
file_para = file_search(path+'SPHER*_paral.fits', count=ntable)

if (ntable ne nright) then stop

;================================================================================

;read in dummy to get array sizes and filter
dum = mrdfits(file_left[0], 0, hdr, /silent)
sz = size(dum)

if (n_elements(sz) eq 5) then sz3 = 1 else sz3 = sz[3]

para = dblarr(sz3*nleft)
ha = dblarr(sz3*nleft)

;output file names for later
filt = get_eso_keyword(hdr, 'HIERARCH ESO INS COMB IFLT')

if (strmatch(filt, '*J23*') eq 1) then begin

  fn_left = path+'img_J2_dc.fits'
  fn_right = path+'img_J3_dc.fits'
  fn_para_left = path+'vec_J2_paral.fits'
  fn_para_right = path+'vec_J3_paral.fits'

  l1 = 1.1895d-6
  l2 = 1.2698d-6

endif

if (strmatch(filt, '*H23*') eq 1) then begin

  fn_left = path+'img_H2_dc.fits'
  fn_right = path+'img_H3_dc.fits'
  fn_para_left = path+'vec_H2_paral.fits'
  fn_para_right = path+'vec_H3_paral.fits'

  l1 = 1.5888d-6
  l2 = 1.6671d-6

endif

if (strmatch(filt, '*K12*') eq 1) then begin

  fn_left = path+'img_K1_dc.fits'
  fn_right = path+'img_K2_dc.fits'
  fn_para_left = path+'vec_K1_paral.fits'
  fn_para_right = path+'vec_K2_paral.fits'

  l1 = 2.1025d-6
  l2 = 2.2550d-6

endif

if (strmatch(filt, '*BB_Y*') eq 1) then begin

  fn_left = path+'img_YL_dc.fits'
  fn_right = path+'img_YR_dc.fits'
  fn_para_left = path+'vec_YL_paral.fits'
  fn_para_right = path+'vec_YR_paral.fits'

  l1 = 1.0425d-6
  l2 = l1

endif

if (strmatch(filt, '*BB_J*') eq 1) then begin

  fn_left = path+'img_JL_dc.fits'
  fn_right = path+'img_JR_dc.fits'
  fn_para_left = path+'vec_JL_paral.fits'
  fn_para_right = path+'vec_JR_paral.fits'

  l1 = 1.2575d-6
  l2 = l1

endif

if (strmatch(filt, '*BB_H*') eq 1) then begin

  fn_left = path+'img_HL_dc.fits'
  fn_right = path+'img_HR_dc.fits'
  fn_para_left = path+'vec_HL_paral.fits'
  fn_para_right = path+'vec_HR_paral.fits'

  l1 = 1.6255d-6
  l2 = l1

endif

if (strmatch(filt, '*BB_K*') eq 1) then begin

  fn_left = path+'img_KL_dc.fits'
  fn_right = path+'img_KR_dc.fits'
  fn_para_left = path+'vec_KL_paral.fits'
  fn_para_right = path+'vec_KR_paral.fits'

  l1 = 2.1813d-6
  l2 = l1

endif


plsc = 0.01225d0
diam = 8.2
fwhm1 = l1/diam*206265.d0/plsc
fwhm2 = l2/diam*206265.d0/plsc

;================================================================================

im_left = dblarr(sz[1], sz[2], sz3*nleft)
; im_left = fltarr(sz[1], sz[2], sz3*nleft)	;in case of low memory

time = strarr(nleft)	;extract time stamp to sync with ascii files containing parallactic angle
for i=0,nleft-1 do begin

  pos = strpos(file_left[i], '_', /reverse_search)
  time[i] = strmid(file_left[i], 0, strlen(file_left[i])-pos)

  st = mrdfits(file_left[i], 0, hdr, /silent)
  if (i eq 0) then hdrscirefleft = hdr
  im_left[*,*,i*sz3:sz3*(i+1)-1] = st

  tmp = mrdfits(file_para[i], 0, /silent)
  para[i*sz3:sz3*(i+1)-1] = tmp
  tmp = mrdfits(file_para[i], 1, /silent)
  ha[i*sz3:sz3*(i+1)-1] = tmp

;   readcol, file_para[i], t1, t2, t3, t4, t5, t6, format='d,d,d,d,d,d', numline=sz3, skipline=8, count=nl, /silent
;   if (nl ne sz3) then begin
; 
;     print, ''
;     print, 'Check if header of text file containing parallactic angle changed. Stop.'
;     stop 
; 
;   endif else begin
; 
;     para[i*sz3:sz3*(i+1)-1] = t6
; 
;   endelse

endfor

;================================================================================

if keyword_set(qc) then begin

  dim = (size(im_left))[1]
  diff = dblarr(dim,dim,n_elements(im_left[0,0,*]))
  medim = median(im_left, dim=3)
  totd = dblarr(n_elements(im_left[0,0,*]))
  totall = totd
  for i=0,n_elements(im_left[0,0,*])-1 do begin

    scale = total(im_left[*,*,i])/total(medim)
    diff[*,*,i] = im_left[*,*,i]-medim*scale
    totd[i] = total(abs(diff[*,*,i]))
    totall[i] = total(im_left[*,*,i])/scale

  endfor

  dum = robust_mean(totd, 2., sigma, rej, goodind=idxgood0)

  ;================================================================================

  ;quality control

  dim = (size(im_left))[1]
  xc = dim/2. & yc = xc
  nframes = (size(im_left))[3]
  sdev1 = dblarr(nframes)
  sdev2 = dblarr(nframes)
  flux = dblarr(nframes)
  mflux = dblarr(nframes)
  annflux1 = dblarr(nframes)
  annflux2 = dblarr(nframes)

  mask_t = shift(dist(dim), xc, yc)

  ;left detector

  mask1 = mask_t ge 3.*fwhm1 and mask_t le 8.*fwhm1
  mask1 = mask1[*]	;convert to 1D vector
  idxmask1 = where(mask1 eq 1)

  mask2 = mask_t ge 28.*fwhm1 and mask_t le 40.*fwhm1
  mask2 = mask2[*]	;convert to 1D vector
  idxmask2 = where(mask2 eq 1)

  for i=0,nframes-1 do begin

    ;measure stddev between 1-4 lambda/D, see Absil 2013
    tmp = im_left[*,*,i]
    tmp = tmp[*]
    sdev1[i] = stddev(tmp[idxmask1])
    sdev2[i] = stddev(tmp[idxmask2])

    annflux1[i] = total(tmp[idxmask1])
    annflux2[i] = total(tmp[idxmask2])

    ;measure flux, additional criteria, is more robust
    tmp = im_left[*,*,i]
    circint_MOD, tmp, xc, yc, 2.*fwhm1, tot, mtot, meantot, maxpx, sdtot, npix, totdif, npixdif, t8
    flux[i] = tot
    mflux[i] = mtot

  endfor

  ;reject bad frames
  idxgood1 = where(sdev1 lt median(sdev1)+stddev(sdev1))
  idxgood2 = where(sdev2 lt median(sdev2)+stddev(sdev2))
  idxgood3 = where(mflux lt median(mflux)+stddev(mflux))

  flag = intarr(nframes)

  for i=0,nframes-1 do begin

    idx0 = where(i eq idxgood0)
    if (idx0[0] eq -1) then flag[i] = 1
    idx1 = where(i eq idxgood1)
    if (idx1[0] eq -1) then flag[i] = 1
    idx2 = where(i eq idxgood2)
    if (idx2[0] eq -1) then flag[i] = 1
    idx3 = where(i eq idxgood3)
    if (idx3[0] eq -1) then flag[i] = 1

  endfor

  idxgood = where(flag eq 0)
  idxbad = where(flag eq 1)

  used_frames = [n_elements(idxgood), nframes]

  print, strcompress(uint(nframes-n_elements(idxgood)),/rem)+' / '+strcompress(uint(nframes),/rem)+' frames rejected'
  print, idxbad

  ;   nframes = n_elements(im_left[0,0,*])
    window, 0, xs=800, ys=1000
    !p.multi=[0,1,3]
    plot, dindgen(nframes)+1, sdev1, /yn, psym=2, xtitle='Frame', ytitle='SDEV 1st annuli', charsize=2
      oplot, !x.crange, median(sdev1)*[1,1]
      oplot, !x.crange, (median(sdev1)+stddev(sdev1))*[1,1], color=cgcolor('red')
      oplot, idxgood+1, sdev1[idxgood], psym=2, color=cgcolor('green')
    plot, dindgen(nframes)+1, sdev2, /yn, psym=2, xtitle='Frame', ytitle='SDEV 2nd annuli', charsize=2
      oplot, !x.crange, median(sdev2)*[1,1]
      oplot, !x.crange, (median(sdev2)+stddev(sdev2))*[1,1], color=cgcolor('red')
      oplot, idxgood+1, sdev2[idxgood], psym=2, color=cgcolor('green')
    plot, dindgen(nframes)+1, mflux, /yn, psym=2, xtitle='Frame', ytitle='Median Flux', charsize=2
      oplot, !x.crange, median(mflux)*[1,1]
      oplot, !x.crange, (median(mflux)+stddev(mflux))*[1,1], color=cgcolor('red')
      oplot, idxgood+1, mflux[idxgood], psym=2, color=cgcolor('green')
    !p.multi=[0,1,0]

    window, 2
    plot, dindgen(n_elements(im_left[0,0,*]))+1., totd, /yn
      oplot, dindgen(n_elements(im_left[0,0,*]))+1., totd, psym=2
      oplot, idxgood0+1, totd[idxgood0], psym=2, color=cgcolor('yellow')
      oplot, idxgood+1, totd[idxgood], psym=2, color=cgcolor('green')

  ;output
;   mwrfits, im_left[*,*,idxgood], fn_left, /silent
;   mwrfits, para[idxgood], fn_para_left, /silent
;   mwrfits, im_left, fn_left+'.orig', /silent
;   mwrfits, para, fn_para_left+'.orig', /silent
;   mwrfits, im_left[*,*,idxbad], fn_left+'.bad', /silent
;   mwrfits, para[idxbad], fn_para_left+'.bad', /silent
  writefits, fn_left, im_left[*,*,idxgood]
  writefits, fn_para_left, para[idxgood]
  writefits, fn_para_left, ha[idxgood], /append
  writefits, fn_left+'.orig', im_left
  writefits, fn_para_left+'.orig', para
  writefits, fn_para_left+'.orig', ha, /append
  writefits, fn_left+'.bad', im_left[*,*,idxbad]
  writefits, fn_para_left+'.bad', para[idxbad]
  writefits, fn_para_left+'.bad', ha[idxbad], /append

endif else begin

  writefits, fn_left, im_left, hdrscirefleft
  writefits, fn_para_left, para
  writefits, fn_para_left, ha, /append

endelse



im_left = 0	;free up memory

;================================================================================

im_right = dblarr(sz[1], sz[2], sz3*nright)
; im_right = fltarr(sz[1], sz[2], sz3*nright)	;in case of low memory

time = strarr(nright)	;extract time stamp to sync with ascii files containing parallactic angle
for i=0,nleft-1 do begin

  pos = strpos(file_right[i], '_', /reverse_search)
  time[i] = strmid(file_right[i], 0, strlen(file_right[i])-pos)

  st = mrdfits(file_right[i], 0, hdr, /silent)
  if (i eq 0) then hdrscirefright = hdr
  im_right[*,*,i*sz3:sz3*(i+1)-1] = st

  tmp = mrdfits(file_para[i], 0, /silent)
  para[i*sz3:sz3*(i+1)-1] = tmp
  tmp = mrdfits(file_para[i], 1, /silent)
  ha[i*sz3:sz3*(i+1)-1] = tmp

;   readcol, file_para[i], t1, t2, t3, t4, t5, t6, format='d,d,d,d,d,d', numline=sz3, skipline=8, count=nl, /silent
;   if (nl ne sz3) then begin
; 
;     print, ''
;     print, 'Check if header of text file containing parallactic angle changed. Stop.'
;     stop 
; 
;   endif else begin
; 
;     para[i*sz3:sz3*(i+1)-1] = t6
; 
;   endelse

endfor

;output

if keyword_set(qc) then begin

;   mwrfits, im_right, fn_right+'.orig', /silent
;   mwrfits, para, fn_para_right+'.orig', /silent
;   mwrfits, im_right[*,*,idxgood], fn_right, /silent
;   mwrfits, para[idxgood], fn_para_right, /silent
  writefits, fn_right+'.orig', im_right
  writefits, fn_para_right+'.orig', para
  writefits, fn_para_right+'.orig', ha, /append
  writefits, fn_right, im_right[*,*,idxgood]
  writefits, fn_para_right, para[idxgood]
  writefits, fn_para_right, ha[idxgood], /append

endif else begin

;   mwrfits, im_right, fn_right, /silent
;   mwrfits, para, fn_para_right, /silent
  writefits, fn_right, im_right, hdrscirefright
  writefits, fn_para_right, para
  writefits, fn_para_right, ha, /append

endelse

im_right = 0	;free up memory

;================================================================================

stop
end
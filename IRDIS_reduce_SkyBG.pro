@la_cosmic_MOD.pro
@headfits.pro
@fxposit.pro
@mrd_hread.pro
@get_eso_keyword.pro
@strsplit.pro
@mrdfits.pro
@fxpar.pro
@gettok.pro
@valid_num.pro
@mrd_skip.pro
@fixpix_mod.pro
@strc.pro
@dist_circle.pro
; @fshift.pro
@mwrfits.pro
@fxaddpar.pro
@reverse.pro
@readcol.pro
@remchar.pro
@sigfig.pro
@fxhmake.pro
@fxwrite.pro
@check_fits.pro
@get_date.pro
@daycnv.pro
@host_to_ieee.pro
@writefits.pro
@mkhdr.pro
@sxaddpar.pro
@get_date.pro
@daycnv.pro
@sxdelpar.pro
@sxpar.pro

pro IRDIS_reduce_SkyBG

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

;selection of SKY frame

idx = where(catg eq 'SCIENCE' and type eq 'SKY')
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
read, 'SKy BG (first and last number): ', sel
skf = tfile[indgen(sel[1]-sel[0]+1)+sel[0]-1]
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

;====================================================================================================

;cut images down to 800x800px
;using bad pixel map of FF too and all the edges are full with bad pixels -> takes a very long time for fixpix_mod
; xc = 512
; yc = 512
dim = 800

;====================================================================================================

; bg = mrdfits(bgf, 0, /silent)
; ff = mrdfits(fff, 0, /silent)
bp = mrdfits(bpf, 0, /silent)

; bgl = bg[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1,*]
; ffl = ff[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1,*]
bpl = bp[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1,*]

; bgr = bg[1024:2047,*]
; ffr = ff[1024:2047,*]
bpr = bp[1024:2047,*]

; bgr = bgr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1,*]
; ffr = ffr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]
bpr = bpr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1]

;====================================================================================================

nfiles = n_elements(skf)

; pacx = dblarr(nfiles)
; pacy = dblarr(nfiles)
sk = dblarr(2048, 1024, nfiles)

;read in dummy file in case NDIT>1
dum = mrdfits(skf[0], 0, /silent)
if (n_elements(size(dum)) gt 5) then nditflag = 1 else nditflag = 0

for i=0,nfiles-1 do begin

  if (nditflag eq 1) then begin

    tmp = mrdfits(skf[i], 0, hdrsk, /silent)
    sk[*,*,i] = median(tmp, dimension=3, /even)


  endif else begin

    sk[*,*,i] = mrdfits(skf[i], 0, hdrsk, /silent)

  endelse


  if (i eq 0) then begin

    dit = double(get_eso_keyword(hdrsk,'HIERARCH ESO DET SEQ1 DIT'))
    exptime = strcompress(get_eso_keyword(hdrsk, 'EXPTIME'),/rem)

  endif

;   pacx[i] = double(get_eso_keyword(hdrsk, 'HIERARCH ESO INS1 PAC X'))/18.d0
;   pacy[i] = double(get_eso_keyword(hdrsk, 'HIERARCH ESO INS1 PAC Y'))/18.d0

endfor

skl = sk[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1,*]

skr = sk[1024:2047,*,*]
skr = skr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1,*]

nim = nfiles	;currently ignoring if there are nexp>1 per file

;====================================================================================================

;bad pixel correction of images and division by FF

skcl = skl
skcr = skr

if (dit gt 10.) then readn = 10. else readn = 4.	;DIT of 10s chosen abritrarily

for i=0,nim-1 do begin

;   tmp = (skl[*,*,i]-bgl)/ffl
;   tmp = (skl[*,*,i])/ffl
  tmp = skl[*,*,i]
  fixpix_mod, tmp, bpl, outim, npix=24, /weight;, /silent
  skcl[*,*,i] = outim
;   la_cosmic_MOD, skcl[*,*,i], outim, path, gain=1.75, readn=readn
;   skcl[*,*,i] = outim

;   tmp = (skr[*,*,i]-bgr)/ffr
;   tmp = (skr[*,*,i])/ffr
  tmp = skr[*,*,i]
  fixpix_mod, tmp, bpr, outim, npix=24, /weight;, /silent
  skcr[*,*,i] = outim
;   la_cosmic_MOD, skcr[*,*,i], outim, path, gain=1.75, readn=readn
;   skcr[*,*,i] = outim

;   sl = shift_sub(skcl[*,*,i], pacx[i], pacy[i])
;   skcl[*,*,i] = sl
;   sr = shift_sub(skcr[*,*,i], pacx[i], pacy[i])
;   skcr[*,*,i] = sr


endfor

;====================================================================================================

;put them back in the oringial array

sk_new = dblarr(2048, 1024, nim)

sk_new[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1,*] = skcl
sk_new[1024+dx[1]-dim/2:1024+dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1,*] = skcr

;====================================================================================================

;combine them and save file

sk_final = dblarr(2048,1024,nim)

if (nim gt 1) then sk_final = median(sk_new, dimension=3, /even) else sk_final = sk_new

exptime = sigfig(exptime,3)
; mwrfits, sk_final, path+'sky_background_'+exptime+'s.fits', /silent
writefits, path+'sky_background_'+exptime+'s.fits', sk_final


stop
end

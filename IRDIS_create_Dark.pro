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
@mwrfits.pro
@fxaddpar.pro
@sigfig.pro
@arrdelete.pro
@writefits.pro
@mkhdr.pro
@sxaddpar.pro
@get_date.pro
@daycnv.pro
@sxdelpar.pro
@sxpar.pro
@check_fits.pro
@fits_add_checksum.pro
@checksum32.pro
@n_bytes.pro
@is_ieee_big.pro
@host_to_ieee.pro
@fits_ascii_encode.pro

pro IRDIS_create_Dark

;===============================================================================================

readcol, '/home/amueller/work/IDLlibs/AO/SPHERE/datapaths.txt', tmp, format='a', /silent
print, ''
for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]

print, ''
read, 'Select Path: ', selp
path = tmp[selp-1]+'/RAW/IRDIS/'

;path = '/home/amueller/Downloads/betaPic/'
;path = '/home/amueller/Downloads/Sirius/'

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

;selection of Dark frames

idx = where(catg eq 'CALIB' and strmatch(type, '*DARK*') eq 1)
if (idx[0] eq -1) then begin

  print, 'No dark found.'
  return

endif

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
read, 'Raw Dark: ', sel

file = tfile[sel-1]

;===============================================================================================



im = mrdfits(file, 0, hdr, /silent)
sz = size(im)
if (n_elements(sz) eq 5) then sz3 = 1 else sz3 = sz[3]

exptime = strcompress(get_eso_keyword(hdr, 'EXPTIME'),/rem)

;exptime = strmid(exptime, 0, strpos(exptime, '.'))
exptime = sigfig(exptime,3)

if (sz3 eq 1) then begin

;   mwrfits, im, path+'master_dark_'+exptime+'s.fits', hdr
  writefits, path+'master_dark_'+exptime+'s.fits', im, hdr

endif else begin

  dark = median(im, dimension=3)

  for j=0,n_elements(hdr)-1 do begin

    if (strmatch(hdr[j], '*NAXIS3*') eq 1) then begin

      newhdr = arrdelete(hdr, at=j, length=1)
      flag = 1

    endif

  endfor

  if (flag eq 1) then hdr = newhdr

;   mwrfits, dark, path+'master_dark_'+exptime+'s.fits', hdr
  writefits, path+'master_dark_'+exptime+'s.fits', dark, hdr

endelse


stop
end
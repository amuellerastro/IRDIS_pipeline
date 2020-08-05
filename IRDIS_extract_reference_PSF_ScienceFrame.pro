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
@arrdelete.pro
@mwrfits.pro
@fxaddpar.pro

pro IRDIS_extract_reference_PSF_ScienceFrame

;=============================================================================================

readcol, '/home/amueller/work/IDLlibs/AO/SPHERE/datapaths.txt', tmp, format='a', /silent
print, ''
for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]

print, ''
read, 'Select Path: ', selp
path = tmp[selp-1]+'/RAW/IRDIS/'

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

;selection of OBJECT,FLUX frames

idx = where(catg eq 'SCIENCE' and type eq 'OBJECT')
;idx = where(catg eq 'CALIB' and type eq 'OBJECT,ASTROMETRY')
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
read, 'Select ONE science file: ', sel
rawf = tfile[sel-1]

;=============================================================================================

;read in science file for the header only

dum = mrdfits(rawf, 0, hdrsci, /silent)
dum = 0

;have to remove naxis keyword, othwerwise mrdfits complains
for i=0,n_elements(hdrsci)-1 do begin

  if (strmatch(hdrsci[i], '*NAXIS3*') eq 1) then newhdr = arrdelete(hdrsci, at=i, length=1)

endfor

;=============================================================================================

;file selection

file = file_search(path+'img_*_dc.fits', count=nfiles)

if (nfiles ne 2) then stop	;assuming 2 filters

rad = 25	;should be always a odd value

;=============================================================================================

;1st file

im = mrdfits(file[0], 0, /silent)

;extract file name
pos = strpos(file[0], '/img_', /reverse_search)
fn = strmid(file[0], pos1+5, 2)

dim = (size(im))[1]
xc = dim/2. & yc = xc

im = im[xc-rad:xc+rad, yc-rad:yc+rad,*]

if (n_elements(size(im)) eq 5) then mim = im $
  else mim = median(im, dimension=3)

mwrfits, mim, path+'PSF_'+fn+'.fits', newhdr

;=============================================================================================

;2nd file

im = mrdfits(file[1], 0, /silent)

;extract file name
pos = strpos(file[1], '/img_', /reverse_search)
fn = strmid(file[1], pos1+5, 2)

dim = (size(im))[1]
xc = dim/2. & yc = xc

im = im[xc-rad:xc+rad, yc-rad:yc+rad,*]
if (n_elements(size(im)) eq 5) then mim = im $
  else mim = median(im, dimension=3)

mwrfits, mim, path+'PSF_'+fn+'.fits', newhdr

;=============================================================================================

stop
end
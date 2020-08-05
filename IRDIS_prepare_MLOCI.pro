@bandpass_filter_MOD.pro
@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@headfits.pro
@fxposit.pro
@mrd_hread.pro
@get_eso_keyword.pro
@sph_ird_transmission.pro
@strnumber.pro
@interpol.pro
@filter_image.pro
@factor.pro
@prime.pro
@psf_gaussian.pro
@gaussian.pro
@convolve.pro
@writefits.pro
@mkhdr.pro
@sxaddpar.pro
@get_date.pro
@daycnv.pro
@sxdelpar.pro
@sxpar.pro
@real_part.pro
@mrdfits.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@mfilter_am.pro

pro IRDIS_prepare_MLOCI

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

;========================================================================

;file selection needed to get header information regarding dit, NDs, BBs for FCP

readcol, '/home/amueller/work/IDLlibs/AO/SPHERE/datapaths.txt', tmp, format='a', /silent
print, ''
for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]

print, ''
read, 'Select Path: ', selp
path = tmp[selp-1]+'/RAW/IRDIS/'

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

idx = where(catg eq 'SCIENCE')
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
read, 'ONE Science Frame: ', sel
;sel = '6'
scf = tfile[sel-1]
read, 'ONE Object,Flux or PSF if none: ', sel
fluxf = tfile[sel-1]

;extract needed keywords

hdrsc = headfits(scf, exten=0, /silent)
  sopti2 = strcompress(get_eso_keyword(hdrsc,'HIERARCH ESO INS1 OPTI2 NAME'),/rem)
  scombind = strcompress(get_eso_keyword(hdrsc, 'HIERARCH ESO INS4 FILT2 NAME'),/rem)	;Assembly for infrared neutral density
  sfilt = strcompress(get_eso_keyword(hdrsc, 'HIERARCH ESO INS1 FILT NAME'),/rem)	;IRDIS filter unique ID
  dit_sci = double(get_eso_keyword(hdrsc, 'HIERARCH ESO DET SEQ1 DIT'))
hdrfl = headfits(fluxf, exten=0, /silent)
  fopti2 = strcompress(get_eso_keyword(hdrfl,'HIERARCH ESO INS1 OPTI2 NAME'),/rem)
  fcombind = strcompress(get_eso_keyword(hdrfl, 'HIERARCH ESO INS4 FILT2 NAME'),/rem)	;Assembly for infrared neutral density
  ffilt = strcompress(get_eso_keyword(hdrfl, 'HIERARCH ESO INS1 FILT NAME'),/rem)	;IRDIS filter unique ID
  dit_fl = double(get_eso_keyword(hdrfl, 'HIERARCH ESO DET SEQ1 DIT'))

;extract filter names for sph_ird_transmission.pro
db_sci = sopti2
bb_sci = sfilt
if (scombind eq 'OPEN') then nd_sci = '0.0' $
  else nd_sci = strmid(scombind, strpos(scombind, '_', /reverse_search)+1, 3)

db_fl = fopti2
bb_fl = ffilt
if (fcombind eq 'OPEN') then nd_fl = '0.0' $
  else nd_fl = strmid(fcombind, strpos(fcombind, '_', /reverse_search)+1, 3)

transm_sci = 10.^sph_ird_transmission(bb_sci, db_sci, nd_sci)
transm_fl = 10.^sph_ird_transmission(bb_fl, db_fl, nd_fl)

ditfact = dit_sci/dit_fl

if (bb_sci ne bb_fl and db_sci ne db_fl) then begin

  print, ''
  print, 'You screwed up with your data.'
  stop

endif


;========================================================================

stars = ['TWHya', 'betaPic', 'HIP60459', 'BlindTest', 'GJ1046', 'HD95086']
nstars = n_elements(stars)
print, ''
for i=0,nstars-1 do print, strcompress(i+1, /rem), ' ', stars[i]
sel = ''
read, 'Select Star?: ', sel

;========================================================================

;TWHya
if (sel eq '1') then begin
  datadir = '/home/amueller/work/SPHERE/data/TWHya/IRDIS/'
  tg_name = ['H2', 'H3']
  camera_filter = ['H2', 'H3']
  filepsf = ['PSF_H2.fits', 'PSF_H3.fits']	;des not exist yet
endif

;------------------------------------------------------------

;betaPic
if (sel eq '2') then begin
  datadir = '/home/amueller/work/SPHERE/data/betaPic/IRDIS/'
  tg_name = ['K1', 'K2']
  camera_filter = ['K1', 'K2']
  filepsf = ['PSF_K1.fits', 'PSF_K2.fits']	;des not exist yet
endif

;------------------------------------------------------------

;2nd Blind test, HIP60459
if (sel eq '3') then begin
  datadir = '/home/amueller/work/SPHERE/data/HIP60459/IRDIS/'
  tg_name = ['H2', 'H3']
  camera_filter = ['H2', 'H3']
  filepsf = ['PSF_H2.fits', 'PSF_H3.fits']
endif

;------------------------------------------------------------

;1st Blind test, star unknown
if (sel eq '4') then begin
  datadir = '/home/amueller/work/SPHERE/data/BlindTest/IRDIS/'
  tg_name = ['K1', 'K2']
  camera_filter = ['K1', 'K2']
  filepsf = ['PSF_K1.fits', 'PSF_K2.fits']
endif

;------------------------------------------------------------

;1st Blind test, star unknown
if (sel eq '5') then begin
  datadir = '/home/amueller/Downloads/GJ1046/'
  tg_name = ['J2', 'J3']
  camera_filter = ['J2', 'J3']
  filepsf = ['PSF_J2.fits', 'PSF_J3.fits']
endif

;------------------------------------------------------------

;1st Blind test, star unknown
if (sel eq '6') then begin
  datadir = '/home/amueller/work/SPHERE/data/HD95086/IRDIS/'
  tg_name = ['K1', 'K2']
  camera_filter = ['K1', 'K2']
  filepsf = ['PSF_K1.fits', 'PSF_K2.fits']
endif

;------------------------------------------------------------

;========================================================================

fquest = ''
read, 'Filter data? y/n: ', fquest
if (fquest eq 'y') then filterflag = '1' else filterflag = '0'


;========================================================================

plsc = 0.01225d0 ; plate scale [arcsec]
diam = 8.2 ; telescope diameter [m]

;========================================================================

procdir = datadir

for xx=0,n_elements(camera_filter)-1 do begin

  resdir = datadir+'MLOCI_'+camera_filter[xx]+'/'
  spawn, 'mkdir -p '+resdir

  tg_name_ori = tg_name[xx]
  procdir = datadir
  tg_name_dyn = tg_name_ori
  tg_name_bas = tg_name_ori


  ;H23 not the same as H32
  if (camera_filter[xx] eq 'K1' and xx eq 0) then lambda = 2.1025d-6
  if (camera_filter[xx] eq 'H2' and xx eq 0) then lambda = 1.5888d-6
  if (camera_filter[xx] eq 'H3' and xx eq 0) then lambda = 1.6653d-6
  if (camera_filter[xx] eq 'J2' and xx eq 0) then lambda = 1.1895d-6

  if (camera_filter[xx] eq 'K2' and xx eq 1) then lambda = 2.2550d-6
  if (camera_filter[xx] eq 'H2' and xx eq 1) then lambda = 1.5890d-6
  if (camera_filter[xx] eq 'H3' and xx eq 1) then lambda = 1.6671d-6
  if (camera_filter[xx] eq 'J3' and xx eq 1) then lambda = 1.2698d-6

  fwhm = lambda/diam*206265.d0/plsc ; fhwm of the psf in px, assuming no Lyot stop
  cutoff_l = 2.d0*fwhm ; parameter of high pass filter

  ;--------------------------------------------------------------------

  ;read in and scale PSF
  psf = mrdfits(datadir+filepsf[xx], 0, /silent)
  psf = psf*ditfact*transm_fl[xx]/transm_sci[xx]
  tg_name_dyn = 'dummy'
  tmp = psf
  if (filterflag eq '1') then begin

    psf = mfilter_am(tmp, highpass=2.d0*fwhm, lowpass=0.5d0*fwhm, /CA)
    spawn, 'rm '+datadir+'*dummy*.fits'

  endif

  psfpeak = (max(psf))[0]

  ;--------------------------------------------------------------------

  ;read in parallactic angle

  paral = mrdfits(datadir+'vec_'+camera_filter[xx]+'_paral.fits', 0, /silent)
  paral = -1.d0*paral

  ;--------------------------------------------------------------------

  ;read in data cube

  im = mrdfits(datadir+'img_'+tg_name[xx]+'_dc.fits', 0, /silent)
  tg_name_dyn = 'scidum'
  tmp = im
  if (filterflag eq '1') then begin

    filtim = mfilter_am(tmp, highpass=2.d0*fwhm, lowpass=0.5d0*fwhm, /CA)
    spawn, 'rm '+datadir+'*scidum*.fits'
    im = filtim

  endif

  ;--------------------------------------------------------------------
  ;output

  for i=0,n_elements(im[0,0,*])-1 do writefits, resdir+'img_'+camera_filter[xx]+'_frame'+'_'+strcompress(i+1,/rem)+'.fits', im[*,*,i]

  openw, lun, resdir+camera_filter[xx]+'_MLOCI_filelist.txt', width=2000, /get_lun

    for i=0,n_elements(im[0,0,*])-1 do printf, lun, 'img_'+camera_filter[xx]+'_frame'+'_'+strcompress(i+1,/rem)+'.fits', paral[i], psfpeak, format='(a22, f16.10, f25.10)'

  close, lun
  free_lun, lun

endfor


stop
end
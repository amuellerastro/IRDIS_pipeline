@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@headfits.pro
@fxposit.pro
@mrd_hread.pro
@get_eso_keyword.pro
@sigfig.pro

pro IRDIS_overview

;===============================================================================================

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
tech = strarr(nfiles)
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
  tech[i] = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO DPR TECH'),/rem)
  catg[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DPR CATG'),/rem)
  dit[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DET SEQ1 DIT'),/rem)
;   dit[i] = strmid(dit[i], 0, strpos(dit[i], '.'))
  dit[i] = sigfig(dit[i],2)
  icor[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS COMB ICOR'),/rem)		;coronograph
;   combind[i] = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO INS4 COMB IND'),/rem)	;ND filter
  combind[i] = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO INS4 FILT2 NAME'),/rem)	;ND filter
  opti2[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS1 OPTI2 NAME'),/rem)	;DB filter
;   filt2[i] = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO INS4 FILT2 NAME'),/rem)	;NIR neutral density filters
  filt[i] = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO INS1 FILT NAME'),/rem)	;BB filter

endfor


openw, lun, path+'files.txt', width=2000, /get_lun

  printf, lun, '                         File       Category             Type  DIT BB-Filt DB-Filt   ND-Filt  Coronograph                          Tech'
  for i=0,nfiles-1 do printf, lun, out[i], catg[i], type[i], dit[i], filt[i], opti2[i], combind[i], icor[i], tech[i], format='(a29,a15,a17,a5,a8,a8,a10,a13,a30)'

close, lun
free_lun, lun

stop
end
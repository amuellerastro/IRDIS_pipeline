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
@proceeding_text.pro
@mwrfits.pro
@fxaddpar.pro 
@sixlin.pro
@readcol.pro
@remchar.pro
@detabify.pro
@fxparpos.pro
@cghistoplot.pro
@convert_to_type.pro
@cgcheckforsymbols.pro
@cgscalevector.pro
@fpufix.pro
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
@cgresizeimage.pro
@congrid.pro
@writefits.pro
@mkhdr.pro
@sxaddpar.pro
@get_date.pro
@daycnv.pro
@sxdelpar.pro
@sxpar.pro

function func_linear, x, p

  fit = p[0]*x+p[1]

  return, fit

end

pro IRDIS_create_InstrumentFlat

;===============================================================================================

readcol, '/home/amueller/work/IDLlibs/AO/SPHERE/IRDIS_mask.txt', dx, dy, format='d,d', /silent

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

;selection of FF frames and BP

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
sel = intarr(2)
read, 'FFs (first and last number): ', sel
fff = tfile[indgen(sel[1]-sel[0]+1)+sel[0]-1]
sel = 0
;read, 'Bad Pixel Mask: ', sel
;bpf = tfile[sel-1]
bpf = path+'static_badpixels.fits'

;===============================================================================================

nx = 2048
ny = 1024
; xc = 512
; yc = 512
dim = 800

bp = mrdfits(bpf, 0, /silent)

nff = n_elements(fff)
ff = dblarr(nx,ny,nff)
dit = dblarr(nff)


;sort frames w.r.t. ascending exposure time
for i=0,nff-1 do begin

  hdr = headfits(fff[i], exten=0)
  dit[i] = double(get_eso_keyword(hdr, 'HIERARCH ESO DET SEQ1 DIT'))

endfor

idxsort = sort(dit)
fff = fff[idxsort]
dit = dit[idxsort]

;for i=0,nff-1 do ff[*,*,i] = mrdfits(fff[i], 0, /silent)
dum = mrdfits(fff[0],0,/silent)
nframes = n_elements(dum[0,0,*])
for i=0,nff-1 do begin

  if (nframes eq 1) then ff[*,*,i] = mrdfits(fff[i], 0, /silent)

  if (nframes gt 1) then begin

    tmp = mrdfits(fff[i], 0, /silent)
    ff[*,*,i] = median(tmp,dim=3)

  endif


endfor

;49Cet test
; for i=0,nff-1 do begin
; 
;   ff1 = mrdfits(fff[i], 0, /silent)
;   ff1 = median(ff1, dim=3)
; 
;   ff[*,*,i] = ff1
; 
; endfor


;===============================================================================================

;only consider 800x800pixel of left and right detector respectively

ffl = ff[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1,*]
ffr = ff[1024:2047,*,*]
ffr = ffr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1,*]

bpl = bp[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1,*]
bpr = bp[1024:2047,*,*]
bpr = bpr[dx[1]-dim/2:dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1,*]

;===============================================================================================

xa = dindgen(nff)
slopel = dblarr(dim, dim)
sloper = dblarr(dim, dim)

dumerr = dblarr(nff)
dumerr[*] = 1.

; window, 0
for xx=0,dim-1 do begin

  for yy=0,dim-1 do begin

    if (bpl[xx,yy] ne 1.) then begin

      ;meanff = ff[xx,yy,*]/dit
      ;start_params = [1., 10.]
      ;result = mpfitfun('func_linear', dit, ff[xx,yy,*], dumerr, start_params, yfit=yfit, /quiet)	;too slow
      ;slope[xx,yy] = result[0]
      sixlin, dit, ffl[xx,yy,*], a, siga, b, sigb
      slopel[xx,yy] = b[0]

;       plot, dit, meanff, psym=2, xr=[-1,15], xst=1
;       oplot, dit, yfit, color=cgcolor('red')

    endif

    if (bpr[xx,yy] ne 1.) then begin

      sixlin, dit, ffr[xx,yy,*], a, siga, b, sigb
      sloper[xx,yy] = b[0]

    endif

  endfor

  proceeding_text,loop=(dim), i=xx, prompt='> x direction   '+string(xx+1,form='(I4)')

endfor

;===============================================================================================

;normalisation

;median slope
idx = where(bpl eq 0.d0)
msll = median(slopel[idx], /even)
idx = where(bpr eq 0.d0)
mslr = median(sloper[idx], /even)

ff_finall = slopel/msll
ff_finalr = sloper/mslr

;===============================================================================================

;put sub frames back together

ff_final = dblarr(nx, ny)
ff_final[dx[0]-dim/2:dx[0]+dim/2-1, dy[0]-dim/2:dy[0]+dim/2-1] = ff_finall
ff_final[1024+dx[1]-dim/2:1024+dx[1]+dim/2-1, dy[1]-dim/2:dy[1]+dim/2-1] = ff_finalr

; mwrfits, ff_final, path+'instrument_flat.fits', /silent
writefits, path+'instrument_flat.fits', ff_final

;===============================================================================================

window, 0
!p.multi=[0,1,2]
idx = where(bpl eq 0.d0)
plot, slopel[idx], psym=3, /yn, xst=1
idx = where(bpr eq 0.d0)
plot, sloper[idx], psym=3, /yn, xst=1
!p.multi=[0,1,0]

window, 2, xs=1024, ys=512
cgimage, ff_final, stretch=2

window, 1, xs=500, ys=500
cghistoplot, ff_final, xr=[0.8,1.2], bin=0.003, yr=[0,2.d5]

;===============================================================================================

;===============================================================================================
;===============================================================================================
;===============================================================================================

; ;2nd method
; 
; xa = dindgen(nff)
; wmff = dblarr(nx,ny)
; 
; dumerr = dblarr(nff)
; dumerr[*] = 1.
; 
; for i=0,nff-1 do ff[*,*,i] = ff[*,*,i]/dit[i]
; 
; 
; ; window, 0
; for xx=0,nx-1 do begin
; 
;   for yy=0,ny-1 do begin
; 
;     if (bp[xx,yy] ne 1.) then begin
; 
;       wmff[xx,yy] = (total((sqrt(dit))*ff[xx,yy,*]))/total(sqrt(dit))
; 
;     endif
; 
;   endfor
; 
;   proceeding_text,loop=(nx), i=xx, prompt='> x direction   '+string(xx+1,form='(I4)')
; 
; endfor
; 
; ;===============================================================================================
; 
; ;normalisation and output
; 
; idx = where(bp eq 0.d0)
; ;median slope
; mff = median(wmff[idx], /even)
; 
; ff_final = wmff/mff
; 
; mwrfits, ff_final, path+'instrument_flat2.fits', /silent


stop
end

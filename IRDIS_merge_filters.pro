@rot.pro
@mwrfits.pro
@fxaddpar.pro
@mrdfits.pro
@fxposit.pro
@mrd_hread.pro
@fxpar.pro
@gettok.pro
@valid_num.pro
@mrd_skip.pro
@transmatch9.pro
@get_w.pro
@dist_circle.pro
@whereprams.pro
@nint.pro
@avg.pro
@amoeba.pro
@shift_sub.pro
@robust_sigma.pro
@dist.pro
@array_indices.pro
@proceeding_text.pro

pro IRDIS_merge_filters

;=====================================================================================

path = '/home/amueller/work/SPHERE/data/TWHya/IRDIS/'
f1 = 'img_H2_dc.fits'
f2 = 'img_H3_dc.fits'
l1 = 1.5888d-6
l2 = 1.6671d-6
v1 = 'vec_H2_paral.fits'
v2 = 'vec_H3_paral.fits'

imout = 'img_H2H3_dc.fits'
paout = 'vec_H2H3_paral.fits'

;=====================================================================================

;read in

im1 = mrdfits(path+f1, 0, hdr1, /silent)
im2 = mrdfits(path+f2, 0, hdr2, /silent)
pa1 = mrdfits(path+v1, 0, /silent)
pa2 = mrdfits(path+v2, 0, /silent)

;for testing
; im1=im1[*,*,0:2]
; im2=im2[*,*,0:2]
; pa1=pa1[0:2]
; pa2=pa2[0:2]


sz = size(im1)
dim = sz[1]

;=====================================================================================

;Radii for matching both filters
rin = [10,23,81]
rout = [22,80,dim/2-30]

;=====================================================================================

;adjust for lambda, always magnify

imcube = dblarr(sz[1], sz[2], 2.*sz[3])
pacube = dblarr(2.*sz[3])

nrad = n_elements(rin)
diff = dblarr(nrad,sz[1],sz[1],sz[3])

if (l1 lt l2) then begin

;   tmp = im1
; 
;   zoom = l2/l1
; 
;   for i=0,sz[3]-1 do tmp[*,*,i] = rot(im1[*,*,i], 0.0, zoom, sz[1]/2, sz[1]/2, c=-0.5, /pivot)
; 
;   imcube[*,*,0:sz[3]-1] = tmp
;   imcube[*,*,sz[3]:*] = im2
;   pacube[0:sz[3]-1] = pa1
;   pacube[sz[3]:*] = pa2

  for i=0,sz[3]-1 do begin

    for j=0,nrad-1 do begin

      diff[j,*,*,i] = transmatch9(im2[*,*,i],im1[*,*,i],mode=2,pr=2.,sr=0.1,ar=2.,zr=0.05, rin=rin[j],rout=rout[j], ftol=1.e-6,/nodisp)

    endfor

  endfor

endif

;=====================================================================================

if (l2 lt l1) then begin

;   tmp = im2
;   zoom = l1/l2
; 
;   for i=0,sz[3]-1 do tmp[*,*,i] = rot(im2[*,*,i], 0.0, zoom, sz[1]/2, sz[1]/2, c=-0.5, /pivot)
; 
;   imcube[*,*,0:sz[3]-1] = tmp
;   imcube[*,*,sz[3]:*] = im1
;   pacube[0:sz[3]-1] = pa2
;   pacube[sz[3]:*] = pa1

  for i=0,sz[3]-1 do begin

    for j=0,nrad-1 do begin

      diff[j,*,*,i] = transmatch9(im1[*,*,i],im2[*,*,i],mode=2,pr=2.,sr=0.1,ar=2.,zr=0.05, rin=rin[j],rout=rout[j], ftol=1.e-6,/nodisp)

    endfor

  endfor

endif

;=====================================================================================

;combine the different annuli in one image

im = reform(diff[0,*,*,*])	;dblarr(dim,dim,sz[3])
mask_t = shift(dist(dim),dim/2,dim/2)

for i=0,sz[3]-1 do begin

  for j=0,nrad-1 do begin

    mask = mask_t ge rin[j] and mask_t le rout[j]

    idx0 = where(mask gt 0.)

    ;wheretomulti, mask, idx0, col, row, frame
    idx = array_indices(mask, idx0)

    for k=0L,n_elements(idx[0,*])-1 do im[idx[0,k],idx[1,k],i] = diff[j,idx[0,k],idx[1,k],i]

  endfor

  proceeding_text,loop=sz[3], i=i, prompt='> Combining annuli for frame        '+string(i+1,form='(I4)')

endfor

if (l1 lt l2) then im1 = im2-diff
if (l2 lt l1) then im2 = im1-diff

imcube[*,*,0:sz[3]-1] = im1
imcube[*,*,sz[3]:*] = im2
pacube[0:sz[3]-1] = pa1
pacube[sz[3]:*] = pa2

;=====================================================================================

;free some memory
im1 = 0
im2 = 0
tmp = 0

;=====================================================================================

;write out cubes

mwrfits, imcube, path+imout, hdr1
mwrfits, pacube, path+paout

stop
end



; 'K1' lambda = 2.1025d-6
; 'H2' lambda = 1.5888d-6
; 'H3' lambda = 1.6653d-6
; 'J2' lambda = 1.1895d-6
; 
; 'K2' lambda = 2.2550d-6
; 'H2' lambda = 1.5890d-6
; 'H3' lambda = 1.6671d-6
; 'J3' lambda = 1.2698d-6
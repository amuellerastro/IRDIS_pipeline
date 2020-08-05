;in case of IRDIS DBI 2 images exist next to each other

pro IRDIS_extract_bands_BlindTest

path = '/home/amueller/work/SPHERE/data/BlindTest/IRDIS/orig/'

; file = 'center_im.fits'
; filter1 = 'K1'
; filter2 = 'K2'
; 
; ;-------------------------------------------------
; 
; st = mrdfits(path+file, 0, hdr, /silent)
; 
; im1 = st[*,*,*,0]
; im2 = st[*,*,*,1]
; 
; for i=0,n_elements(hdr)-1 do begin
; 
;   if (strmatch(hdr[i], '*NAXIS3*') eq 1) then newhdr = arrdelete(hdr, at=i, length=1)
;   if (strmatch(hdr[i], '*NAXIS4*') eq 1) then newhdr = arrdelete(hdr, at=i, length=1)
; 
; endfor
; 
; 
; mwrfits, path+im1, 'img_'+filter1+'_dc.fits', hdr, /silent
; mwrfits, path+im2, 'img_'+filter2+'_dc.fits', hdr, /silent



;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------

;extract median unsat PSF as a fake companion for PCA by Mawet
; 
file = 'median_unsat.fits'
filter1 = 'K1'
filter2 = 'K2'

;-------------------------------------------------

st = mrdfits(path+file, 0, hdr, /silent)

im1 = st[*,*,0]
im2 = st[*,*,1]

for i=0,n_elements(hdr)-1 do begin

  if (strmatch(hdr[i], '*NAXIS3*') eq 1) then newhdr = arrdelete(hdr, at=i, length=1)
  if (strmatch(hdr[i], '*NAXIS4*') eq 1) then newhdr = arrdelete(hdr, at=i, length=1)

endfor


; im1 = im1;/max(im1)
; im2 = im2;/max(im2)

mwrfits, im1, path+'PSF_'+filter1+'.fits', newhdr;, /silent
mwrfits, im2, path+'PSF_'+filter2+'.fits', newhdr;, /silent



stop
end
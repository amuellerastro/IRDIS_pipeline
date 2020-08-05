;assumption:
;DIT nd ND filter identical between nights

pro IRDIS_combine_nights

base = '/data/beegfs/astro-storage/groups/henning/amueller/SPHERE_Miriam/'
id = 'IPTau'

spawn, 'mkdir -p '+base+id
spawn, 'mkdir -p '+base+id+'/RAW'
spawn, 'mkdir -p '+base+id+'/RAW/IRDIS/'
spawn, 'mkdir -p '+base+id+'/IRDIS'
resdir = base+id+'/IRDIS/'
filter = 'K'

im11 = mrdfits(base+id+'_1/RAW/IRDIS/img_'+filter+'1_dc.fits',0,hdr11,/silent)
im12 = mrdfits(base+id+'_1/RAW/IRDIS/img_'+filter+'2_dc.fits',0,hdr12,/silent)
im21 = mrdfits(base+id+'_2/RAW/IRDIS/img_'+filter+'1_dc.fits',0,hdr21,/silent)
im22 = mrdfits(base+id+'_2/RAW/IRDIS/img_'+filter+'2_dc.fits',0,hdr22,/silent)

v1 = mrdfits(base+id+'_1/RAW/IRDIS/vec_'+filter+'1_paral.fits',0,/silent)
ha1 = mrdfits(base+id+'_1/RAW/IRDIS/vec_'+filter+'2_paral.fits',1,/silent)
v2 = mrdfits(base+id+'_2/RAW/IRDIS/vec_'+filter+'1_paral.fits',0,/silent)
ha2 = mrdfits(base+id+'_2/RAW/IRDIS/vec_'+filter+'2_paral.fits',1,/silent)

fpsf11 = file_search(base+id+'_1/RAW/IRDIS/PSF_'+filter+'1_*.fits', count=n11)
fpsf12 = file_search(base+id+'_1/RAW/IRDIS/PSF_'+filter+'2_*.fits', count=n12)
fpsf21 = file_search(base+id+'_2/RAW/IRDIS/PSF_'+filter+'1_*.fits', count=n21)
fpsf22 = file_search(base+id+'_2/RAW/IRDIS/PSF_'+filter+'2_*.fits', count=n22)

nx = n_elements(im11[*,0,0])
ny = n_elements(im11[0,*,0])
n1 = n_elements(im11[0,0,*])
n2 = n_elements(im21[0,0,*])

tmp = mrdfits(fpsf11[0], 0, /silent)
pnx = n_elements(tmp[*,0,0])
pny = n_elements(tmp[0,*,0])
; pn1 = n_elements(fpsf11[0,0,*])+n_elements(fpsf21[0,0,*])
; pn2 = n_elements(fpsf12[0,0,*])+n_elements(fpsf22[0,0,*])
psfhdr1 = headfits(base+id+'_1/RAW/IRDIS/PSF_'+filter+'1.fits', exten=0, /silent)
psfhdr2 = headfits(base+id+'_1/RAW/IRDIS/PSF_'+filter+'2.fits', exten=0, /silent)


im1 = dblarr(nx,ny,n1+n2)
im2 = dblarr(nx,ny,n1+n2)
v = dblarr(n1+n2)
ha = dblarr(n1+n2)
; psf1 = dblarr(pnx,pny,pn1)
; psf2 = dblarr(pnx,pny,pn2)

for i=0,n1-1 do im1[*,*,i] = im11[*,*,i]
for i=0,n2-1 do im1[*,*,n1+i] = im21[*,*,i]
for i=0,n1-1 do im2[*,*,i] = im12[*,*,i]
for i=0,n2-1 do im2[*,*,n1+i] = im22[*,*,i]

for i=0,n1-1 do v[i] = v1[i]
for i=0,n2-1 do v[n1+i] = v2[i]

for i=0,n1-1 do ha[i] = ha1[i]
for i=0,n2-1 do ha[n1+i] = ha2[i]

idxs = sort(ha)
im1 = im1[*,*,idxs]
im2 = im2[*,*,idxs]
v = v[idxs]
ha = ha[idxs]



;psf

;first filter
tmp1 = dblarr(pnx,pny,1000)
for i=0,n11-1 do begin

    t1 = mrdfits(fpsf11[i],0,/silent)
    t3 = total(tmp1,1)
    t4 = total(t3,1)
    idx1 = where(t4 eq 0.)
    for j=0,n_elements(t1[0,0,*])-1 do tmp1[*,*,idx1[0]+j] = t1[*,*,j]

endfor
tmp2 = dblarr(pnx,pny,1000)
for i=0,n21-1 do begin

    t2 = mrdfits(fpsf21[i],0,/silent)
    t5 = total(tmp2,1)
    t6 = total(t5,1)
    idx2 = where(t6 eq 0.)

    for j=0,n_elements(t2[0,0,*])-1 do tmp2[*,*,idx2[0]+j] = t2[*,*,j]

endfor

tmp = [[[tmp1]],[[tmp2]]]
t1 = total(tmp,1)
t2 = total(t1,1)
idx = where(t2 ne 0.)
psf1 = tmp[*,*,idx]

;second filter
tmp1 = dblarr(pnx,pny,1000)
for i=0,n12-1 do begin

    t1 = mrdfits(fpsf12[i],0,/silent)
    t3 = total(tmp1,1)
    t4 = total(t3,1)
    idx1 = where(t4 eq 0.)
    for j=0,n_elements(t1[0,0,*])-1 do tmp1[*,*,idx1[0]+j] = t1[*,*,j]

endfor
tmp2 = dblarr(pnx,pny,1000)
for i=0,n22-1 do begin

    t2 = mrdfits(fpsf22[i],0,/silent)
    t5 = total(tmp2,1)
    t6 = total(t5,1)
    idx2 = where(t6 eq 0.)

    for j=0,n_elements(t2[0,0,*])-1 do tmp2[*,*,idx2[0]+j] = t2[*,*,j]

endfor

tmp = [[[tmp1]],[[tmp2]]]
t1 = total(tmp,1)
t2 = total(t1,1)
idx = where(t2 ne 0.)
psf2 = tmp[*,*,idx]


;output
writefits, resdir+'img_'+filter+'1_dc.fits', im1, hdr11
writefits, resdir+'img_'+filter+'2_dc.fits', im2, hdr12
writefits, resdir+'vec_'+filter+'1_paral.fits', v
writefits, resdir+'vec_'+filter+'1_paral.fits', ha, /append
writefits, resdir+'vec_'+filter+'2_paral.fits', v
writefits, resdir+'vec_'+filter+'2_paral.fits', ha, /append
writefits, base+id+'/RAW/IRDIS/'+'PSF_'+filter+'1_1.fits', psf1, psfhdr1
writefits, base+id+'/RAW/IRDIS/'+'PSF_'+filter+'2_1.fits', psf2, psfhdr2
writefits, resdir+'PSF_'+filter+'1.fits', median(psf1, dim=3), psfhdr1
writefits, resdir+'PSF_'+filter+'2.fits', median(psf2, dim=3), psfhdr2


stop
end

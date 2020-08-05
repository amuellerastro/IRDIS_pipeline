@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@mrdfits.pro
@fxposit.pro
@mrd_hread.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@mwrfits.pro
@fxaddpar.pro
@arrdelete.pro
@writefits.pro
@mkhdr.pro
@sxaddpar.pro
@get_date.pro
@daycnv.pro
@sxdelpar.pro
@sxpar.pro


pro IRDIS_remove_frame

;================================================================================

readcol, '/home/amueller/work/IDLlibs/AO/SPHERE/datapaths.txt', tmp, format='a', /silent
print, ''
for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]

print, ''
read, 'Select Path: ', selp
path = tmp[selp-1]+'/RAW/IRDIS/'

;================================================================================

; read, 'Number of bad frames: ', nbad
; bad = intarr(nbad)
; read, 'Enter number of bad frames (index starts with 1!): ', bad

; ;TCha
; bad = [6,7,8,13,15,16,21,22,23,24,25,26,27,28,29,32,33,34,35,36,37,38,41,42,43,45,46,47,48,49,50,64,65,66,67,76,77,78,79,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96]
; ; PROVIDE BAD FRAMES
; ; bad = [50]

;V1032Cen
; bad = [21,22,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64]

;PDS70_H_201505
;bad = [46,48,49,50,51,52,53,58,59,60]
;PDS70_H_201506
bad = [51,52,53,54,55,56,57,58,59,60,61,62,63,64]


nbad = n_elements(bad)

bad = bad-1	;bring the numbers back to IDL counting (1st index = 0)

; nbad = bad[1]-bad[0]+1
; nbad = nbad[0]

;================================================================================

;extract file names

file = file_search(path+'img_*_dc.fits', count=nfiles)
if (nfiles ne 2) then stop
pos1 = strpos(file[0], '/', /reverse_search)
file1 = strmid(file[0], pos1+1, strlen(file[0])-pos1-1)
pos2 = strpos(file[1], '/', /reverse_search)
file2 = strmid(file[1], pos2+1, strlen(file[1])-pos2-1)

file = file_search(path+'vec_*_paral.fits', count=nfiles)
if (nfiles ne 2) then stop
pos1 = strpos(file[0], '/', /reverse_search)
afile1 = strmid(file[0], pos1+1, strlen(file[0])-pos1-1)
pos2 = strpos(file[1], '/', /reverse_search)
afile2 = strmid(file[1], pos2+1, strlen(file[1])-pos2-1)

;================================================================================

spawn, 'mv '+path+file1+' '+path+file1+'.orig'
spawn, 'mv '+path+file2+' '+path+file2+'.orig'

spawn, 'mv '+path+afile1+' '+path+afile1+'.orig'
spawn, 'mv '+path+afile2+' '+path+afile2+'.orig'

cpim1 = file1+'.orig'
cpim2 = file2+'.orig'

acpim1 = afile1+'.orig'
acpim2 = afile2+'.orig'

;================================================================================

im1 = mrdfits(path+cpim1, 0, hdr1, /silent)
im2 = mrdfits(path+cpim2, 0, hdr2, /silent)

a1 = mrdfits(path+acpim1, 0, hdra1, /silent)
a2 = mrdfits(path+acpim2, 0, hdra2, /silent)
ha1 = mrdfits(path+acpim1, 1, /silent)
ha2 = mrdfits(path+acpim2, 1, /silent)


sz = size(im1)

;================================================================================

good = indgen(sz[3])
for i=0,nbad-1 do good = arrdelete(good, at=bad[i]-i, length=1, /overwrite)

new1 = im1[*,*,good]
new2 = im2[*,*,good]

newa1 = a1[good]
newa2 = a2[good]
newha1 = ha1[good]
newha2 = ha2[good]


; stop

; new1 = dblarr(sz[1], sz[2], sz[3]-nbad)
; new2 = dblarr(sz[1], sz[2], sz[3]-nbad)
; 
; ;assuming 1st one is OK
; if (bad[0] gt 0 and bad[1] lt sz[3]-1) then begin
; 
;   new1[*,*,0:bad[0]-1] = im1[*,*,0:bad[0]-1]
;   new1[*,*,bad[0]:*] = im1[*,*,bad[1]+1:*]
; 
;   new2[*,*,0:bad[0]-1] = im2[*,*,0:bad[0]-1]
;   new2[*,*,bad[0]:*] = im2[*,*,bad[1]+1:*]
; 
; endif else begin
; 
;   print, ''
;   print, 'Code does currently not allow the first and last frame to be bad.'
; 
; endelse

; mwrfits, new1, path+file1, hdr1
; mwrfits, new2, path+file2, hdr2
writefits, path+file1, new1, hdr1
writefits, path+file2, new2, hdr2

;================================================================================

; newa1[0:bad[0]-1] = im1[0:bad[0]-1]
; newa1[bad[0]:*] = im1[bad[1]+1:*]
; 
; newa2[0:bad[0]-1] = im2[0:bad[0]-1]
; newa2[bad[0]:*] = im2[bad[1]+1:*]
; 

; newa1 = arrdelete(a1, at=bad, length=1)
; newa2 = arrdelete(a2, at=bad, length=1)

; mwrfits, newa1, path+afile1, hdra1, /silent
; mwrfits, newa2, path+afile2, hdra2, /silent
writefits, path+afile1, newa1, hdra1
writefits, path+afile1, newha1, /append
writefits, path+afile2, newa2, hdra2
writefits, path+afile2, newha2, /append

;================================================================================

im1 = 0
im2 = 0
new1 = 0
new2 = 0

stop
end
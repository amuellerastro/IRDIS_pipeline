pro IRDIS_prepare_data_RDI_amueller, extract=extract

filter = ['H2', 'H3']
cut = 110    ;final cube 50x50px

read, 'Average X frames (0 = no avergae): ', xf

if keyword_set(extract) then begin

    quest = ''
    ;read, 'Average frames? y/n: ', quest
    ;if (quest eq 'y') then 

    base = '/data/beegfs/astro-storage/groups/feldt/sphere_shared/automatic_reduction/esorexer/complete_reduction/'
    resdir1 = '/data/beegfs/astro-storage/groups/launhardt/amueller/IRDIS_RDI/data/'+filter[0]+'_'+strcompress(fix(xf),/rem)+'ave/'
    resdir2 = '/data/beegfs/astro-storage/groups/launhardt/amueller/IRDIS_RDI/data/'+filter[1]+'_'+strcompress(fix(xf),/rem)+'ave/'
    
    spawn, 'mkdir -p '+resdir1
    spawn, 'mkdir -p '+resdir2
    
    datafile = '/data/beegfs/astro-storage/groups/launhardt/amueller/IRDIS_RDI/stars.txt'
    readcol, datafile, datapath, format='a', /silent
    npath = n_elements(datapath)
    
    for xx=0,npath-1 do begin

        file = file_search(base+datapath[xx]+'/DB_H23/*/converted_mastercubes/center_im.fits', count=nc)
        pfile = file_search(base+datapath[xx]+'/DB_H23/*/converted_mastercubes/rotnth.fits', count=np)
        ffile = file_search(base+datapath[xx]+'/DB_H23/*/converted_mastercubes/median_unsat.fits', count=nf)
    
;         paral = dblarr(nc)
;         mjd = dblarr(nc)
;         airmass = dblarr(nc)
;         seeing = dblarr(nc)
;         r0 = dblarr(nc)
;         t0 = dblarr(nc)
;         WindSpeed = dblarr(nc)
;         Humidity = dblarr(nc)
        
;         hdr0 = headfits(file[0], exten=0, /silent)
;         naxis = get_eso_keyword(hdr0, 'NAXIS1')
;         date = get_eso_keyword(hdr0, 'DATE-OBS')
;         night = strmid(date, 0, 10)
;         target = strcompress(get_eso_keyword(hdr0, 'HIERARCH ESO OBS TARG NAME'),/rem)
;         cube = dblarr(naxis,naxis,nc)
        
        for i=0,nc-1 do begin
        
            cube = mrdfits(file[i],0,hdr,/silent)
            paral = mrdfits(pfile[i],0,/silent)
            paral = -1.*paral   ;need to do this because of SPHERE DC angles
            
            flux = mrdfits(ffile[i],0,fhdr,/silent)
            fdim = double(get_eso_keyword(fhdr, 'NAXIS1'))
            flux = flux[fdim/2-13:fdim/2+13,fdim/2-13:fdim/2+13,*,0]
            ;here I want to take the median but some files have are missing 2 filters 2exposures. So I just take the first image of each filter, assuming these are 2 filters
            flux1 = (flux[*,*,0])
            flux2 = (flux[*,*,1])
            
            naxis = double(get_eso_keyword(hdr, 'NAXIS1'))
            naxis3 = double(get_eso_keyword(hdr, 'NAXIS3'))
            target = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO OBS TARG NAME'),/rem)
            date = get_eso_keyword(hdr, 'DATE-OBS')
            night = strmid(date, 0, 10)
            
            cube1 = cube[naxis/2-cut/2:naxis/2+cut/2-1,naxis/2-cut/2:naxis/2+cut/2-1,*,0]
            cube2 = cube[naxis/2-cut/2:naxis/2+cut/2-1,naxis/2-cut/2:naxis/2+cut/2-1,*,1]
            
            if (quest eq 'y') then begin
            
                if (naxis3 ge 2.*xf) then begin
            
                    rest = naxis3 mod xf
                    tnf = naxis3-rest
                    idxrf = round(linspace(xf-1,tnf-1,tnf/xf))

                    tc1 = dblarr(n_elements(cube1[*,0,0]),n_elements(cube1[*,0,0]),n_elements(idxrf))
                    tc2 = dblarr(n_elements(cube2[*,0,0]),n_elements(cube2[*,0,0]),n_elements(idxrf))
                    tp = dblarr(n_elements(idxrf))
                    
                    for j=0,n_elements(idxrf)-1 do begin
                    
                        if (j eq 0) then tc1[*,*,j] = median(cube1[*,*,0:idxrf[j]], dim=3) else tc1[*,*,j] = median(cube1[*,*,idxrf[j-1]+1:idxrf[j]], dim=3)  
                        
                        if (j eq 0) then tc2[*,*,j] = median(cube2[*,*,0:idxrf[j]], dim=3) else tc2[*,*,j] = median(cube2[*,*,idxrf[j-1]+1:idxrf[j]], dim=3)  
                        
                        if (j eq 0) then tp[j] = mean(paral[0:idxrf[j]]) else tp[j] = mean(paral[idxrf[j-1]+1:idxrf[j]])
                    
                    endfor
                    
                    cube1 = tc1
                    cube2 = tc2
                    paral = tp

                endif
                    
            endif
            
            
            fn1 = resdir1+target+'_'+night+'_'+filter[0]+'.fits'
            writefits, fn1, cube1, hdr
            writefits, fn1, paral, /append
            writefits, fn1, flux1, /append
            
            fn2 = resdir2+target+'_'+night+'_'+filter[1]+'.fits'
            writefits, fn2, cube2, hdr
            writefits, fn2, paral, /append
            writefits, fn2, flux2, /append


;             cube[*,*,i] = mrdfits(file[i],0,hdr,/silent)
;             paral[i] = mrdfits(pfile[i],0,/silent)
;             mjd[i] = get_eso_keyword(hdr,'MJD-OBS')
;             t1seeing = double(get_eso_keyword(hdr,'HIERARCH ESO TEL AMBI FWHM END'))
;             t2seeing = double(get_eso_keyword(hdr,'HIERARCH ESO TEL AMBI FWHM START'))
;             seeing[i] = mean([t1seeing,t2seeing])
;             airmass[i] = double(get_eso_keyword(hdr, 'AIRMASS'))
;             r0[i] = double(get_eso_keyword(hdr,'HIERARCH ESO AOS RTC DET DST R0MEAN'))
;             Humidity[i] = double(get_eso_keyword(hdr,'HIERARCH ESO TEL AMBI RHUM'))
;             t0[i] = double(get_eso_keyword(hdr,'HIERARCH ESO AOS RTC DET DST T0MEAN'))
;             WindSpeed[i] = double(get_eso_keyword(hdr,'HIERARCH ESO TEL AMBI WINDSP'))
            
        endfor
        
;         val = dblarr(nc)
;         medim = median(cube, dim=3)
;         for i=0,nc-1 do val[i] = ms_ssim(medim, cube[*,*,i])
;         resistant_mean, val, 3, t1, t2, nbad, /double, goodvec=idxg, badvec=idxb
;         
;         cube_cut = cube[naxis/2-cut/2:naxis/2+cut/2-1,naxis/2-cut/2:naxis/2+cut/2-1,idxg]
;         paral = paral[idxg]
;         mjd = mjd[idxg]
;         airmass = airmass[idxg]
;         seeing = seeing[idxg]
;         r0 = r0[idxg]
;         t0 = t0[idxg]
;         WindSpeed = WindSpeed[idxg]
;         Humidity = Humidity[idxg]
;     
; ;         sequenceinfo = {mjd:mjd,airmass:airmass,seeing:seeing,r0:r0,t0:t0,WindSpeed:WindSpeed,Humidity:Humidity}
; 
;         sequenceinfo = [[mjd],[airmass],[seeing],[r0],[t0],[WindSpeed],[Humidity]]
;         
;         fn = resdir+target+'_'+night+'.fits'
;         writefits, fn, cube_cut, hdr0
;         writefits, fn, paral, /append
;         writefits, fn, transpose(sequenceinfo), /append
    
        proceeding_text, loop=npath, i=xx, prompt='> Targets        '+string(xx+1,form='(I4)')
  
    endfor

endif

dim = cut
xa = dindgen(dim*dim)
errdum = dblarr(dim^2.)
errdum[*] = 1.


;-------------------------------------------------------------------------

for ifilt=0,n_elements(filter)-1 do begin

    datapath = '/data/beegfs/astro-storage/groups/launhardt/amueller/IRDIS_RDI/data/'+filter[ifilt]+'_'+strcompress(fix(xf),/rem)+'ave/'
    file = file_search(datapath+'*fits', count=nfiles)

    for xx=0,nfiles-1 do begin


        im = mrdfits(file[xx],0,hdr,/silent)
        paral = mrdfits(file[xx],1,/silent)
;         data = mrdfits(file[xx],2,/silent)
        
        target = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO OBS TARG NAME'),/rem)
        naxis3 = double(get_eso_keyword(hdr, 'NAXIS3'))

        ;file name for output
        pos1 = strpos(file[xx], '/', /reverse_search)
        pos2 = strpos(file[xx], '.fits', /reverse_search)
        fn = strmid(file[xx], pos1+1, pos2-pos1-1)

;         if (naxis3 ne (size(data))[2]) then stop

        flux = dblarr(naxis3)
        sdevflux = dblarr(naxis3)
        meanflux = dblarr(naxis3)
        
   		mask_t = shift(dist(dim), dim/2., dim/2.)
   		mask = mask_t ge 4.; and mask_t le 15.
   		idxmask = where(mask ne 0.)
   		npxmask = double(n_elements(idxmask))
        ;mask = 1.
        ;npxmask = dim*dim
        ;idxmask = indgen(dim*dim)
            
        im1d = dblarr(naxis3, npxmask)
            
        for i=0,naxis3-1 do begin
        
            t1 = (im[*,*,i])[*]
            
            im1d[i,*] = t1[idxmask]
                
            flux[i] = total(im1d[i,*])
            sdevflux[i] = stddev(im1d[i,*])
            meanflux[i] = mean(im1d[i,*])

        endfor
        
        fn = datapath+fn+'.sav'
        save, target, im1d, naxis3, paral, flux, sdevflux, meanflux, mask, npxmask, idxmask, filename=fn
        
        proceeding_text, loop=nfiles, i=xx, prompt='> Target        '+string(xx+1,form='(I4)')

    endfor
    
endfor

stop
end

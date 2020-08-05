IRDIS
=====

Requirements: 
-------------

  -coyote library
  -astrolib
  -library from Craig Markwardt (mpfit etc.)

  -update datapaths.txt with the main path, structure must be ...../<star>/IRDIS and ..../<star>/RAW/IRDIS
  -update paths for IRDIS_mask.txt and datapaths.txt in every IRDIS_<xxx>.pro file
  -IRDIS_reduce_DBI: provide star name with coordinates, plx, PM, RV

  for PCA
  Science frames must have an even number of rows/columns e.g. 800x800
  Star must be centered on dim/2, e.g. 400,400
  PSF must have a odd dimension, e.g. 23x23 and centered on 11,11 ((dim-1)/2)

  Files:
    3 dark,background with same DIT for object, sky, distortion
    5 flats with different DIT
    object,center
    object,fluc
    sky
    science object
    lamp,distort


Reduction:
----------

  IRDIS_overview - optional
    -lists properties of all frames

  compute_field_rotation - optional
    -computes and displays the amount of field rotation during 1 DIT as a function of angular distance

  IRDIS_create_Dark
    -create also a dark for Object,Flux

  IRDIS_create_BPmap
    -creation of BP mask based on FFs
    -two FFs with low and high DIT should be provided, if difference in DIT is too small the ratio is small and we might miss bad dark pixels
    -another set of BP is defined through the amount of flux present
    -for a considered dimension of 800x800 pixel a typical value of the percentage of BPs are 42% to 43%

  IRDIS_create_InstrumentFlat

  IRDIS_reduce_SkyBG

  OBSOLETE?
  simple_check_which_starcenter
    -in case of several star center frames check which one fits better
    -super simple, probably wrong anyway ;)

  IRDIS_distortion - !!!  DO NOT USE AT THE MOMENT !!!
    -assumptions:
      -holes in pinhole mask are evenly spaced
      -field distortions are smallest in the center of the image

  IRDIS_star_center
    - USE ONLY 1st STAR CENTER AFTER FLUX,CALIB AT THE BEGINNING.

  IRDIS_reduce_DBI
    coordinates of target and target name
    check if star is centered around 400,400px

  IRDIS_extract_reference_PSF /  not yet updated: IRDIS_extract_reference_PSF_ScienceFrame
    -reduces only one cube at a time (rename file manually), one could interpolate between several observations
    -PSF must have odd number of elements (e.g. 51x51)
    -can have a different DIT than science, careful with master dark
    -IRDIS_extract_reference_PSF_ScienceFrame: for data recorded without coronograph simply extract PSF from the already reduced science data

  IRDIS_extract_bands, qc=qc
    -creates cube and parallactic angle file
    -removes bad frams based on some criteria (mainly halo statistics)
      -creates .orig (all frames), .bad (bad frames only), and .fits (good frames)
      -check .bad for supposedly good frames and .fits for additional bad frames
      -makes problems in case of instable AO and clouds

  IRDIS_remove_frame - optional
    -removes bad frames after checking manually
    -bad frames have to be provided inside the code. frame 1 is 1 for index.

  OBSOLETE(?)
  IRDIS_merge_filters - optional
    -!! check frames for quality !!
    -manually provide file names and wavelength
    -combine 2 filters into one cube if desired
    -can take a while depending on size of input frames

  COPY IMAGE CUBE, ANGLES, PSF TO WORKING DIR (..../star/IRDIS/)

  SDI_SPHERE
  PCA_LOCI_ADI, /sphere
    -PCA, LLSG, LOCI, ADI
      !!! NOT set up for H32 yet !!!
    -PSF must have odd number of elements (e.g. 51x51)
    -path and file name has to be provided in code
    -currently eigen values are overwritten when doing fake planet injection in PCA

  planet_mass_detection_limit
    -??? filter handling correct ???
    -currently rounding the age, i.e. 2.3 Myr is 2Myr, 2.6Myr is 3 Myr in the model
    -currently only for PCA reduction

  OBSOLETE - make_overview_plot

  IRDIS_SADI
    -dimension must be even (e.g. 800x800)
    -im1-im2 after matching, then ADI
    -probably a lot of self subtraction for close-in objects

  IRDIS_prepare_MLOCI
    -bring data to the required format for MLOCI


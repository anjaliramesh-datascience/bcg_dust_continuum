# Imaging and Image Properties of ALMA BCGs

- params.py contains the choice of TCLEAN parameters
- These scripts are to be run on CASA with execfile(file)

## alma_imaging.py

## Natural Weighting

1. Auto-multithresh algorithm
    - Run in interactive mode
    - Alter mask if required
    - Save mask
2. Convert stored pbcor file to fits file

## Briggs Weighting

1. Automatic imaginf
2. Use mask generated with natural weighting
3. Convert to fits

## UVtapering

1. Tapered imaging of BCGs


## alma_flux.py

1. IMFIT on all the images generated
2. Plot flux vs beam size
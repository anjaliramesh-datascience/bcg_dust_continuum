# Compute Sersic profiles of BCGs

## read_fits_cutout(file):

1. Reads the fits image and header of the bcg
2. Generate 10"x10" cutout
3. Generate PSF image
4. Generate weightmap - rms map image (same size as the input cutout image)

## compute_sersic(image, weightmap, psf):

1. Runs statmorph on the images and returns morphology parameters

## plot_profiles()

1. Takes in the results generated from statmorph to generate profile plots
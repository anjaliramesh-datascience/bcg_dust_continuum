import sys
sys.path.append("/Users/arames52/bcg_dust_continuum/src/")
import bcg_parameter_file as p
import helper_functions as hf
import glob
import pandas as pd

# Create folders for output storage
folder_names_list = ["natural_imaging", "briggs_imaging", "tapered_imaging"]
hf.folder_creation(folder_names_list)
# get the list of weights and tapering requirements
weights, tapers = hf.get_weights_and_tapers()


def alma_imaging(file:str):
    # extracting object id to for image storage purposes
    bcg_name = file.split("/")[-1].split("_calibrated.ms")[0]
    # starting the loop to start the imaging of ALMA source with different image parameters like image size, cell size, 
    # deconvolving method, object point spread function, noise distribution, pixel weighting scheme, etc.
    for i in range(len(weights)):
        if tapers[i] != None:
            image_name = p.data_path + folder_names_list[2] + bcg_name + "_" + tapers[i]
        else:
            image_name = p.data_path + weights[i] + "_imaging/" + bcg_name
        # Retrieving saved mask
        mask = glob.glob(p.mask_path + bcg_name + "*.mask")[0]

        # Calling CASA's TCLEAN function with all the defined parameters
        tclean(vis = file, imagename = image_name, imsize = p.imsize, intent = p.intent, cell = p.cell,
        specmode = p.specmode, gridder = p.gridder, deconvolver = p.deconvolver, usemask = p.usemask, mask = mask,
         weighting = weights[i], niter = p.niter, field = p.field, pbcor = p.pbcor, robust = p.robust, 
         threshold = p.threshold, interactive = p.interactive, uvtaper = tapers[i])

        # Exporting .pbcor images to .fits image format
        exportfits(imagename = image_name + ".image.pbcor", fitsimage = image_name + ".fits")

    return None


for file in sorted(glob.glob(p.data_path + "Measurement_Sets/*.ms")):
    alma_imaging(file)


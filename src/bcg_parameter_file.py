# Path 
import sys
import os
code_path = "/Users/arames52/bcg_dust_continuum/src/"
data_path = '/Users/arames52/bcg_dust_continuum/notebook/data/'
measurement_sets_path = data_path + "Measurement_Sets/"
alma_fits_path = data_path + "ALMA_images/"
catalogs_path = data_path + 'catalogs/'
mask_path = data_path + 'Masks/'
natural_imaging_path = data_path + "natural_imaging/"
briggs_imaging_path = data_path + "briggs_imaging/"
tapered_imaging_path = data_path + "tapered_imaging/"
sys.path.append(code_path)
region_parameter_file_path = "/Users/arames52/bcg_dust_continuum/notebook/data/bcg_regions.txt"
swire_outpath = catalogs_path + "SWIRE_catalog/"
hermes_outpath = catalogs_path + "herschel_catalog/"

# TCLEAN
 
imsize = [864,864]
cell = '0.045arcsec'
gridder = 'standard'
deconvolver = 'hogbom'
specmode = 'mfs'
niter = 5000
pbcor = True
usemask = 'user'
field = '0'
threshold = '0.06mJy'
interactive = False
robust = 0.5
intent = 'OBSERVE_TARGET#ON_SOURCE'

# auto-masking parameters
sidelobethreshold= 3.0
noisethreshold= 5.0
lownoisethreshold= 1.5
negativethreshold= 0.0
minbeamfrac= 0.3



# data_ingestion.py

# Creating output directories
if not os.path.exists(swire_outpath):
    os.makedirs(swire_outpath)
if not os.path.exists(hermes_outpath):
    os.makedirs(hermes_outpath)



hermes_urls = ["https://hedam.lam.fr/HerMES/data/DR3/packages/L6-XMM-LSS-SWIRE_xID24_DR3.tar.bz2",
"https://hedam.lam.fr/HerMES/data/DR3/packages/L6-ELAIS-S1-SWIRE_xID24_DR3.tar.bz2",
"https://hedam.lam.fr/HerMES/data/DR3/packages/L5-CDFS-SWIRE_xID24_DR3.tar.bz2"]


# bcg categories

non_detections = ['ES1-12', 'ES1-26', 'XMM-19' , 'XMM-27']
weak_detections = ['CDFS19', 'ES1_z_0.88', 'ES1-35', 'XMM_z_0.81']
lensed_detections = ['XMM-11', 'ES1-34']
good_detections = ["CDFS-18", "ES1-18", "ES1-25", "ES1_z_0.99","ES1_z_0.99b","ES1_z_1.04","ES1_z_1.38","ES1_z_1.40",
"ES1_z_1.60", "ES1_z_1.65", "ES1_z_1.70", "XMM-113", "XMM-29", "XMM-30", "XMM_z_0.9", "XMM_z_1.0"]
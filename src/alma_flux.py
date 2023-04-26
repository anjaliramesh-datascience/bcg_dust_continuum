import sys
sys.path.append("/Users/arames52/Research/Analysis/")
import open_fits
import os
import glob
import photutils
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.io import fits
import scipy.ndimage as ndi
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import pickle
from photutils.aperture import SkyCircularAperture, CircularAperture, CircularAnnulus, SkyCircularAnnulus
from photutils.aperture import aperture_photometry
from photutils.aperture import ApertureStats

# CASA Images path

briggs_path = "/Users/arames52/Research/Data/CASA_Outputs/Briggs/fits/"
natural_path = "/Users/arames52/Research/Data/CASA_Outputs/Natural/fits/"
taper_07_path = "/Users/arames52/Research/Data/CASA_Outputs/UVtaper/0.7arcsec/fits/"
taper_1_path = "/Users/arames52/Research/Data/CASA_Outputs/UVtaper/1arcsec/fits/"
taper_2_path = "/Users/arames52/Research/Data/CASA_Outputs/UVtaper/2arcsec/fits/"

# Output path
output_path = "/Users/arames52/Research/Analysis/IMFIT/Residuals/"

bcg_detections = ["CDFS-18", "ES1-18","ES1-25", "ES1_z_0.99","ES1_z_0.99b","ES1_z_1.04","ES1_z_1.38","ES1_z_1.40",
"ES1_z_1.60", "ES1_z_1.65", "ES1_z_1.70", "XMM-113", "XMM-11", "XMM-29", "XMM-30", "XMM_z_0.9", "XMM_z_1.0"]
bcg_coords = pd.read_csv("/Users/arames52/Research/Data/BCG_coords.txt", names = ['bcg', 'ra','dec'], delim_whitespace=True)
path_dict = {0: briggs_path, 1: natural_path, 2: taper_07_path, 3: taper_1_path, 4: taper_2_path}
region_radius = {0:2, 1:2, 2:3, 3:4, 4:5}
img_type = {0:{}, 1:{}, 2: {}, 3: {}, 4: {}}
output_file_name = {0: "briggs_imfit", 1: "natural_imfit", 2: "taper_07_imfit", 3: "taper_1_imfit",
4: "taper_2_imfit"}

def alma_flux(image, region, residual_filename):
    imfit_result = imfit(image, region = region, residual = residual_filename)
    return imfit_result

def casa_imfit(bcg):

    ra = bcg_coords[bcg_coords['bcg'] == bcg]['ra'].values[0]
    dec = bcg_coords[bcg_coords['bcg'] == bcg]['dec'].values[0]
    position = SkyCoord(ra, dec, unit = 'deg', frame = 'fk5')
    imfit_details = {}

    for i in range(5):
        file = glob.glob(path_dict[i] + bcg + "_*.fits")[0]
        d,h,w = open_fits.read_alma_fits(file)
        pix_position = position.to_pixel(w)
        x, y = str(ra) + "deg", str(dec) + "deg"
        r = str(region_radius[i]) + "arcsec"
        imfit_region = "circle[["+ x +","+y+ "]," + r+ "]"
        residual_filename = output_path + str(i) + "/" + bcg + "_res.image"
        imfit_dict = alma_flux(file, imfit_region, residual_filename)
        img_type[i][bcg] = imfit_dict
        exportfits(imagename = residual_filename, fitsimage = residual_filename + '.fits')


for bcg in bcg_detections:
    casa_imfit(bcg)

for i in range(5):
    with open("/Users/arames52/Research/Analysis/IMFIT/" + output_file_name[i] + ".pkl", "wb") as f:
        pickle.dump(img_type[i], f)
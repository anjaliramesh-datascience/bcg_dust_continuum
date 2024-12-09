import glob
import pandas as pd
import pdb
import pickle
import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import ast

data_path = "/Users/arames52/Research/Data/Images/ALMA/"
residual_path = "/Users/arames52/Research/Data/GnS_profile/Gaussian/Residuals/"
imfit_results = {}
images = sorted(glob.glob(data_path + "*.fits"))
bcg_coords_df = pd.read_csv("/Users/arames52/Research/Data/BCG_coords.txt", delim_whitespace = True,
names = ['bcg', 'ra', 'dec'])
bcg_coords = {}
for ind, row in bcg_coords_df.iterrows():
    bcg_coords[row['bcg']] = [row['ra'], row['dec']]

for file in images:
    bcg_name = file.split("_natural.fits")[0].split(data_path)[-1]
    print(bcg_name)
    residual_file = residual_path + bcg_name + "_residual.image"
    ra = str(bcg_coords[bcg_name][0]) + "deg"
    dec = str(bcg_coords[bcg_name][1]) + "deg"
    if bcg_name == "XMM-113":
        width = "1.7arcsec"
    else:
        width = "3arcsec"
    region = "centerbox[[" + ra + "," + dec + "]," + "["+ width +","+width+ "]]"
    imfit_results[bcg_name] = imfit(file, region = region, residual = residual_file)
    exportfits(imagename = residual_file, fitsimage = residual_file + ".fits")

with open(residual_path + "imfit_results.pkl", "wb") as f:
    pickle.dump(imfit_results, f)

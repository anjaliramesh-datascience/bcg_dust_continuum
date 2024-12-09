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

path = "/Users/arames52/Desktop/CASA Imaging/MS_files/tclean_output/"
df = pd.read_csv("imfit_regions.csv")

def read_fits(path):

    hdu = fits.open(path)[0]
    data = hdu.data[0,0,:,:]
    header = hdu.header
    wcs = WCS(header, naxis = 2)
    data[np.where(np.isnan(data))] = 0

    return data, header

bcgs = ['CDFS-18','CDFS19','ES1-18','ES1-25','ES1_z_0.99','ES1_z_0.99b','ES1_z_1.04','ES1_z_1.38','ES1_z_1.60','ES1_z_1.65','ES1_z_1.70','XMM-113','XMM-11','XMM-29','XMM-30','XMM_z_0.81','XMM_z_0.9', 'XMM_z_1.0']
imfit_res = {}
imfit_results = {}
imgs = ['natural/', 'briggs_05/']
tapered = ["tapered_07arcsec/", "tapered/", "tapered_2arcsec/", "tapered_3arcsec/"]
images = imgs + tapered
for img_type in images:
    for bcg in bcgs:
        img_path = path + img_type + bcg + "_" +img_type.split("/")[0] + ".image.pbcor"
        residual_file = path + img_type + bcg + "_imfit_residual.image"
        imfit_data = ast.literal_eval(df[df['BCG'] == bcg][img_type.split("/")[0]].values[0])
        ra = str(imfit_data[0]) + "deg"
        dec = str(imfit_data[1]) + "deg"
        radius = str(imfit_data[2]) + "arcsec"
        region = "circle[[" + ra + "," + dec + "]," + radius + "]"
        print(bcg, img_type)
        imfit_res[bcg] = imfit(img_path, region = region, residual = residual_file)
        exportfits(imagename = residual_file, fitsimage = residual_file + ".fits")
    imfit_results[img_type.split("/")[0]] = imfit_res
    imfit_res = {}

with open("imfit_results_.pkl", "wb") as f:
    pickle.dump(imfit_results, f)

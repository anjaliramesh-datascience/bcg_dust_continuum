from astropy.io import fits
import pandas as pd
from astropy.table import Table
import numpy as np
import sys
import pickle
import glob
import astropy.units as u
from astropy.coordinates import SkyCoord
sys.path.append("/Users/arames52/bcg_dust_continuum/src/")
import bcg_parameter_file as p
import data_ingestion as dload

alma_images_path = "/Users/arames52/bcg_dust_continuum/notebook/data/ALMA_images/"

def signal_rms(i):

    bcg_image_file = glob.glob(alma_images_path + bcgs[i] + "*.fits")[0]

    x, y = str(ras[i]) + "deg", str(decs[i]) + "deg"
    r = str(round(radius[i],2)) + "arcsec"
    r_out = str(round(radius[i] + 3,2)) + "arcsec"
    region = "circle[["+ x +","+y+ "]," + r+ "]"
    annulus = "annulus[[" + x + "," + y + "], [" + r + "," + r_out + "]]"

    signal = imstat(imagename = bcg_image_file, region = region)
    max_sig = signal['max'][0]
    noise = imstat(imagename = bcg_image_file, region = annulus)
    rms = noise['rms'][0]

    return max_sig, rms

bcgs = list(dload.bcg_regions_load()['id'])
ras = list(dload.bcg_regions_load()['ra'])
decs = list(dload.bcg_regions_load()['dec'])
radius = list(dload.bcg_regions_load()['radius'])

s_n_dict = {}

for i in range(len(bcgs)):
    s_n_dict[bcgs[i]] = (signal_rms(i))

with open(p.data_path + "signal_noise_dict.pkl", "wb") as f:
    pickle.dump(s_n_dict, f)


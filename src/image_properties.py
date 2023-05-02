from astropy.io import fits
import pandas as pd
from astropy.table import Table
import numpy as np
import sys
import pickle
import glob
import astropy.units as u
from astropy.coordinates import SkyCoord

alma_images_path = "/Users/arames52/bcg_dust_continuum/notebook/data/MW_images/ALMA_images/"
bcg_regions = "/Users/arames52/bcg_dust_continuum/notebook/data/Derived_Data/bcg_regions.txt"

bcg_regions_df = pd.read_csv(bcg_regions, sep = ',')

bcgs = list(bcg_regions_df['id'])
ras = list(bcg_regions_df['ra'])
decs = list(bcg_regions_df['dec'])
radius = list(bcg_regions_df['radius'])
contour_steps = {'CDFS-18':25, 'CDFS19':1, 'ES1-12':1,'ES1-26':1, 'ES1-35':1, 'XMM-19':1, 'XMM-27':1,
                'ES1-18':1, 'ES1-25':1, 'ES1_z_0.88':1, 'ES1_z_0.99':1, 'ES1_z_0.99b':1,
                'ES1_z_1.60':1, 'XMM-11':1, 'XMM-29':1, 'XMM-30':1,'XMM-113':1, 'XMM_z_0.81':1,
                'XMM_z_0.9':1, 'XMM_z_1.0':1, 'ES1-34':2, 'ES1_z_1.04':2, 'ES1_z_1.38':2, 'ES1_z_1.40':2,
                'ES1_z_1.70':4, 'ES1_z_1.65':10}

def file_from_string(path, bcg):

    return glob.glob(path + bcg + "*.fits")[0]

def signal_rms(bcg_name):

    bcg_image_file = file_from_string(alma_images_path, bcg_name)
    ind = bcgs.index(bcg_name)

    # Defining region parameters
    x, y = str(ras[ind]) + "deg", str(decs[ind]) + "deg"
    r = str(round(radius[ind],1)) + "arcsec"
    r_out = str(round(radius[ind] + 3,1)) + "arcsec"

    # Region and annulus string
    region = "circle[["+ x +","+y+ "]," + r+ "]"
    annulus = "annulus[[" + x + "," + y + "], [" + r + "," + r_out + "]]"

    # Calling IMSTAT task to calculate image statistics
    signal = imstat(imagename = bcg_image_file, region = region)
    max_sig = round(signal['max'][0],7)
    noise = imstat(imagename = bcg_image_file, region = annulus)
    rms = round(noise['rms'][0],7)
    s_n = round(max_sig/rms,2)
    dex = contour_steps[bcg_name]

    return s_n, rms, dex

def run_all():
    stat_table_list = []
    for bcg in bcgs:
        print(bcg)
        s_n, rms,dex = signal_rms(bcg)
        stat_table_list.append({"id" : bcg, "S/N" : s_n, "RMS" : rms, "contour_steps":dex})
    stat_df = pd.DataFrame(stat_table_list)
    return stat_df

# df = run_all()
# df.to_csv("/Users/arames52/bcg_dust_continuum/notebook/data/Derived_Data/bcg_properties.csv", index = False)


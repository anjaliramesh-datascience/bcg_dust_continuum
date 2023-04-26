import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import pdb

def bcg_regions_load():
    file = "/Users/arames52/bcg_dust_continuum/notebook/data/bcg_regions.txt"
    bcg_regions = pd.read_csv(file, delimiter=',')
    return bcg_regions

def read_alma_fits(file):

    hdu = fits.open(file)[0]
    header = hdu.header
    data = hdu.data[0,0,:,:]
    wcs = WCS(header, naxis = 2)

    data = np.nan_to_num(data)

    return data,header,wcs

def read_any_fits(file):
    hdu = fits.open(file)[0]
    header = hdu.header
    data = hdu.data
    wcs = WCS(header)
    data = np.nan_to_num(data)
    return data, header, wcs

    







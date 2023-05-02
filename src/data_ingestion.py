import pandas as pd
import numpy as np
from dl import queryClient as qc
from PyAstronomy import pyasl
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys
sys.path.append('/Users/arames52/bcg_dust_continuum/src/')
import bcg_parameter_file as p
import pickle
import pycurl
import urllib.request
import tarfile

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

    
# bcg_coordinates = pd.read_csv(p.data_path + 'BCG_coords.txt', delim_whitespace=True, names = ['bcg', 'ra', 'dec'])

# DES catalog

# des_flux_catalog = {}

# def download_des_catalog(ra:float, dec:float, fov = 0.005):

#     """
#     ra - right ascention of target object
#     dec - declination of target object
#     fov - field of view

#     Returns a dataframe df with DES grizy band magnitudes and magnitude errors
#     """

#     query = "SELECT ra, dec, mag_auto_g, magerr_auto_g, mag_auto_i, magerr_auto_i, mag_auto_r, magerr_auto_r, mag_auto_y, magerr_auto_y, mag_auto_z, magerr_auto_z FROM des_dr2.mag WHERE 't' = Q3C_RADIAL_QUERY(ra,dec," + str(ra) + "," + str(dec) + "," + str(fov) + ")"
#     df = qc.query(sql=query,fmt='pandas')
#     # df['angular_distance'] = pyasl.getAngDist(row['ra'], row['dec'], np.array(df['ra']), np.array(df['dec']))
#     # df = df.sort_values(by = ['angular_distance'], ascending = True).reset_index(drop = True)
    
#     return df

# # Spitzer SWIRE catalog

# def download_swire_catalog(catalog_name, ra, dec, bcg):
    
#     """
#     catalog_name - IRSA catalog name
#     """
#     swire_outfile = swire_outpath + bcg + ".csv"
#     query = "https://irsa.ipac.caltech.edu/SCS?table="+catalog_name+"&RA="+str(ra)+"&DEC="+str(dec)+"&SR=0.005&format=csv"
#     with open(swire_outfile, "wb") as fp:
#         curl = pycurl.Curl()
#         curl.setopt(pycurl.URL, query)
#         curl.setopt(pycurl.WRITEDATA, fp)
#         curl.perform()
#         curl.close()
    
#     return None

# # Herschel HerMES Catalog

# def download_hermes_catalog(thetarfile):
#     """
#     download tarfile and untar using urlib.request and tarfile module
#     input - tar file download link
#     """
#     ftpstream = urllib.request.urlopen(thetarfile)
#     thetarfile = tarfile.open(fileobj=ftpstream, mode="r|bz2")
#     thetarfile.extractall(hermes_outpath)
#     return None


# # Querying for each source and downloading

# for ind, row in bcg_coordinates.iterrows():
#     des_flux_catalog[row['bcg']] = download_des_catalog(row['ra'], row['dec'], 0.005) # DES magnitudes extraction

#     # IRSA catalog name for the different sky areas
#     if row['bcg'].startswith('CDFS'):
#         irsa_catalog_name = "chandra_cat_f05"
#     elif row['bcg'].startswith('ES1'):
#         irsa_catalog_name = "elaiss1_cat_f05"
#     else:
#         irsa_catalog_name = "xmm_cat_s05"
    
#     download_swire_catalog(irsa_catalog_name, row['ra'], row['dec'], row['bcg']) # SWIRE fluxes download as csv files


# # Saving DES catalog as a pickle file
# with open(p.output_path + 'des_catalog.pkl', 'wb') as f:
#     pickle.dump(des_flux_catalog, f)

# # Downloading HerMES data for 24 micron positions
# for thetarfile in p.hermes_urls:
#     download_hermes_catalog(thetarfile)



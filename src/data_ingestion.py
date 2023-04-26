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

    
bcg_coordinates = pd.read_csv(p.data_path + 'BCG_coords.txt', delim_whitespace=True, names = ['bcg', 'ra', 'dec'])

# DES catalog

des_flux_catalog = {}

def download_des_catalog(ra:float, dec:float, fov = 0.005):

    """
    ra - right ascention of target object
    dec - declination of target object
    fov - field of view

    Returns a dataframe df with DES grizy band magnitudes and magnitude errors
    """

    query = "SELECT ra, dec, mag_auto_g, magerr_auto_g, mag_auto_i, magerr_auto_i, mag_auto_r, magerr_auto_r, mag_auto_y, magerr_auto_y, mag_auto_z, magerr_auto_z FROM des_dr2.mag WHERE 't' = Q3C_RADIAL_QUERY(ra,dec," + str(ra) + "," + str(dec) + "," + str(fov) + ")"
    df = qc.query(sql=query,fmt='pandas')
    # df['angular_distance'] = pyasl.getAngDist(row['ra'], row['dec'], np.array(df['ra']), np.array(df['dec']))
    # df = df.sort_values(by = ['angular_distance'], ascending = True).reset_index(drop = True)
    
    return df

# Spitzer SWIRE catalog

def download_swire_catalog(catalog_name, ra, dec, bcg):
    
    """
    catalog_name - IRSA catalog name
    """
    swire_outfile = swire_outpath + bcg + ".csv"
    query = "https://irsa.ipac.caltech.edu/SCS?table="+catalog_name+"&RA="+str(ra)+"&DEC="+str(dec)+"&SR=0.005&format=csv"
    with open(swire_outfile, "wb") as fp:
        curl = pycurl.Curl()
        curl.setopt(pycurl.URL, query)
        curl.setopt(pycurl.WRITEDATA, fp)
        curl.perform()
        curl.close()
    
    return None

# Herschel HerMES Catalog

def download_hermes_catalog(thetarfile):
    """
    download tarfile and untar using urlib.request and tarfile module
    input - tar file download link
    """
    ftpstream = urllib.request.urlopen(thetarfile)
    thetarfile = tarfile.open(fileobj=ftpstream, mode="r|bz2")
    thetarfile.extractall(hermes_outpath)
    return None


# Querying for each source and downloading

for ind, row in bcg_coordinates.iterrows():
    des_flux_catalog[row['bcg']] = download_des_catalog(row['ra'], row['dec'], 0.005) # DES magnitudes extraction

    # IRSA catalog name for the different sky areas
    if row['bcg'].startswith('CDFS'):
        irsa_catalog_name = "chandra_cat_f05"
    elif row['bcg'].startswith('ES1'):
        irsa_catalog_name = "elaiss1_cat_f05"
    else:
        irsa_catalog_name = "xmm_cat_s05"
    
    download_swire_catalog(irsa_catalog_name, row['ra'], row['dec'], row['bcg']) # SWIRE fluxes download as csv files


# Saving DES catalog as a pickle file
with open(p.output_path + 'des_catalog.pkl', 'wb') as f:
    pickle.dump(des_flux_catalog, f)

# Downloading HerMES data for 24 micron positions
for thetarfile in p.hermes_urls:
    download_hermes_catalog(thetarfile)



# Read DES catalog
with open(p.catalogs_path + 'des_catalog.pkl', 'rb') as f:
    des_catalog = pickle.load(f)

mips_nn = {}
irac_nn = {}

des_catalog_all_sources = pd.DataFrame()
swire_catalog_all_sources = pd.DataFrame()
hermes_catalog_all_sources = pd.DataFrame()

mag_cols = ['mag_auto_g', 'mag_auto_r', 'mag_auto_i', 'mag_auto_z','mag_auto_y']
magerr_cols = ['magerr_auto_g', 'magerr_auto_r', 'magerr_auto_i','magerr_auto_z', 'magerr_auto_y']
swire_columns = ['ra','dec','flux_kr_36', 'uncf_kr_36','flux_kr_45',
    'uncf_kr_45', 'flux_kr_58', 'uncf_kr_58', 'flux_kr_80', 'uncf_kr_80',
    'flux_kr_24','uncf_kr_24']
hermes_columns = ['F24', 'e_F24', 'F250', 'et_F250', 'F350', 'et_F350', 'F500', 'et_F500']

def catalog_matching(c, df):

    ra = np.array(df['ra']) * u.degree
    dec = np.array(df['dec']) * u.degree
    catalog = SkyCoord(ra, dec, frame = 'icrs')
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)
    idx_all = c.separation(catalog) < 0.007*u.deg
    match = pd.DataFrame(df.iloc[idx]).T.reset_index(drop = True)
    nearest_neighbours = pd.DataFrame(df.iloc[idx_all]).reset_index(drop=True)

    return match, nearest_neighbours

def mag_to_mJy(mag):
    return 10**((8.90 - mag)/2.5)*1000

for ind, row in bcg_coordinates.iterrows():

    des_df = des_catalog[row['bcg']]
    c = SkyCoord(row['ra'] * u.degree, row['dec']*u.degree, frame = 'icrs')
    des_match, des_nn = catalog_matching(c, des_df)
    des_fluxes = des_match[['ra', 'dec']].copy()
    for err_col,col in zip(magerr_cols, mag_cols):
        des_fluxes[col] = mag_to_mJy(des_match[col])
        des_fluxes[err_col] = des_match[err_col]*des_fluxes[col]
    des_catalog_all_sources = pd.concat([des_catalog_all_sources, des_fluxes])

# Read SWIRE catalog

    swire_df = pd.read_csv(swire_outpath + row['bcg'] + '.csv')
    swire_df = swire_df[swire_columns]
    if len(swire_df) > 0:
        swire_match, swire_nn = catalog_matching(c, swire_df)
        swire_catalog_all_sources = pd.concat([swire_catalog_all_sources, pd.concat([swire_match[swire_columns[:2]], swire_match[swire_columns[2:]].mul(0.001)], axis = 1)])
    else:
        swire_match = pd.DataFrame([[0]*swire_catalog_all_sources.shape[1]],columns=swire_catalog_all_sources.columns)
        swire_catalog_all_sources = pd.concat([swire_catalog_all_sources, swire_match])
    irac_nn[row['bcg']] = swire_nn[['ra', 'dec', 'flux_kr_36', 'uncf_kr_36']]
# Read HerMES catalog
    if row['bcg'].startswith('CDFS'):
        hermes_file = "/Users/arames52/bcg_dust_continuum/notebook/data/catalogs/herschel_catalog/L5-CDFS-SWIRE_xID24-DR3/L5-CDFS-SWIRE_xID24_DR3.fits"
    elif row['bcg'].startswith('ES1'):
        hermes_file = "/Users/arames52/bcg_dust_continuum/notebook/data/catalogs/herschel_catalog/L6-ELAIS-S1-SWIRE_xID24-DR3/L6-ELAIS-S1-SWIRE_xID24_DR3.fits"
    else:
        hermes_file = "/Users/arames52/bcg_dust_continuum/notebook/data/catalogs/herschel_catalog/L6-XMM-LSS-SWIRE_xID24-DR3/L6-XMM-LSS-SWIRE_xID24_DR3.fits"
    hermes_hdu = fits.open(hermes_file)
    hermes_df = Table(hermes_hdu[1].data).to_pandas()
    hermes_df = hermes_df.rename(columns={'RA':'ra', 'Dec':'dec'})
    
    hermes_match, hermes_nn = catalog_matching(c, hermes_df)
    hermes_match = hermes_match[hermes_columns]
    hermes_catalog_all_sources = pd.concat([hermes_catalog_all_sources, hermes_match], ignore_index=True)
    mips_objects = hermes_nn[["ra", "dec", "F24", "e_F24"]]
    mips_nn[row['bcg']] = mips_objects

swire_catalog_all_sources['id'] = list(bcg_coordinates['bcg'])
des_catalog_all_sources['id'] = list(bcg_coordinates['bcg'])
hermes_catalog_all_sources['id'] = list(bcg_coordinates['bcg'])

swire_catalog_all_sources = swire_catalog_all_sources.reset_index(drop = True)
des_catalog_all_sources = des_catalog_all_sources.reset_index(drop = True)

# Rename columns suitable for Cigale filter names
des_catalog_all_sources = des_catalog_all_sources.rename(columns = {"mag_auto_g":"des_g", "magerr_auto_g": "des_g_err", "mag_auto_r": "des_r", "magerr_auto_r": "des_r_err",
                     "mag_auto_i": "des_i", "magerr_auto_i":"des_i_err", "mag_auto_z":"des_z", "magerr_auto_z":"des_z_err",
                     "mag_auto_y":"des_Y", "magerr_auto_y": "des_Y_err"})
swire_catalog_all_sources = swire_catalog_all_sources.rename(columns = {"flux_kr_36":"spitzer.irac.ch1", "uncf_kr_36": "spitzer.irac.ch1_err",
    "flux_kr_45": "spitzer.irac.ch2", 'uncf_kr_45': "spitzer.irac.ch2_err", 'flux_kr_58': "spitzer.irac.ch3", 'uncf_kr_58': "spitzer.irac.ch3_err", 
    'flux_kr_80': "spitzer.irac.ch4", 'uncf_kr_80': "spitzer.irac.ch4_err"})

hermes_catalog_all_sources['F24'] = hermes_catalog_all_sources['F24'].apply(lambda x: x*0.001)
hermes_catalog_all_sources['e_F24'] = hermes_catalog_all_sources['e_F24'].apply(lambda x: x*0.001)
hermes_catalog_all_sources = hermes_catalog_all_sources.rename(columns = {"F24":"spitzer.mips.24", "e_F24": "spitzer.mips.24_err", "F250": "herschel.spire.PSW", 
    "et_F250": "herschel.spire.PSW_err","F350": "herschel.spire.PMW", "et_F350":"herschel.spire.PMW_err", 
    "F500": "herschel.spire.PLW", "et_F500":"herschel.spire.PLW_err"})

with open(p.data_path + "mips_prior_positions.pkl", "wb") as f:
    pickle.dump(mips_nn, f)
with open(p.data_path + "irac_prior_posisions.pkl", "wb") as f:
    pickle.dump(irac_nn, f)
# ALMA flux



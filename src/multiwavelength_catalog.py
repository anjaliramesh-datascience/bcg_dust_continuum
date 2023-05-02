import sys
sys.path.append("/Users/arames52/bcg_dust_continuum/src/")
# import bcg_parameter_file as p
# import helper_functions as hf
# import data_ingestion as di
import glob
import statmorph
import photutils
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.io import fits
import scipy.ndimage as ndi
from astropy.wcs import WCS
from scipy import stats
from math import *
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from matplotlib.patches import Rectangle, Ellipse, Circle
from astropy.visualization import simple_norm, ZScaleInterval, MinMaxInterval
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization.stretch import LinearStretch, LogStretch, SqrtStretch
from statmorph.utils.image_diagnostics import make_figure
import warnings
warnings.filterwarnings("ignore")
from astropy.wcs import FITSFixedWarning
warnings.filterwarnings('ignore', category=FITSFixedWarning)
from astropy.coordinates import Angle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.cosmology import FlatLambdaCDM
from regions import PixCoord
from regions import CircleAnnulusSkyRegion, CircleAnnulusPixelRegion
from regions import CircleSkyRegion, CirclePixelRegion
from ast import literal_eval
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
import matplotlib
plt.style.use(['science'])
import pickle
from astropy.wcs.utils import skycoord_to_pixel
import ast
from matplotlib import animation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pylab import figure, cm
from photutils.aperture import SkyCircularAperture
from photutils.aperture import aperture_photometry
from os.path import exists
plt.rcParams.update({'figure.max_open_warning': 1000})

# Paths

des_catalog_path = "/Users/arames52/bcg_dust_continuum/notebook/data/Catalogs/des_catalog.pkl"
swire_catalog_path = "/Users/arames52/bcg_dust_continuum/notebook/data/Catalogs/SWIRE_catalog/"
herschel_catalog_path = "/Users/arames52/bcg_dust_continuum/notebook/data/Catalogs/herschel_catalog/"
alma_catalog_path = "/Users/arames52/bcg_dust_continuum/notebook/data/Derived_Data/natural_imfit_results.pkl"
bcg_redshift_path = "/Users/arames52/bcg_dust_continuum/notebook/data/Derived_Data/BCG_redshifts.xlsx"
bcg_coordinates_path = "/Users/arames52/bcg_dust_continuum/notebook/data/Derived_Data/bcg_regions.txt"
sparcs_catalog_path = "/Users/arames52/bcg_dust_continuum/notebook/data/Catalogs/sparcs_catalog/"

# Column names
mag_cols = ['mag_auto_g', 'mag_auto_r', 'mag_auto_i', 'mag_auto_z','mag_auto_y']
magerr_cols = ['magerr_auto_g', 'magerr_auto_r', 'magerr_auto_i','magerr_auto_z', 'magerr_auto_y']
swire_columns = ['ra','dec','flux_kr_36', 'uncf_kr_36','flux_kr_45',
    'uncf_kr_45', 'flux_kr_58', 'uncf_kr_58', 'flux_kr_80', 'uncf_kr_80',
    'flux_kr_24','uncf_kr_24']
hermes_columns = ['ra', 'dec', 'F24', 'e_F24', 'F250', 'et_F250', 'F350', 'et_F350', 'F500', 'et_F500']

cigale_columns = ['id', 'redshift', 'des_g', 'des_r', 'des_i', 'des_z', 'des_Y',
       'des_g_err', 'des_r_err', 'des_i_err', 'des_z_err', 'des_Y_err',
       'spitzer.irac.ch1', 'spitzer.irac.ch1_err', 'spitzer.irac.ch2',
       'spitzer.irac.ch2_err', 'spitzer.irac.ch3', 'spitzer.irac.ch3_err',
       'spitzer.irac.ch4', 'spitzer.irac.ch4_err', 'spitzer.mips.24',
       'spitzer.mips.24_err', 'herschel.spire.PSW', 'herschel.spire.PSW_err',
       'herschel.spire.PMW', 'herschel.spire.PMW_err', 'herschel.spire.PLW',
       'herschel.spire.PLW_err', 'ALMA', 'ALMA_err']
sparcs_columns = ["ra", "dec", "z", "z_err", "irac1", "irac1_err"]

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

def vega_to_mjy(arr,err):
    arr = arr + 2.78
    mjy_flux = 10**((8.90 - arr)/2.5)*1000
    mjy_err = mjy_flux * np.sqrt((err/arr)**2)
    return mjy_flux, mjy_err

def bcg_regions_load():
    file = bcg_coordinates_path
    bcg_regions = pd.read_csv(file, delimiter=',')
    return bcg_regions


def multiwavelength_catalog(bcg_name):

    bcg_coordinates = bcg_regions_load()
    bcg_info = bcg_coordinates[bcg_coordinates['id'] == bcg_name]
    with open(des_catalog_path, "rb") as f:
        des_catalog = pickle.load(f)
    des_df = des_catalog[bcg_name]
    ra = bcg_info['ra'].values[0]
    dec = bcg_info['dec'].values[0]

    c = SkyCoord(ra * u.degree, dec * u.degree, frame  ='icrs')
    des_match, des_nn = catalog_matching(c, des_df)
    des_fluxes = des_match[['ra', 'dec']].copy()
    for err_col,col in zip(magerr_cols, mag_cols):
        des_fluxes[col] = mag_to_mJy(des_match[col])
        des_fluxes[err_col] = des_match[err_col]*des_fluxes[col]
    
    des_fluxes = des_fluxes.rename(columns = {"mag_auto_g":"des_g", "magerr_auto_g": "des_g_err", "mag_auto_r": "des_r", "magerr_auto_r": "des_r_err",
                     "mag_auto_i": "des_i", "magerr_auto_i":"des_i_err", "mag_auto_z":"des_z", "magerr_auto_z":"des_z_err",
                     "mag_auto_y":"des_Y", "magerr_auto_y": "des_Y_err"})

    swire_df = pd.read_csv(swire_catalog_path + bcg_name + ".csv")
    swire_df = swire_df[swire_columns]
    swire_match, swire_nn = catalog_matching(c, swire_df)
    swire_match = pd.concat([swire_match[swire_columns[:2]], swire_match[swire_columns[2:]].mul(0.001)], axis = 1)
    swire_match = swire_match.rename(columns = {"flux_kr_36":"spitzer.irac.ch1", "uncf_kr_36": "spitzer.irac.ch1_err",
    "flux_kr_45": "spitzer.irac.ch2", 'uncf_kr_45': "spitzer.irac.ch2_err", 'flux_kr_58': "spitzer.irac.ch3", 'uncf_kr_58': "spitzer.irac.ch3_err", 
    'flux_kr_80': "spitzer.irac.ch4", 'uncf_kr_80': "spitzer.irac.ch4_err"})
    
    if bcg_name.startswith('CDFS'):
        sparcs_file = "/Users/arames52/bcg_dust_continuum/notebook/data/Catalogs/sparcs_catalog/CDFS_zIRAC_simple.lst"
        hermes_file = "/Users/arames52/bcg_dust_continuum/notebook/data/catalogs/herschel_catalog/L5-CDFS-SWIRE_xID24-DR3/L5-CDFS-SWIRE_xID24_DR3.fits"
    elif bcg_name.startswith('ES1'):
        sparcs_file = "/Users/arames52/bcg_dust_continuum/notebook/data/Catalogs/sparcs_catalog/ELAISS1_zIRAC_simple.lst"
        hermes_file = "/Users/arames52/bcg_dust_continuum/notebook/data/catalogs/herschel_catalog/L6-ELAIS-S1-SWIRE_xID24-DR3/L6-ELAIS-S1-SWIRE_xID24_DR3.fits"
    else:
        sparcs_file = "/Users/arames52/bcg_dust_continuum/notebook/data/Catalogs/sparcs_catalog/XMM-LSS_zIRAC_simple.lst"
        hermes_file = "/Users/arames52/bcg_dust_continuum/notebook/data/catalogs/herschel_catalog/L6-XMM-LSS-SWIRE_xID24-DR3/L6-XMM-LSS-SWIRE_xID24_DR3.fits"
    
    sparcs_df = pd.read_csv(sparcs_file, delim_whitespace = True, comment = '#',names = ['ra','dec','APMAG_176','EAPMAG_176','APMAG_366','EAPMAG_366','MAG_BEST','EMAG_BEST','CH1BF','C1EBAP','C1_EPXAPC','C13PXERR','CH2BF','C2EBAP','C2_3PXAPC','C2_3PXERR' ])
    sparcs_match, sparcs_nn = catalog_matching(c, sparcs_df)
    sparcs_match = sparcs_match.rename(columns={"APMAG_366":"z", "EAPMAG_366":"z_err", "CH1BF":"irac1", "C1EBAP":"irac1_err"})
    sparcs_match = sparcs_match[sparcs_columns]
    sparcs_match['z'], sparcs_match['z_err'] = vega_to_mjy(np.array(sparcs_match['z']), np.array(sparcs_match['z_err']))
    sparcs_match['irac1'], sparcs_match['irac1_err'] = vega_to_mjy(np.array(sparcs_match['irac1']), np.array(sparcs_match['irac1_err']))

    hermes_hdu = fits.open(hermes_file)
    hermes_df = Table(hermes_hdu[1].data).to_pandas()
    hermes_df = hermes_df.rename(columns={'RA':'ra', 'Dec':'dec'})
    hermes_match, hermes_nn = catalog_matching(c, hermes_df)
    hermes_match = hermes_match[hermes_columns]
    hermes_match['F24'] = hermes_match['F24'].apply(lambda x: x*0.001)
    hermes_match['e_F24'] = hermes_match['e_F24'].apply(lambda x: x*0.001)
    hermes_match = hermes_match.rename(columns = {"F24":"spitzer.mips.24", "e_F24": "spitzer.mips.24_err", "F250": "herschel.spire.PSW", 
    "et_F250": "herschel.spire.PSW_err","F350": "herschel.spire.PMW", "et_F350":"herschel.spire.PMW_err", 
    "F500": "herschel.spire.PLW", "et_F500":"herschel.spire.PLW_err"})

    # print(des_fluxes)

    with open(alma_catalog_path, "rb") as f:
        alma_flux_results = pickle.load(f)

    redshift_info = pd.read_excel(bcg_redshift_path)
    redshift = redshift_info[redshift_info['bcg'] == bcg_name]['redshift'].values[0]
    mw_flux_dict = {"id":bcg_name, "redshift":redshift}
    for col in cigale_columns[2:12]:
        mw_flux_dict[col] = des_fluxes[col].values[0]
    for col in cigale_columns[12:20]:
        mw_flux_dict[col] = swire_match[col].values[0]
    for col in cigale_columns[20:28]:
        mw_flux_dict[col] = hermes_match[col].values[0]

    # if bcg_name in ["CDFS19", "ES1"]

    alma_flux_bcg = alma_flux_results[bcg_name]['deconvolved']['component0']['flux']
    mw_flux_dict['ALMA'] = alma_flux_bcg['value'][0] * 1000
    mw_flux_dict['ALMA_err'] = alma_flux_bcg['error'][0] * 1000

    mw_nn = {"id" : bcg_name, "des_nn" : des_nn[['ra', 'dec']], "irac_nn" : swire_nn[['ra', 'dec']], "mips_nn": hermes_nn[['ra', 'dec']], "sparcs_nn":sparcs_nn[['ra', 'dec']]}

    return mw_flux_dict, mw_nn, sparcs_match


# good_detections = ["CDFS-18", "ES1-18", "ES1-25", "ES1_z_0.99","ES1_z_0.99b","ES1_z_1.04","ES1_z_1.38","ES1_z_1.40",
# "ES1_z_1.60", "ES1_z_1.65", "ES1_z_1.70", "XMM-113", "XMM-29", "XMM-30", "XMM_z_0.9", "XMM_z_1.0"]
# weak_detections = ['ES1-12', 'ES1-26', 'XMM-19' , 'XMM-27','CDFS19', 'ES1_z_0.88', 'ES1-35', 'XMM_z_0.81']

# for i in range(len(good_detections)):
#     mw_catalog = multiwavelength_catalog(good_detections[i])
#     df = pd.DataFrame(mw_catalog, index = [i])
#     df.to_csv("/Users/arames52/bcg_dust_continuum/notebook/data/Derived_Data/" + good_detections[i] + "_flux_catalog.csv")

# for bcg in weak_detections:
#     mw_catalog = multiwavelength_catalog(bcg)
#     df = pd.DataFrame(mw_catalog)
#     df.to_csv("/Users/arames52/bcg_dust_continuum/notebook/data/Derived_Data/" + bcg + "_flux_catalog.csv")
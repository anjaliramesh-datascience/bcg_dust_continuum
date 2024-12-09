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
from image_processing import snr_rms
from image_processing import read_fits
from os.path import exists
plt.rcParams.update({'figure.max_open_warning': 1000})

# Path and column definitions
catalog_path = "/Users/arames52/Research/Data/Catalogs/"
sparcs_columns = ['RA','DEC','APMAG_176','EAPMAG_176','APMAG_366','EAPMAG_366','MAG_BEST',
'EMAG_BEST','CH1BF','C1EBAP','C1_EPXAPC','C13PXERR','CH2BF','C2EBAP','C2_3PXAPC','C2_3PXERR']
mag_cols = ['mag_auto_g', 'mag_auto_r', 'mag_auto_i', 'mag_auto_z','mag_auto_y']
magerr_cols = ['magerr_auto_g', 'magerr_auto_r', 'magerr_auto_i','magerr_auto_z', 'magerr_auto_y']
flux_cols = ['des_g', 'des_r', 'des_i', 'des_z','des_y']
fluxerr_cols = ['des_g_err', 'des_r_err', 'des_i_err', 'des_z_err','des_y_err']

# Optical - DES

def des_catalog_matching(bcg, ra, dec):

    des_catalog = pd.read_csv(catalog_path + "DES/" + bcg + ".txt", delimiter=",")

    c1 = SkyCoord(ra, dec, frame = 'fk5', unit = 'deg')
    c2 = SkyCoord(np.array(des_catalog['ra']), np.array(des_catalog['dec']), frame = 'fk5', unit = 'deg')

    des_catalog['dist'] = c1.separation(c2).arcsecond

    des_catalog = des_catalog[des_catalog['dist'] < 1].reset_index(drop = True)

    def mag_to_mJy(mag):
        return 10**((8.90 - mag)/2.5)*1000

    for err_col,col, new_col, new_err in zip(magerr_cols, mag_cols, flux_cols, fluxerr_cols):
        des_catalog[new_col] = mag_to_mJy(des_catalog[col])
        des_catalog[new_err] = mag_to_mJy(des_catalog[err_col])

    return des_catalog

def sparcs_catalog_matching(bcg, ra, dec):

    if bcg.startswith("CDFS"):
        sparcs_catalog = pd.read_csv(catalog_path + "Sparcs/" + "CDFS.lst", comment = '#', names = sparcs_columns,
        delim_whitespace = True)
    elif bcg.startswith("ES1"):
         sparcs_catalog = pd.read_csv(catalog_path + "Sparcs/" + "ES1.lst", comment = '#', names = sparcs_columns,
         delim_whitespace = True)
    elif bcg.startswith("XMM"):
        sparcs_catalog = pd.read_csv(catalog_path + "Sparcs/" + "XMM.lst", comment = '#', names = sparcs_columns,
        delim_whitespace = True)

    c1 = SkyCoord(ra, dec, frame = 'fk5', unit = 'deg')
    c2 = SkyCoord(np.array(sparcs_catalog['RA']), np.array(sparcs_catalog['DEC']), frame = 'fk5', unit = 'deg')

    sparcs_catalog['dist'] = c1.separation(c2).arcsecond
    sparcs_catalog = sparcs_catalog[sparcs_catalog['dist'] < 1].reset_index(drop = True)

    return sparcs_catalog



def irac_catalog_matching(bcg, ra, dec):

    if bcg.startswith("CDFS"):
        irac_catalog = pd.read_csv(catalog_path + "Swire/" + "CDFS.csv")
    elif bcg.startswith("ES1"):
         irac_catalog = pd.read_csv(catalog_path + "Swire/" + "ES1.csv")
    elif bcg.startswith("XMM"):
        irac_catalog = pd.read_csv(catalog_path + "Swire/" + "XMM.csv")

    c1 = SkyCoord(ra, dec, frame = 'fk5', unit = 'deg')
    c2 = SkyCoord(np.array(irac_catalog['ra']), np.array(irac_catalog['dec']), frame = 'fk5', unit = 'deg')

    irac_catalog['dist'] = c1.separation(c2).arcsecond

    irac_catalog = irac_catalog[irac_catalog['dist'] <1].reset_index(drop = True)

    return irac_catalog

def hermes_catalog_matching(bcg, ra, dec):

    if bcg.startswith("CDFS"):
        fits_table = fits.open(catalog_path + "Hermes/" + "CDFS.fits")
        hermes_catalog = Table(fits_table[1].data).to_pandas()
    elif bcg.startswith("ES1"):
         fits_table = fits.open(catalog_path + "Hermes/" + "ES1.fits")
         hermes_catalog = Table(fits_table[1].data).to_pandas()
    elif bcg.startswith("XMM"):
        fits_table = fits.open(catalog_path + "Hermes/" + "XMM.fits")
        hermes_catalog = Table(fits_table[1].data).to_pandas()

    c1 = SkyCoord(ra, dec, frame = 'fk5', unit = 'deg')
    c2 = SkyCoord(np.array(hermes_catalog['RA']), np.array(hermes_catalog['Dec']), frame = 'fk5', unit = 'deg')

    hermes_catalog['dist'] = c1.separation(c2).arcsecond

    hermes_catalog = hermes_catalog[hermes_catalog['dist'] < 1].reset_index(drop = True)
    hermes_catalog['ra'] = hermes_catalog['RA']
    hermes_catalog['dec'] = hermes_catalog['Dec']
    return hermes_catalog

def normalize_2d(matrix):
    norm = np.linalg.norm(matrix)
    matrix = matrix/norm  # normalized matrix
    return matrix


def plot_bcg_with_coords(bcg, ra, dec):

    image_path = "/Users/arames52/Research/Data/Images/"

    """
    DES and Sparcs z band, IRAC 3.6 micron, MIPS 24 micron and ALMA images postage stamps plotted
    with ALMA wcs
    """
    des_file = glob.glob(image_path + "DES/" + bcg + "/*_z_*.fits")[0]
    sparcs_file = glob.glob(image_path + "Sparcs_z/" + bcg + ".fits")[0]
    irac_file = glob.glob(image_path + "IRAC1/" + bcg + "*.fits")[0]
    mips_file = glob.glob(image_path + "MIPS1/" + bcg + "*.fits")[0]
    alma_file = glob.glob(image_path + "ALMA/" + bcg + "*.fits")[0]

    des_z_h, des_z_d, des_z_wcs = read_fits(des_file)
    sparcs_z_h, sparcs_z_d, sparcs_z_wcs = read_fits(sparcs_file)
    irac_h, irac_d, irac_wcs = read_fits(irac_file)
    mips_h, mips_d, mips_wcs = read_fits(mips_file)
    alma_h, alma_d, alma_wcs = read_fits(alma_file, 'alma')

    c = SkyCoord(ra, dec, frame = 'fk5', unit = 'deg')
    centerpix = alma_wcs.world_to_pixel(c)

    xlim = (centerpix[0]-9/0.045, centerpix[0]+9/0.045)
    ylim = (centerpix[1]-9/0.045, centerpix[1]+9/0.045)

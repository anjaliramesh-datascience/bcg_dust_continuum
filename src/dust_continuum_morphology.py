import sys
sys.path.append("/Users/arames52/bcg_dust_continuum/src/")
import bcg_parameter_file as p
import helper_functions as hf
import data_ingestion as di
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

image_path = "/Users/arames52/bcg_dust_continuum/notebook/data/MW_images/ALMA_images/"

def read_alma_fits(file):
    hdu = fits.open(file)[0]
    header = hdu.header
    data = hdu.data[0,0,:,:]
    wcs = WCS(header, naxis = 2)
    data = np.nan_to_num(data)
    return data,header,wcs

def image_cutout(file, size_in_arcsec):
    data,header,wcs = read_alma_fits(file)
    ra, dec = header['CRVAL1'], header['CRVAL2']
    position = SkyCoord(ra, dec, frame = 'fk5', unit = 'deg')
    s = size_in_arcsec
    size = u.Quantity((s,s), u.arcsec)
    cutout = Cutout2D(data, position = position, size = size, wcs = wcs)
    return cutout

def compute_sersic(bcg_name):

    file = glob.glob(image_path + bcg_name +"*.fits")[0]
    data, header, wcs = read_alma_fits(file)
    bmaj = header['BMAJ']*u.deg
    bmin = header['BMIN']*u.deg
    bpa = u.Quantity(header['BPA'], unit = "deg")
    cdelt = header['CDELT2']*u.deg
    x_sigma = bmin/cdelt
    y_sigma = bmaj/cdelt
    psf = Gaussian2DKernel(x_stddev=x_sigma*gaussian_fwhm_to_sigma, y_stddev=y_sigma*gaussian_fwhm_to_sigma, theta = bpa)
    image = image_cutout(file, 8).data   
    weightmap = np.full(image.shape, 2e-5)
    threshold = 2*2e-5
    segm = photutils.detect_sources(image, threshold, npixels = 5)
    label = np.argmax(segm.areas) + 1
    segmap = segm.data == label
    segmap_float = ndi.uniform_filter(np.float64(segmap), size=10)
    segmap = np.int64(segmap_float > 0.5)
    source_morphs = statmorph.source_morphology(image, segmap, weightmap = weightmap)
    morph = source_morphs[0]

    return morph

def save_sersic_results():
    sersic_profile = {}
    good_detections = ["CDFS-18", "ES1-18", "ES1-25", "ES1_z_0.99","ES1_z_0.99b","ES1_z_1.04","ES1_z_1.38","ES1_z_1.40",
"ES1_z_1.60", "ES1_z_1.65", "ES1_z_1.70", "XMM-113", "XMM-29", "XMM-30", "XMM_z_0.9", "XMM_z_1.0"]
    for bcg in good_detections:
        morph = compute_sersic(bcg)
        fig = make_figure(morph)
        fig.suptitle(bcg, x = 0.8)
        fig.savefig("/Users/arames52/bcg_dust_continuum/notebook/plots/" + bcg + "_sersic.png", dpi = 300)
        sersic_profile[bcg] = (morph.sersic_n, morph.sersic_rhalf * 0.045)
    return sersic_profile

sersic_profile_results = save_sersic_results()
with open("/Users/arames52/bcg_dust_continuum/notebook/data/Derived_Data/sersic_results.pkl", "wb") as f:
    pickle.dump(sersic_profile_results, f)
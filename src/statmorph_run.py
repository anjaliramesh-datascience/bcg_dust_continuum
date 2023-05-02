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


#ALMA Images
alma_images_path = "/Users/arames52/Research/Data/Images/ALMA/"
alma_images = sorted(glob.glob(alma_images_path + "*.fits"))
# Gaussian Residuals
residual_path = "/Users/arames52/Research/Data/GnS_profile/Gaussian/Residuals/"
residual_postfix = "_residual.image.fits"
#IMFIT Results
with open("/Users/arames52/Research/Data/GnS_profile/Gaussian/Residuals/imfit_results.pkl", "rb") as f:
    imfit_results = pickle.load(f)
#BCG Coordinates
bcg_coords_df = pd.read_csv("/Users/arames52/Research/Data/BCG_coords.txt", delim_whitespace = True,
names = ['bcg', 'ra', 'dec'])
bcg_coords = {}
for ind, row in bcg_coords_df.iterrows():
    bcg_coords[row['bcg']] = [row['ra'], row['dec']]
#RMS and SNR dict
with open("/Users/arames52/Research/Data/rms.pkl", "rb") as f:
    rms_dict = pickle.load(f)
with open("/Users/arames52/Research/Data/snr.pkl", "rb") as f:
    snr_dict = pickle.load(f)
# Function to read the fits image and make 3"x3" cutout image of the same
def read_fits_cutout(file, bcg):

    hdu = fits.open(file)[0]
    header = hdu.header
    data = hdu.data[0,0,:,:]
    wcs = WCS(header, naxis = 2)
    data = np.nan_to_num(data)
    data[data<0] = 0

    ra, dec = bcg_coords[bcg][0], bcg_coords[bcg][1]
    position = SkyCoord(ra, dec, frame = 'fk5', unit = 'deg')
    s = u.Quantity((3,3), u.arcsec)
    image_cutout = Cutout2D(data, position = position, size = s, wcs = wcs).data
    image_cutout[image_cutout<0] = 0

    bmaj = header['BMAJ']*u.deg.to(u.arcsec)/0.045
    bmin = header['BMIN']*u.deg.to(u.arcsec)/0.045
    pa = (header['BPA']  + 90)*u.deg.to(u.rad)
    psf = Gaussian2DKernel(x_stddev = gaussian_fwhm_to_sigma*bmin, y_stddev = gaussian_fwhm_to_sigma*bmaj, theta = pa, x_size = image_cutout.shape[0], y_size = image_cutout.shape[1]).array
    return data, header, wcs, image_cutout, psf
#Compute Sersic Profile
def compute_sersic(image, weightmap, psf, bcg):

    threshold = 3*rms_dict[bcg]
    segm = photutils.detect_sources(image, threshold, npixels = 5)
    label = argmax(segm.areas) + 1
    segmap = segm.data == label
    segmap_float = ndi.uniform_filter(float64(segmap), size=10)
    segmap = segmap_float > 0.5
    source_morphs = statmorph.source_morphology(image, segmap, weightmap = weightmap, psf = psf)
    morph = source_morphs[0]

    return morph

for file in alma_images:

    bcg = file.split(alma_images_path)[-1].split("_natural.fits")[0]
    data, header, wcs, image_cutout, psf = read_fits_cutout(file)
    weightmap = np.full(image_cutout.shape, rms_dict[bcg])
    morph = compute_sersic(image_cutout, weightmap, psf,bcg)


def plot_profiles():
    pass
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
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from matplotlib.patches import Rectangle, Ellipse, Circle
from astropy.visualization import simple_norm, ZScaleInterval
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization.stretch import LinearStretch, LogStretch
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
plt.style.use(['science',"default"])
import pickle
from pylab import figure, cm
from decimal import Decimal


with open("imfit_results.pkl", "rb") as f:
    imfit_results = pickle.load(f)

detections = ['CDFS-18','CDFS19','ES1-18','ES1-25','ES1_z_0.99','ES1_z_0.99b','ES1_z_1.04','ES1_z_1.38','ES1_z_1.60','ES1_z_1.65','ES1_z_1.70','XMM-113','XMM-11','XMM-29','XMM-30','XMM_z_0.81','XMM_z_0.9', 'XMM_z_1.0']


images_type = ["briggs_05", "natural","tapered_07arcsec", "tapered", "tapered_2arcsec", "tapered_3arcsec"]

for bcg in detections:
    fwhm = []
    flux = []
    error = []
    fig,ax = plt.subplots(1,1, figsize = (5,5))
    for w in  images_type:
        imfit_dict = imfit_results[w][bcg]['deconvolved']['component0']
        # if w in ["tapered_2arcsec", "tapered_3arcsec"]:
        #     flux.append(imfit_results[w][bcg]['results']['component0']['peak']['value'])
        #     error.append(imfit_results[w][bcg]['results']['component0']['peak']['error'])
        # else:
        flux.append(imfit_dict['flux']['value'][0])
        error.append(imfit_dict['flux']['error'][0])
        fwhm.append(imfit_results[w][bcg]['results']['component0']['beam']['beamarcsec']['major']['value'])

    flux_diff = []
    flux = np.array(flux)*1000
    error = np.array(error)*1000
    ax.scatter(fwhm, flux, color = 'red')
    ax.set_title(bcg)
    ax.errorbar(fwhm, flux, yerr = error,ls = 'none', ecolor = 'black', alpha = 0.5,solid_capstyle='projecting', capsize=5, elinewidth = 1)
    flux_difference = []
    for i in range(0,6):
        if i+1 >= 6:
            continue
        else:
            quad_err = np.sqrt(error[i]**2 + error[i+1]**2)
            fd = np.abs((flux[i+1]-flux[i]))/quad_err
            flux_diff.append(fd)

        str_i = str(i)
        text = str(round(fd,2)) + "$\sigma_" + str_i + "$"
        flux_difference.append(text)
        for text, x_pos in zip(flux_difference, np.array([0.05, 0.18, 0.32, 0.5, 0.8])):
            ax.text(x = x_pos, y = 0.1, s = text,transform=ax.transAxes)
    ax.set_xlabel('Beam FWHM [arcsec]', size = 15)
    ax.set_ylabel('Flux Density [mJy]',size = 15)
    # plt.show()
    fig.savefig("/Users/arames52/Desktop/Plots/check/" + bcg + "_fb.jpg", dpi = 300)

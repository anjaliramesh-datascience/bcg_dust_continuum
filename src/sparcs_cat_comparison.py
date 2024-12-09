import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from astropy.visualization import simple_norm
from astropy.modeling import models
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization.stretch import LinearStretch
from astropy.visualization import ZScaleInterval
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import units as u
import astropy
import photutils
import time
from astropy.stats import sigma_clipped_stats
import statmorph
from statmorph.utils.image_diagnostics import make_figure
import glob
from PIL import Image
import os
import pandas as pd
import numpy as np
import pyregion
from regions import Regions
from math import *
import sep
import scipy
import warnings
from astropy.wcs import FITSFixedWarning
warnings.filterwarnings('ignore')

sparcs_catalog = "/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/sparcs_cats/"
bcg_coords = pd.read_csv("/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/BCG_coords.txt")
flux_table = pd.read_csv("all_fluxes.csv")

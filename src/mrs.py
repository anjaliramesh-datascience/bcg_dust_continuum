import glob
# import statmorph
# import photutils
# import numpy as np
# import pandas as pd
# import astropy.units as u
# from astropy.io import fits
# import scipy.ndimage as ndi
# from astropy.wcs import WCS
# from scipy import stats
# import matplotlib.pyplot as plt
# import matplotlib.colors as colors
# from astropy.nddata import Cutout2D
# from astropy.coordinates import SkyCoord
# from astropy.convolution import Gaussian2DKernel
# from astropy.stats import gaussian_fwhm_to_sigma
# from matplotlib.patches import Rectangle, Ellipse, Circle
# from astropy.visualization import simple_norm, ZScaleInterval
# from astropy.visualization.mpl_normalize import ImageNormalize
# from astropy.visualization.stretch import LinearStretch, LogStretch
# from statmorph.utils.image_diagnostics import make_figure
# import warnings
# from astropy.wcs import FITSFixedWarning
# warnings.filterwarnings('ignore', category=FITSFixedWarning)
# from astropy.coordinates import Angle
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from astropy.cosmology import FlatLambdaCDM
# from regions import PixCoord
# from regions import CircleAnnulusSkyRegion, CircleAnnulusPixelRegion
# from ast import literal_eval
# cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
# plt.style.use(['science',"default"])
import sys
sys.path.append("/Users/arames52/Downloads/analysis_scripts/") # Add the path to where you have downloaded analysis utils in your computer
import analysisUtils as au
import pickle


measurement_sets = sorted(glob.glob("/Users/arames52/Desktop/CASA Imaging/MS_files/measurement_sets/*.ms"))
mrs_dict = {}
for file in measurement_sets:
    mrs = au.estimateMRS(file)
    name = file.split("/Users/arames52/Desktop/CASA Imaging/MS_files/measurement_sets/")[-1].split("_calibrated.ms")[0]
    mrs_dict[name] = mrs

with open("mrs.pkl","wb") as f:
    pickle.dump(mrs_dict, f)

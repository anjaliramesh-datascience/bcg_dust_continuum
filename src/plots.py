import matplotlib
import matplotlib.pyplot as plt
import astropy.units as u
from matplotlib.patches import Rectangle, Ellipse, Circle
from astropy.modeling.functional_models import Ellipse2D
from matplotlib.lines import Line2D
import matplotlib.colors as colors
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm, ZScaleInterval
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization.stretch import LinearStretch, LogStretch
import pandas as pd
import sys
sys.path.append("/Users/arames52/bcg_dust_continuum/src/")
import image_properties as img_prop
import data_ingestion as dload
import multiwavelength_catalog as mwcat
import cmocean
import numpy as np
import bcg_parameter_file as bpf

def figure_settings():
    pass
    

def plot_sfr_mstar():
    pass


# Image Path 

des_images_path = "/Users/arames52/bcg_dust_continuum/notebook/data/MW_images/DES_images/"
sparcs_images_path = "/Users/arames52/bcg_dust_continuum/notebook/data/MW_images/Sparcs_z_images/"
irac_images_path = "/Users/arames52/bcg_dust_continuum/notebook/data/MW_images/IRAC_images/"
mips_images_path = "/Users/arames52/bcg_dust_continuum/notebook/data/MW_images/MIPS_images/"
alma_images_path = "/Users/arames52/bcg_dust_continuum/notebook/data/MW_images/ALMA_images/"


# Function to plot bcg stamps - Sparcs z band, IRAC 3.6 microns, MIPS 24 microns and ALMA 1.2mm image

def multiwavelength_postage_stamps(bcg_name, size):

    # Get image file path from bcg name and main folder path
    z_img = img_prop.file_from_string(sparcs_images_path, bcg_name)
    irac_img = img_prop.file_from_string(irac_images_path, bcg_name)
    mips_img = img_prop.file_from_string(mips_images_path, bcg_name)
    alma_img = img_prop.file_from_string(alma_images_path, bcg_name)


    # Get fits data, header and wcs
    data_z, header_z, wcs_z = dload.read_any_fits(z_img)
    data_irac, header_irac, wcs_irac = dload.read_any_fits(irac_img)
    data_mips, header_mips, wcs_mips = dload.read_any_fits(mips_img)
    data_alma, header_alma, wcs_alma = dload.read_alma_fits(alma_img)

    # Initialize a figure
    fig = plt.figure(dpi = 300)

    # Add 4 subplots corresponding to the 4 multi-wavelength image
    ax1 = fig.add_subplot(141, projection=wcs_alma)
    ax2 = fig.add_subplot(142, projection=wcs_alma,sharey=ax1)
    ax3 = fig.add_subplot(143, projection=wcs_alma,sharey=ax1)
    ax4 = fig.add_subplot(144, projection=wcs_alma,sharey=ax1)

    # Get axes transformation
    z_trans = ax1.get_transform(wcs_z)
    irac_trans = ax2.get_transform(wcs_irac)
    mips_trans = ax3.get_transform(wcs_mips)

    # Define FWHM of PSF
    z_fwhm = 0.67 * u.arcsec
    irac_fwhm = 2 * u.arcsec
    mips_fwhm = 6 * u.arcsec

    # Define ALMA beam
    bmaj = round(header_alma['BMAJ']*u.degree.to(u.arcsec)/bpf.pixel_scale) # in units of pixels
    bmin = round(header_alma['BMIN']*u.degree.to(u.arcsec)/bpf.pixel_scale)
    bpa = round(header_alma['BPA']) + 90
    alma_fwhm_maj = header_alma['BMAJ'] * u.degree
    alma_fwhm_min = header_alma['BMIN'] * u.degree

    # Setting ALMA beam properties
    ra = header_alma['CRVAL1']
    dec = header_alma['CRVAL2']
    center = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5')
    extent = SkyCoord(ra*u.deg - alma_fwhm_maj, dec*u.deg + alma_fwhm_min, frame='fk5')
    centerpix = wcs_alma.world_to_pixel(center)
    extpix = wcs_alma.world_to_pixel(extent)
    arad, brad = (extpix[0] - centerpix[0]), (extpix[1] - centerpix[1])

    # Setting stamps cutout size - 10" x 10"
    stamp_size = size/bpf.pixel_scale # this is the pixel size - change for different pixel
    xlim = (centerpix[0]-stamp_size,centerpix[0]+stamp_size)
    ylim = (centerpix[1]-stamp_size,centerpix[1]+stamp_size)

    # ALMA contour levels
    bcg_properties = pd.read_csv("/Users/arames52/bcg_dust_continuum/notebook/data/Derived_Data/bcg_properties.csv")
    img_prop_df = bcg_properties[bcg_properties['id'] == bcg_name]
    s_n = img_prop_df['S/N'].values[0]
    rms = img_prop_df['RMS'].values[0]
    dex = img_prop_df['contour_steps'].values[0]
    levels = np.arange(3*rms, s_n*rms, dex*rms, dtype = None)

    # Plotting all images
    im1 = ax1.imshow(data_z, transform = z_trans, origin = 'lower',vmin = 0, vmax=200, cmap = cmocean.cm.thermal)
    im2 = ax2.imshow(data_irac, transform=irac_trans, origin = 'lower', vmin =0, vmax=1,cmap = cmocean.cm.thermal)
    im3 = ax3.imshow(data_mips, transform=mips_trans, origin='lower', vmin = 0, vmax = 0.5,cmap = cmocean.cm.thermal)
    im4 = ax4.imshow(data_alma, origin='lower', vmin = 0, vmax = 0.0001, cmap = cmocean.cm.thermal)

    # Plotting contours
    ax1.contour(data_alma, levels=levels, colors='red', linewidths = 0.2, zorder = 1)
    ax2.contour(data_alma, levels=levels, colors='red', linewidths = 0.2, zorder = 1)
    ax3.contour(data_alma, levels=levels, colors='red', linewidths = 0.2, zorder = 1)
    ax4.contour(data_alma, levels=levels, colors = 'red', linewidths = 0.2, zorder = 1, alpha = 1)

    # Source coordinates and plotting object centroids
    _, mw_cat, _ = mwcat.multiwavelength_catalog(bcg_name)
    for ind,row in mw_cat['sparcs_nn'].iterrows():
        c_z = SkyCoord(ra = row['ra']*u.degree, dec = row['dec']*u.degree, frame = 'icrs')
        c_pix_z = wcs_z.world_to_pixel(c_z)
        ax1.scatter(c_pix_z[0], c_pix_z[1], marker = '+', color = 'black', zorder = 10, transform = irac_trans, s = 2, lw= 0.5, label = 'sparcs sources')
    for ind,row in mw_cat['irac_nn'].iterrows():
        c_irac = SkyCoord(ra = row['ra']*u.degree, dec = row['dec']*u.degree, frame = 'icrs')
        c_pix_irac = wcs_irac.world_to_pixel(c_irac)
        ax2.scatter(c_pix_irac[0], c_pix_irac[1], marker = '+', color = 'black', zorder = 5, transform = irac_trans, s = 2, lw= 0.5, label = 'irac sources')
    for ind,row in mw_cat['mips_nn'].iterrows():
        c_mips = SkyCoord(ra = row['ra']*u.degree, dec = row['dec']*u.degree, frame = 'icrs')
        c_pix_mips = wcs_mips.world_to_pixel(c_mips)
        ax3.scatter(c_pix_mips[0], c_pix_mips[1], marker = '+', color = 'black', zorder = 5, transform = mips_trans, s = 2, lw= 0.5, label = 'mips sources')

    c_alma = SkyCoord(ra = header_alma['CRVAL1'] * u.deg, dec = header_alma['CRVAL2']*u.deg, frame = 'icrs')
    c_pix_alma = wcs_alma.world_to_pixel(c_alma)
    ax4.scatter(c_pix_alma[0], c_pix_alma[1], marker = '+', color = 'black', zorder = 5, s = 2, lw= 0.5, label = 'alma sources')

    # Matplotlib ellipse and circle patch to depict the different beam sizes


    ax1.set_ylabel("Dec")
    ax1.set_xlabel("RA")
    ax2.set_xlabel("RA")
    ax3.set_xlabel("RA")
    ax4.set_xlabel("RA")

    ax1.coords[1].set_ticks_visible(False)
    ax1.coords[1].set_ticklabel_visible(False)
    ax2.coords[1].set_ticks_visible(False)
    ax2.coords[1].set_ticklabel_visible(False)
    ax3.coords[1].set_ticks_visible(False)
    ax3.coords[1].set_ticklabel_visible(False)
    ax4.coords[1].set_ticks_visible(False)
    ax4.coords[1].set_ticklabel_visible(False)
    ax1.coords[0].set_ticks_visible(False)
    ax1.coords[0].set_ticklabel_visible(False)
    ax2.coords[0].set_ticks_visible(False)
    ax2.coords[0].set_ticklabel_visible(False)
    ax3.coords[0].set_ticks_visible(False)
    ax3.coords[0].set_ticklabel_visible(False)
    ax4.coords[0].set_ticks_visible(False)
    ax4.coords[0].set_ticklabel_visible(False)
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax2.set_xlim(xlim)
    ax2.set_ylim(ylim)
    ax3.set_xlim(xlim)
    ax3.set_ylim(ylim)
    ax4.set_xlim(xlim)
    ax4.set_ylim(ylim)
    ax1.set_title("SpARCS z-band", size = 5)
    ax2.set_title("IRAC 3.6 $\mu$m", size = 5)
    ax3.set_title("MIPS 24 $\mu$m", size = 5)
    ax4.set_title("ALMA 1.2 mm", size = 5)
    # fig.suptitle(bcg_name, y = 0.70)


    return fig



def imfit_plot():
    pass

def sfr_mstar():
    df = pd.read_csv("/Users/arames52/bcg_dust_continuum/notebook/data/Derived_Data/results_table.csv")
    pass


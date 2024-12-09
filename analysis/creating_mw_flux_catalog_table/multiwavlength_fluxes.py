from astropy.io import fits
import numpy as np
from astropy.table import Table
from decimal import *
import pandas as pd
from math import *
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import ascii
import glob
import os
from pathlib import Path
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
import matplotlib.pyplot as plt
from astropy.nddata import Cutout2D

def des_flux(bcg,ra, dec):

    des_catalog_fluxes = pd.DataFrame(columns = ['des_g', 'des_r', 'des_i', 'des_z', 'des_Y', 'des_g_err','des_r_err', 'des_i_err', 'des_z_err', 'des_Y_err', 'id'])
    des_df = pd.read_csv("des_catalog.csv")
    cols = ['ra', 'dec', 'mag_auto_g', 'mag_auto_r', 'mag_auto_i', 'mag_auto_z','mag_auto_y','magerr_auto_g', 'magerr_auto_r', 'magerr_auto_i','magerr_auto_z', 'magerr_auto_y']
    des_fluxes = pd.DataFrame(columns = cols)
    def mag_to_mJy(mag):
        return 10**((8.90 - mag)/2.5)*1000
    mag_cols = ['mag_auto_g', 'mag_auto_r', 'mag_auto_i', 'mag_auto_z','mag_auto_y']
    magerr_cols = ['magerr_auto_g', 'magerr_auto_r', 'magerr_auto_i','magerr_auto_z', 'magerr_auto_y']
    c = SkyCoord(ra = ra*u.degree, dec = dec*u.degree, frame = 'icrs')
    des_catalog = SkyCoord(ra=np.array(des_df['ra'])*u.degree, dec=np.array(des_df['dec'])*u.degree, frame = 'icrs')
    idx_all = c.separation(des_catalog) < 0.004*u.deg
    idx, d2d, d3d = c.match_to_catalog_sky(des_catalog)
    match = pd.DataFrame(des_df.iloc[idx]).T.reset_index(drop = True)
    nearest_points = pd.DataFrame(des_df.iloc[idx_all]).reset_index(drop = True)
    match_jy = match[['ra', 'dec']].copy()
    for err_col,col in zip(magerr_cols, mag_cols):
        match_jy[col] = mag_to_mJy(match[col])
        match_jy[err_col] = match[err_col]*match_jy[col]
    match_jy['id'] = bcg
    des_fluxes = pd.concat([des_fluxes, match_jy])
    des_fluxes = des_fluxes.rename(columns = {"mag_auto_g":"des_g", "magerr_auto_g": "des_g_err", "mag_auto_r": "des_r", "magerr_auto_r": "des_r_err",
                     "mag_auto_i": "des_i", "magerr_auto_i":"des_i_err", "mag_auto_z":"des_z", "magerr_auto_z":"des_z_err",
                     "mag_auto_y":"des_Y", "magerr_auto_y": "des_Y_err"})
#    with open(bcg + "/" + bcg+'_allflux.reg','a') as f:
#        region = "# Region file format: DS9 version 4.0 \nglobal color=blue font='helvetica 10 normal roman' edit=1 move=1 delete=1 highlite=1 include=1 wcs=wcs \n" + "fk5; cross point " + str(ra) + " " + str(dec)
#        f.write(region)
    header = "# Region file format: DS9 version 4.0\nglobal color=blue font='helvetica 10 normal roman' edit=1 move=1 delete=1 highlite=1 include=1 wcs=wcs"
    f = open(bcg+'_allflux.reg','a')
    f.write(header)
    for ind,row in nearest_points.iterrows():
        region = "\nfk5; cross point " + str(row['ra']) + " " + str(row['dec'])
        f.write(region)
    des_flux = des_fluxes[['des_g', 'des_r', 'des_i', 'des_z', 'des_Y', 'des_g_err','des_r_err', 'des_i_err', 'des_z_err', 'des_Y_err', 'id']].copy()
    des_catalog_fluxes = des_catalog_fluxes.append(des_flux, ignore_index = True)
#    des_catalog_fluxes.to_csv(bcg + "_desfuxes.csv")

    return des_catalog_fluxes

def irac_flux(bcg, ra, dec):

    irac_catalog_fluxes = pd.DataFrame(columns = ['id','spitzer.irac.ch1',
       'spitzer.irac.ch1_err', 'spitzer.irac.ch2', 'spitzer.irac.ch2_err',
       'spitzer.irac.ch3', 'spitzer.irac.ch3_err', 'spitzer.irac.ch4',
       'spitzer.irac.ch4_err'])
    if bcg.startswith("ES1"):
        irac_es1_coord = pd.read_csv("es1_irac.csv", usecols = ['ra', 'dec'])
        irac_es1_fluxes = pd.read_csv("es1_irac.csv", usecols = ['flux_kr_36', 'uncf_kr_36','flux_kr_45', 'uncf_kr_45', 'flux_kr_58','uncf_kr_58', 'flux_kr_80', 'uncf_kr_80']).fillna(0)
        irac_es1_fluxes = irac_es1_fluxes * 0.001
        irac_catalog = pd.concat([irac_es1_coord, irac_es1_fluxes], axis = 1)
    if bcg.startswith("XMM"):
        irac_xmm_coord = pd.read_csv("xmm_irac.csv", usecols = ['ra', 'dec'])
        irac_xmm_fluxes = pd.read_csv("xmm_irac.csv", usecols = ['flux_kr_36', 'uncf_kr_36','flux_kr_45', 'uncf_kr_45', 'flux_kr_58','uncf_kr_58', 'flux_kr_80', 'uncf_kr_80']).fillna(0)
        irac_xmm_fluxes = irac_xmm_fluxes * 0.001
        irac_catalog = pd.concat([irac_xmm_coord, irac_xmm_fluxes], axis = 1)
    if bcg.startswith("CDFS"):
        irac_cdfs_coord = pd.read_csv("cdfs_irac.csv", usecols = ['ra', 'dec'])
        irac_cdfs_fluxes = pd.read_csv("cdfs_irac.csv", usecols = ['flux_kr_36', 'uncf_kr_36','flux_kr_45', 'uncf_kr_45', 'flux_kr_58','uncf_kr_58', 'flux_kr_80', 'uncf_kr_80']).fillna(0)
        irac_cdfs_fluxes = irac_cdfs_fluxes * 0.001
        irac_catalog = pd.concat([irac_cdfs_coord, irac_cdfs_fluxes], axis = 1)
    c = SkyCoord(ra = ra*u.degree, dec = dec*u.degree, frame = 'icrs')
    irac_catalog_coords = SkyCoord(ra=np.array(irac_catalog['ra'])*u.degree, dec=np.array(irac_catalog['dec'])*u.degree, frame = 'icrs')
    idx_all = c.separation(irac_catalog_coords) < 0.004*u.deg
    idx, d2d, d3d = c.match_to_catalog_sky(irac_catalog_coords)
    match_irac = pd.DataFrame(irac_catalog.iloc[idx]).T.reset_index(drop = True)
    nearest_points_irac = pd.DataFrame(irac_catalog.iloc[idx_all]).reset_index(drop = True)
    match_irac['id'] = bcg
    match_irac = match_irac.rename(columns = {"flux_kr_36":"spitzer.irac.ch1", "uncf_kr_36": "spitzer.irac.ch1_err","flux_kr_45": "spitzer.irac.ch2", 'uncf_kr_45': "spitzer.irac.ch2_err", 'flux_kr_58': "spitzer.irac.ch3", 'uncf_kr_58': "spitzer.irac.ch3_err", 'flux_kr_80': "spitzer.irac.ch4", 'uncf_kr_80': "spitzer.irac.ch4_err"})
#    with open(bcg + "/" + bcg+'_allflux.reg','a') as f:
#        region = "\nglobal color=red font='helvetica 10 normal roman' edit=1 move=1 delete=1 highlite=1 include=1 wcs=wcs \n" + "fk5; cross point " + str(ra) + " " + str(dec)
#        f.write(region)
    header = "\nglobal color=red font='helvetica 10 normal roman' edit=1 move=1 delete=1 highlite=1 include=1 wcs=wcs"
    f = open(bcg+'_allflux.reg','a')
    f.write(header)
    for ind,row in nearest_points_irac.iterrows():
        region = "\nfk5; cross point " + str(row['ra']) + " " + str(row['dec'])
        f.write(region)

    match_irac = match_irac[['id', 'spitzer.irac.ch1', 'spitzer.irac.ch1_err', 'spitzer.irac.ch2',
       'spitzer.irac.ch2_err', 'spitzer.irac.ch3', 'spitzer.irac.ch3_err',
       'spitzer.irac.ch4', 'spitzer.irac.ch4_err']]
    irac_catalog_fluxes = pd.concat([irac_catalog_fluxes, match_irac])

    return irac_catalog_fluxes


def hermes_flux(bcg, ra, dec):
    hermes_catalog = pd.DataFrame()
    hermes_path = "/Users/arames52/Desktop/Primary Project/Multiwavelength_bcg/Herschel/catalog/"
    c = SkyCoord(ra = ra*u.degree, dec = dec*u.degree, frame = 'icrs')
    if bcg.startswith("ES1"):
        es1_f = fits.open(hermes_path + "L6-ELAIS-S1-SWIRE_xID24_DR3.fits")
        es1_hermes = Table(es1_f[1].data).to_pandas()
        es1_h_catalog = SkyCoord(ra=np.array(es1_hermes['RA'])*u.degree, dec=np.array(es1_hermes['Dec'])*u.degree, frame = 'icrs')
        idx, d2d, d3d = c.match_to_catalog_sky(es1_h_catalog)
        es1_h_match =  pd.DataFrame(es1_hermes.iloc[idx]).T.reset_index(drop = True)
        es1_h_match['id'] = [bcg]
        hermes_catalog = pd.concat([hermes_catalog, es1_h_match])
    if bcg.startswith("XMM"):
        xmm_f = fits.open(hermes_path +"L6-XMM-LSS-SWIRE_xID24_DR3.fits")
        xmm_hermes = Table(xmm_f[1].data).to_pandas()
        xmm_h_catalog = SkyCoord(ra=np.array(xmm_hermes['RA'])*u.degree, dec=np.array(xmm_hermes['Dec'])*u.degree, frame = 'icrs')
        idx, d2d, d3d = c.match_to_catalog_sky(xmm_h_catalog)
        xmm_h_match =  pd.DataFrame(xmm_hermes.iloc[idx]).T.reset_index(drop = True)
        xmm_h_match['id'] = [bcg]
        hermes_catalog = pd.concat([hermes_catalog, xmm_h_match])
    if bcg.startswith("CDFS"):
        cdfs_f = fits.open(hermes_path +"L5-CDFS-SWIRE_xID24_DR3.fits")
        cdfs_hermes = Table(cdfs_f[1].data).to_pandas()
        cdfs_h_catalog = SkyCoord(ra=np.array(cdfs_hermes['RA'])*u.degree, dec=np.array(cdfs_hermes['Dec'])*u.degree, frame = 'icrs')
        idx, d2d, d3d = c.match_to_catalog_sky(cdfs_h_catalog)
        cdfs_h_match =  pd.DataFrame(cdfs_hermes.iloc[idx]).T.reset_index(drop = True)
        cdfs_h_match['id'] = [bcg]
        hermes_catalog = pd.concat([hermes_catalog, cdfs_h_match])

    hermes_catalog = hermes_catalog.reset_index(drop = True)
    for ind,row in hermes_catalog.iterrows():
        with open(row['id'] + "_allflux.reg","a") as f:
            region = "\nglobal color=green font='helvetica 10 normal roman' edit=1 move=1 delete=1 highlite=1 include=1 wcs=wcs \nfk5; cross point " + str(row['RA']) + " " + str(row['Dec'])
            bcg_coordinate = "\nglobal color=cyan font='helvetica 10 normal roman' edit=1 move=1 delete=1 highlite=1 include=1 wcs=wcs \nfk5; cross point " + str(ra) + " " + str(dec)
            f.write(region)
            f.write(bcg_coordinate)

    hermes_catalog['id'] == [bcg]
    hermes_catalog = hermes_catalog[['F24', 'e_F24', 'F250', 'et_F250', 'F350', 'et_F350', 'F500', 'et_F500', 'id']]
    hermes_catalog = hermes_catalog.rename(columns = {"F24":"spitzer.mips.24", "e_F24": "spitzer.mips.24_err", "F250": "herschel.spire.PSW", "et_F250": "herschel.spire.PSW_err","F350": "herschel.spire.PMW", "et_F350":"herschel.spire.PMW_err", "F500": "herschel.spire.PLW", "et_F500":"herschel.spire.PLW_err"})

    return hermes_catalog

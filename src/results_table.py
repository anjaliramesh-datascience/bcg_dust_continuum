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
import matplotlib.pyplot as plt
import matplotlib
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
import sys
sys.path.append("/Users/arames52/bcg_dust_continuum/src/")
import dust_continuum_morphology as dcm


non_detections = ['ES1-12', 'ES1-26', 'XMM-19' , 'XMM-27']
weak_detections = ['CDFS19', 'ES1_z_0.88', 'ES1-35', 'XMM_z_0.81']
lensed_detections = ['XMM-11', 'ES1-34']
good_detections = ["CDFS-18", "ES1-18", "ES1-25", "ES1_z_0.99","ES1_z_0.99b","ES1_z_1.04","ES1_z_1.38","ES1_z_1.40",
"ES1_z_1.60", "ES1_z_1.65", "ES1_z_1.70", "XMM-113", "XMM-29", "XMM-30", "XMM_z_0.9", "XMM_z_1.0"]
detections = good_detections + lensed_detections
bad_detections = non_detections + weak_detections
all_bcgs = good_detections + lensed_detections + weak_detections + non_detections
cigale_results_path = "/Users/arames52/Downloads/cigale-v2020.0/pcigale/data/"
bcg_redshift_df = pd.read_excel("/Users/arames52/bcg_dust_continuum/notebook/data/Derived_Data/BCG_redshifts.xlsx")

def sfr_mstar_results(bcg_name):

    results_file = glob.glob(cigale_results_path + bcg_name + "/out/results.fits")[0]
    data = Table(fits.open(results_file)[1].data)
    df = data[['id','bayes.sfh.sfr', 'bayes.sfh.sfr_err', 'bayes.stellar.m_star', 'bayes.stellar.m_star_err', 'best.reduced_chi_square']].to_pandas()
    df = df.rename(columns = {'bayes.sfh.sfr':"SFR", 'bayes.sfh.sfr_err':"SFR_err", 'bayes.stellar.m_star':"Stellar_mass", 'bayes.stellar.m_star_err':"Stellar_mass_err", "best.reduced_chi_square": "reduced_chi_square"})
    df['sSFR'] = df['SFR']/df['Stellar_mass']
    df['sSFR_err'] = (df['SFR']/df['Stellar_mass'])*np.sqrt((df['SFR_err']/df['SFR'])**2 + (df['Stellar_mass_err']/df['Stellar_mass'])**2)
    df['sSFR_gyr'] = np.array(df['sSFR'])*1e9
    df['sSFR_err_gyr'] = np.array(df['sSFR_err'])*1e9
    df['redshift'] = bcg_redshift_df[bcg_redshift_df['bcg'] == bcg_name]['redshift'].values[0]
    df['spec_z'] = bcg_redshift_df[bcg_redshift_df['bcg'] == bcg_name]['spec_z'].values[0]
    df['Age'] = cosmo.age(df['redshift']).value
    df['sSFR_MS'] = 26 * np.array(df['Age'])**(-2.2)
    
    return df

def sersic_results(bcg_name):

    morph = dcm.compute_sersic(bcg_name)
    redshift = bcg_redshift_df[bcg_redshift_df['bcg'] == bcg_name]['redshift'].values[0]
    re_kpc = cosmo.kpc_proper_per_arcmin(redshift)*(morph.sersic_rhalf*0.045 /60)
    sersic_properties = [{"n": round(morph.sersic_n,1), "re": round(re_kpc.value,1)}]
    df = pd.DataFrame(sersic_properties)

    return df

def main_table():

    master_df = pd.DataFrame()

    for bcg in all_bcgs:
        sf_df = sfr_mstar_results(bcg)
        if bcg in bad_detections:
            sersic_df = pd.DataFrame([{"n":None, "re": None}])
        else:
            sersic_df = sersic_results(bcg)
        
        df = pd.concat([sf_df, sersic_df], axis = 1)
        master_df = pd.concat([master_df, df], axis = 0)

    return master_df.reset_index(drop = True)

# results_table = main_table()
# bcg_type = []
# for ind,row in results_table.iterrows():
#     if row['sSFR_gyr'] < row['sSFR_MS']/2:
#         bcg_type.append(0)
#     elif (row['sSFR_gyr'] > row['sSFR_MS']/2) & (row['sSFR_gyr'] < row['sSFR_MS']*2):
#         bcg_type.append(1)
#     elif (row['sSFR_gyr'] > row['sSFR_MS']*2):
#         bcg_type.append(2)
# results_table['bcg_type'] = bcg_type
# results_table.to_csv("/Users/arames52/bcg_dust_continuum/notebook/data/Derived_Data/results_table.csv", index = False)
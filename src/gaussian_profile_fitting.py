from astropy.io import fits
import sys
sys.path.append('/Users/arames52/bcg_dust_continuum/src/')
import glob
import pandas as pd
import pickle

# Setting path of files and folder names
casa_imaging_path = "/Users/arames52/bcg_dust_continuum/notebook/data/CASA_imaging/"
casa_imaging_folders = ["briggs_imaging/", "natural_imaging/", "tapered_imaging/"]
tapered_imaging_folders = ["0.7_arcsec/", "1_arcsec/", "2_arcsec/"]

# BCG region parameters file
region_file_path = "/Users/arames52/bcg_dust_continuum/notebook/data/Derived_Data/bcg_regions.txt"
region_df = pd.read_csv(region_file_path, sep = ",")
bcgs = list(region_df['id'])

def region_params(ra, dec, radius):
    x, y = str(ra) + "deg", str(dec) + "deg"
    r =  str(radius) + "arcsec"
    region = "circle[["+ x +","+y+ "]," + r+ "]"
    return region

def imfit_routine(bcg_name:str):
    imfit_results = {}
    casa_images = {}
    # 5 CASA images produced - Briggs, natural weighting and tapered images with natural weighting
    casa_images[0] = glob.glob(casa_imaging_path + casa_imaging_folders[0] + bcg_name + "*.fits")[0]
    casa_images[1] = glob.glob(casa_imaging_path + casa_imaging_folders[1] + bcg_name + "*.fits")[0]
    casa_images[2] = glob.glob(casa_imaging_path + casa_imaging_folders[2] + tapered_imaging_folders[0]+ bcg_name + "*.fits")[0]
    casa_images[3] = glob.glob(casa_imaging_path + casa_imaging_folders[2] + tapered_imaging_folders[1]+ bcg_name + "*.fits")[0]
    casa_images[4] = glob.glob(casa_imaging_path + casa_imaging_folders[2] + tapered_imaging_folders[2]+ bcg_name + "*.fits")[0]

    for i in range(5):
        img = casa_images[i]
        ra = region_df[region_df['id']== bcg_name]['ra'].values[0]
        dec = region_df[region_df['id']== bcg_name]['dec'].values[0]
        if i == 3:
            radius = region_df[region_df['id']== bcg_name]['radius'].values[0] + 0.3
        elif i == 4:
            radius = region_df[region_df['id']== bcg_name]['radius'].values[0] + 0.6
        else:
            radius = region_df[region_df['id']== bcg_name]['radius'].values[0]
        
        imfit_region = region_params(ra,dec,radius)
        imfit_res = imfit(img, region = imfit_region)
        imfit_results[i] = imfit_res

    return imfit_results

def imfit_single_image(bcg_name:str):
    image = glob.glob(casa_imaging_path + casa_imaging_folders[1] + bcg_name + "*.fits")[0]
    ra = region_df[region_df['id']== bcg_name]['ra'].values[0]
    dec = region_df[region_df['id']== bcg_name]['dec'].values[0]
    radius = region_df[region_df['id']== bcg_name]['radius'].values[0]
    imfit_region = region_params(ra,dec,radius)
    imfit_res = imfit(image, region = imfit_region)
    return imfit_res

# Imfit on non-detections or low S/N detections not considered (Tapering does not improve the detections)
dont_care_about = ["ES1-12", "ES1-26", "ES1-35", "CDFS19", "XMM-19", "XMM-27"]

full_imfit_results = {}
for bcg in bcgs:
    if bcg in dont_care_about:
        continue
    elif bcg in ["ES1-34", "ES1_z_0.88", "XMM_z_0.81"]:
        full_imfit_results[bcg] = imfit_single_image(bcg)
    else:
        full_imfit_results[bcg] = imfit_routine(bcg)
    
with open("/Users/arames52/bcg_dust_continuum/notebook/data/Derived_Data/imfit_results.pkl", "wb") as f:
    pickle.dump(full_imfit_results, f)


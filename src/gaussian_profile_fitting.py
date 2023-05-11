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

def natural_imfit_flux(bcg_name:str):

    data_path = casa_imaging_path + casa_imaging_folders[1] + bcg_name + "_natural.fits"
    ra = region_df[region_df['id']== bcg_name]['ra'].values[0]
    dec = region_df[region_df['id']== bcg_name]['dec'].values[0]
    radius = region_df[region_df['id']== bcg_name]['radius'].values[0]
    region_string = region_params(ra, dec, radius)
    imfit_result_bcg = imfit(data_path, region = region_string)

    # with open(casa_imaging_path + bcg_name + "_natural_imfit_results.pkl", "wb") as f:
    #     pickle.dump(imfit_result_bcg, f)
    
    return imfit_result_bcg

def imfit_routine(bcg_name:str):
    imfit_results = {}
    casa_images = {}
    # 4 CASA images produced - Briggs, natural weighting and tapered images (0.7" and 1") with natural weighting
    casa_images[0] = glob.glob(casa_imaging_path + casa_imaging_folders[0] + bcg_name + "*.fits")[0]
    casa_images[1] = glob.glob(casa_imaging_path + casa_imaging_folders[1] + bcg_name + "*.fits")[0]
    casa_images[2] = glob.glob(casa_imaging_path + casa_imaging_folders[2] + tapered_imaging_folders[0]+ bcg_name + "*.fits")[0]
    casa_images[3] = glob.glob(casa_imaging_path + casa_imaging_folders[2] + tapered_imaging_folders[1]+ bcg_name + "*.fits")[0]
    ra = region_df[region_df['id']== bcg_name]['ra'].values[0]
    dec = region_df[region_df['id']== bcg_name]['dec'].values[0]
    r = region_df[region_df['id']== bcg_name]['radius'].values[0]
    radius = [r-0.3, r, r, r + 0.5]
    for i in range(4):
        imfit_region = region_params(ra,dec,radius[i])
        imfit_res = imfit(casa_images[i], region = imfit_region)
        imfit_results[i] = imfit_res

    return imfit_results

# Imfit on non-detections or low S/N detections not considered (Tapering does not improve the detections)
non_detections = ["ES1-12", "ES1-26", "ES1-35", "CDFS19", "XMM-19", "XMM-27"]
good_detections = ["CDFS-18", "ES1-18", "ES1-25", "ES1_z_0.99", "ES1_z_1.04", "ES1_z_1.38", "ES1_z_1.40", 
"ES1_z_1.60", "ES1_z_1.65", "ES1_z_1.70", "XMM-113", "XMM-11", "XMM-29", "XMM-30",
"XMM_z_0.9", "XMM_z_1.0"]
weak_detections = ["ES1-34", "ES1_z_0.88","ES1_z_0.99b", "XMM_z_0.81"]

# Running imfit on naturally weighted image and then all
imfit_results = {}
for bcg in good_detections + weak_detections:
    print("Natural Weighting IMFIT")
    print(bcg)
    imfit_results[bcg] = natural_imfit_flux(bcg)
all_imfit_results = {}
for bcg in good_detections:
    print("IMFIT on all")
    print(bcg)
    all_imfit_results[bcg] = imfit_routine(bcg)

with open(casa_imaging_path + "natural_imfit_results.pkl", "wb") as f:
    pickle.dump(imfit_results, f)
with open(casa_imaging_path + "all_imfit_results.pkl", "wb") as f:
    pickle.dump(all_imfit_results, f)


# Parsing through saved imfit results

def extract_imfit_results():

    flux_dict = {}
    flux_err_dict = {}
    beam_dict = {}

    with open("/Users/arames52/bcg_dust_continuum/notebook/data/CASA_imaging/all_imfit_results.pkl", "rb") as f:
        imfit_results = pickle.load(f)
    
    for i in range(4):
        flux_dict[i] = [val[i]['deconvolved']['component0']['flux']['value'][0] if val[i] != False else 0 for key, val in imfit_results.items()]
        flux_err_dict[i] = [val[i]['deconvolved']['component0']['flux']['error'][0] if val[i] != False else 0 for key, val in imfit_results.items()]
        beam_dict[i] = [round(val[i]['deconvolved']['component0']['beam']['beamarcsec']['major']['value'],1) if val[i] != False else 0 for key, val in imfit_results.items()]
    
    return flux_dict, flux_err_dict, beam_dict, list(imfit_results.keys())

def plot_imfit_results():

    fig, ax = plt.subplots(4,4, dpi = 300)
    ax = ax.ravel()
    flux, fluxerr, beam, bcgs = extract_imfit_results()
    for i in range(16):
        for j in range(4):
            ax[i].scatter(beam[j][i], flux[j][i] * 1000)
        ax[i].set_title(bcgs[i])

    return fig

# plt.style.use("science")
# fig = plot_imfit_results()
# fig.suptitle("Dust Continuum Flux Measurement", size = 20)
# fig.supxlabel("Beam major FWHM ['']", size = 20, y = 0.05)
# fig.supylabel("Measured Flux [mJy]", size = 20, x = 0.05)
# fig.savefig("/Users/arames52/bcg_dust_continuum/notebook/plots/imfit_plot.png", dpi = 300)
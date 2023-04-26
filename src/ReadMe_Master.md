# The Dust Continuum Morphology of Brightest Cluster Galaxies

## Project Description

In this project, we analyze the physical properties of objects (galaxies) in radio-wavelength images (fits data type) using statistical methods like 2D Gaussian modeling, Chi-square fitting, and object surface profile modeling. For a more detailed description of the science, refer to paper in prep. 

### Some preliminary steps and information

1. git clone https://github.com/anjali-ramesh-1995/bcg_dust_continuum.git
2. In the bcg_parameter_file.py file replace "/Users/arames52/" wherever found to folder path you cloned the repository to
3. Assuming conda is installed - conda activate bcg_project/
4. conda install -r requirements.txt
5. **src** folder has all the codes and readme files
6. **notebook** folder contains data and plots 
7. Activate the environment and run ...

KEYWORDS - 

## Data Description

### Primary Data 

Calibrated measurement sets obtained after executing *scriptforPI.py from the downloaded raw observations from ALMA Science Archive with project id 018.1.00828.S and 2019.1.01027.S (P.I. Noble). 

### Supporting Data

Multi-wavelength Flux Catalog - Provided the coordinates of the objects, some functions in data_ingestion.py help in querying catalog information and image cutouts from several surveys and telescope archive data like Dark Energy Survey (DES), the Spitzer Space Telescope, Herschel Space Observatory. For more information look at data_ingestion.py.

KEYWORDS - 

## Data Analysis

Note - helper_functions.py (hf), bcg_parameter_file.py (bpf), data_ingestion.py (di)   are imported into most scripts

### Data Reduction and Imaging with CASA

* Input - Calibrated measurement sets located at "/Users/arames52/bcg_dust_continuum/notebook/data/Measurement_Sets/"
* Code - data_reduction_and_imaging.py (drai)
* Functions
    1. get_tclean_parameters() {hf} - Return list of weighting and tapering options read from weighting_tapering_values.txt
    2. folder_creation() {hf}  - creates folders for output storage
    3. alma_imaging {drai} - given measurement set, and list of weights and tapering schemes - 
        * run TCLEAN with defined parameters in bpf
        * convert .pbcor image file to .fits for convenient usage in python
    
* Note - If auto-masking algorithm to be performed, change to interactive = True and usemask = "auto-multithresh"

### 2D Gaussian Profile Modeling with CASA

* Input - .fits image data produced from different imaging parameters
* Code - gaussian_profile_fitting.py (gpf)
* Functions - 
    1. region_params() - returns region that encloses the object

* Note - Run this in CASA terminal

### Image Properties - RMS, S/N - with CASA

* Input - .fits image data
* Code - image_properties.py
* Functions - 
    1. signal_rms() - runs CASA's IMSTAT to calculate properties in a given region

* Note - Run on CASA terminal

### Multi-wavelength catalog download and compilation

Download object information from different surveys and telescope and save as csv file 

data_ingestion.py

### Surface Brightness Profile Fitting with Statmorph

dust_continuum_morphology.py

* Functions - make_psf(), compute_sersic(), plot_profiles()


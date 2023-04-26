import pandas as pd
import numpy as np
import sys
sys.path.append("/Users/arames52/bcg_dust_continuum/src/")
import bcg_parameter_file as p
import os
import glob
from pathlib import Path

## TCLEAN ROUTINE

def get_weights_and_tapers():

    # Uploading TCLEAN weighting and tapering parameter file
    tclean_param_file = pd.read_csv(p.code_path + "weighting_tapering_values.txt", delimiter=',')
    tclean_param_file = tclean_param_file.astype(object).replace(np.nan, 'None')
    weights = list(tclean_param_file['weighting'])
    tapers = list(tclean_param_file['uvtaper'])

    return weights, tapers


def folder_creation(folder_names_list):

    # tclean output folder creation
    for folder_name in folder_names_list:
        path_string = p.data_path + folder_name + "/"
        Path(path_string).mkdir(parents=True, exist_ok=True)
    
    return None




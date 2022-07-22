# dommmap_tools.py 
# This file contains helper functions that maintain the ECOD domain dictionary

import os

import pathlib

import requests

from numpy import load as npload

from numpy import save as npsave

from DomainMapper.dommap_io import *


# file name variables
__ecod_domain_txt_url = 'http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt'

__ecod_domain_txt_fn = 'ecod.latest.domains.txt'

__ecod_domain_npy_fn = 'ecod.latest.domains.npy'

# install path variables
__pt = os.path.dirname(os.path.realpath(__file__))

__txt_fn_path = os.path.join(__pt,__ecod_domain_txt_fn)

__npy_fn_path = os.path.join(__pt,__ecod_domain_npy_fn)


# Checking if file is older than 2 months
def __out_of_date(file_path):
    return os.times().elapsed - os.path.getctime(file_path) > 2628000


# Checking if file path exsists
def __ecod_domain_exists(file_path):
    if os.path.exists(file_path):
        return True
    return False


# Using requests to download the ecod domain defintions
def __download_latest_ecod_domain_txt():
    ecod_domain_def = requests.get(__ecod_domain_txt_url, allow_redirects=False)

    with open(__txt_fn_path,"wb") as write_file:
        write_file.write(ecod_domain_def.content)


# Parsing ecod domain definitions
def __parse_ecod_domain_txt(file_path,save_path):
    ecod_domain_dict = dict()
    ecod_line_info_idx = [3, 9, 10, 11, 12, 13] # Necessary column indices from ecod domain definitions

    with open(file_path, 'r') as ecod_latest_file:
        for line in ecod_latest_file:
            if not line.startswith('#'):
                clean_line = [l.strip('"') for l in line.split('\t')]
                f_id, arch, x_group, h_group, t_group, f_group = [clean_line[i] for i in ecod_line_info_idx]
                if f_group not in ecod_domain_dict.keys():
                    if x_group == 'NO_X_NAME':
                        if h_group != 'NO_H_NAME':
                            x_group = h_group
                        else:
                            x_group = t_group
                    ecod_domain_dict[f_group] = [f_id, arch, x_group, t_group]

    npsave(save_path, ecod_domain_dict)
    os.remove(file_path)


# Update workflow
def update():
    notice_msg("Updating 'ecod.latest.domains' from http://prodata.swmed.edu/ecod")
    __download_latest_ecod_domain_txt()
    __parse_ecod_domain_txt(__txt_fn_path,__npy_fn_path)


# Loading parsed ecod domain defintions and return as a dictionary
def load(file_path = False):
    # if user provides their own ecod domains, this part will parse their file
    if file_path:
        if __ecod_domain_exists(file_path):
            # if the domain definitions are out of date, a warning will be presented but the file will be parsed and read
            if __out_of_date(file_path):
                warning_msg("WARNING: '{}' is out of date. Please update your domain definitions file.".format(pathlib.Path(file_path).name))

            __parse_ecod_domain_txt(file_path,__ecod_domain_npy_fn)
            return npload(__ecod_domain_npy_fn, allow_pickle = True).item()
            
        else:
            error_msg("ERROR: '{}' could not be found.".format(pathlib.Path(file_path).name))
    
    # default is to parse the domain definitions provided in this module
    else:
        # if the domain definitions are not present, a warning will be presented and the file will be updated, the  the file will be read
        if __ecod_domain_exists(__npy_fn_path):
            # if the domain definitions are out of date, a warning will be presented but the file will be read
            if __out_of_date(__npy_fn_path):
                warning_msg("'ecod.latest.domains' is out of date. Please re-run program with the --update flag.")
            return npload(__npy_fn_path, allow_pickle = True).item()

        else:
            notice_msg("'ecod.latest.domains' not found. Proceeding with download ...")
            update()
            return npload(__npy_fn_path, allow_pickle = True).item()


if __name__ == "__main__":
    ecod = load()
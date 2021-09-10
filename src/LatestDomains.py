import os
import numpy as np
from wget import download

ecod_domain_txt_fn = 'ecod.latest.domains.txt'
ecod_domain_npy_fn = 'ecod.latest.domains.npy'

# TODO - check_for_ecod_json
# 1. Check a data directory that would come as part of the package
def __check_for_ecod_domain_file():
    if os.path.exists(ecod_domain_npy_fn) and os.times().elapsed - os.path.getctime(ecod_domain_npy_fn) > 2628000:
        print("ecod.latest.domains.npy already exisits")
        return False
    return True

def __download_latest_ecod_domain_txt():
    if not os.path.exists(ecod_domain_txt_fn):
        print("Downloading 'ecod.latest.domains.txt' from ECOD website")
        ecod_domain_url = 'http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt'
        download(ecod_domain_url)

def __parse_ecod_domain_txt(ecod_latest_fn):
    ecod_domain_dict = dict()
    ecod_line_info_idx = [3, 9, 10, 11, 12, 13]
    with open(ecod_latest_fn, 'r') as ecod_latest_file:
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

    np.save(ecod_domain_npy_fn ,ecod_domain_dict)

    os.remove(ecod_latest_fn)

def parse_user_ecod_domain_txt(ecod_latest_fn):
    if os.times().elapsed - os.path.getctime(ecod_latest_fn) > 2628000:
        print('\'{}\'is out of date. Consider fetching newest version at http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt'.format(ecod_latest_fn))
    ecod_domain_dict = dict()
    ecod_line_info_idx = [3, 9, 10, 11, 12, 13]
    with open(ecod_latest_fn, 'r') as ecod_latest_file:
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

    return ecod_domain_dict

def load():
     return np.load(ecod_domain_npy_fn, allow_pickle = True).item()

# TODO - fetch()
# 1. Check to see if an exisiting 'ecod.latest.domains.npy' is out of date or not before downlaoding/parsing
def fetch():
    if __check_for_ecod_domain_file():
        __download_latest_ecod_domain_txt()
        __parse_ecod_domain_txt(ecod_domain_txt_fn)

if __name__ == "__main__":
    fetch()
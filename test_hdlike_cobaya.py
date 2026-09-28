"""Test the initialization of the CMB-HD likelihood, and evaluate it at
the fiducial cosmology.

"""

import os
import argparse
import warnings
import numpy as np
import hdlike
from hd_mock_data import hd_data
from cobaya.yaml import yaml_load_file
from cobaya.run import run


# command-line args:
formatter_class = argparse.ArgumentDefaultsHelpFormatter
description = ("Test the initialization of `HDLike` with Cobaya, and "
               "evaluate the likelihood at the fiducial cosmology.")
parser = argparse.ArgumentParser(formatter_class=formatter_class,
                                 description=description)
parser.add_argument('--lensed',  action='store_true',
                    help="Use lensed CMB data (delensed data is used by default).")
parser.add_argument('--desi', action='store_true',
                    help="Include DESI BAO data (it is not included by default).")
parser.add_argument('--feedback', action='store_true',
                    help=("Use a model with baryonic feedback effects "
                          "(a CDM-only model is used by default)."))
parser.add_argument('--use_class', action='store_true',
                    help=("Use CLASS instead of CAMB (CAMB is used by default)."
                          " NOTE that (1) CLASS does not calculate delensed"
                          " power spectra, so you must also pass `--lensed`,"
                          " and (2) CLASS must be modified to use the `hdlike`"
                          " defaults; see the `hdlike` 'readme' file."))
args = parser.parse_args()

if args.use_class and (not args.lensed):
    raise ValueError("Cannot calculate delensed CMB power spectra with CLASS. "
                     "You must either pass `--lensed` to use lensed power "
                     "spectra, or use CAMB to use delensed power spectra.")

# tell the user what they asked for
cmb_type = 'lensed' if args.lensed else 'delensed'
print(f'Using {cmb_type} CMB-HD data.')
if args.desi:
    print('Using DESI BAO data.')
if args.feedback:
    print('Using a model with baryonic feedback.')
if args.use_class:
    print("Using CLASS.")
print('-----\n')

# read in the YAML file
yaml_file = 'example_cmbhd_CLASS.yaml' if args.use_class else 'example_cmbhd_CAMB.yaml'
info = yaml_load_file(yaml_file)
info['output'] = None # don't save the output

# get the path to the hdlike files:
data_dir = os.path.join(os.path.dirname(hdlike.__file__), 'data/')
data_path = lambda x: os.path.join(data_dir, x)
# set which kind of data is being used:
info['likelihood']['hdlike.hdlike.HDLike'] = {'delensed': not args.lensed,
                                              'baryonic_feedback': args.feedback}
if args.desi:
    desi_data_fname = 'mock_desi_bao_rs_over_DV_data.txt'
    if args.use_class:
        desi_data_fname = f'class_{desi_data_fname}'
    desi_data_file = data_path(desi_data_fname)
    desi_cov_file = data_path('mock_desi_bao_rs_over_DV_cov.txt')
    info['likelihood']['bao.generic'] = {'measurements_file': desi_data_file,
                                         'cov_file': desi_cov_file}

# set the correct non-linear model:
if args.feedback:
    if args.use_class:
        info['theory']['classy']['extra_args']['hmcode_version'] = '2020_baryonic_feedback'
    else:
        info['theory']['camb']['extra_args']['halofit_version'] = 'mead2020_feedback'

# use the `evaluate` sampler at a fixed cosmology:
hd_datalib = hd_data.HDMockData()
if args.use_class:
    all_params = hd_datalib.class_settings()
    # get CAMB params for correct values of Neff, mnu
    camb_params = hd_datalib.camb_settings()
    all_params['N_eff'] = camb_params['nnu']
    all_params['mnu'] = camb_params['mnu']
    param_names = ['theta_s_100', 'omega_b', 'omega_cdm', 'ln_A_s_1e10', 'n_s', 'tau_reio', 'N_eff', 'mnu']
else:
    all_params = hd_datalib.camb_settings()
    all_params['theta_MC_100'] = 100 * all_params['cosmomc_theta']
    all_params['logA'] = np.log(1e10 * all_params['As'])
    param_names = ['theta_MC_100', 'ombh2', 'omch2', 'logA', 'ns', 'tau', 'nnu', 'mnu']
fiducial_params = {param: all_params[param] for param in param_names}
info['sampler'] = {'evaluate': {'override': fiducial_params.copy()}}


# run cobaya
updated_info, sampler = run(info, no_mpi=True)

# get the chi^2 values
products = sampler.products()
sample = products['sample'].data.iloc[0]
chi2_hd = sample['chi2__hdlike.hdlike.HDLike']
if args.desi:
    chi2_desi = sample['chi2__bao.generic']
else:
    chi2_desi = 0
chi2_tot = chi2_hd + chi2_desi

# compare with the expected results:
if args.use_class:
    # TODO: update these
    expected_chi2_desi = 2.60556e-9
    expected_chi2_hd_values = {'lensed': 3.91772e-10,
                               'lensed_feedback': 1.5083e-22}
else:
    expected_chi2_desi = 0.00452987
    expected_chi2_hd_values = {'lensed': 2.69758,
                               'delensed': 3.09648,
                               'lensed_feedback': 2.62660,
                               'delensed_feedback': 3.05534}
chi2key = f'{cmb_type}_feedback' if args.feedback else cmb_type
expected_chi2_hd = expected_chi2_hd_values[chi2key]
expected_chi2_tot = expected_chi2_hd
if args.desi:
    expected_chi2_tot += expected_chi2_desi

chi2_hd_diff = expected_chi2_hd - chi2_hd
if args.desi:
    chi2_desi_diff = expected_chi2_desi - chi2_desi
chi2_tot_diff = expected_chi2_tot - chi2_tot

print('\n-----')
print(f'chi^2 for CMB-HD = {chi2_hd}; expected chi^2 = {expected_chi2_hd} (difference = {chi2_hd_diff})')
if args.desi:
    print(f'chi^2 for DESI = {chi2_desi}; expected chi^2 = {expected_chi2_desi} (difference = {chi2_desi_diff})')
    print(f'Total chi^2 for CMB-HD + DESI = {chi2_tot}; expected chi^2 = {expected_chi2_tot} (difference = {chi2_tot_diff})')

if np.isclose(chi2_tot, expected_chi2_tot, atol=0.05):
    print('Success! The chi^2 value from Cobaya matches the expected value.')
else:
    print('Test failed: the chi^2 value from Cobaya does not match the expected value.')
    print('Note: this may occur if you have modified the `example_cmbhd.yaml` file.')




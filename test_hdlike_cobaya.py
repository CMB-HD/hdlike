"""tests the initialization of the CMB-HD likelihood and evaluate it at the fiducial cosmology"""
import os
import argparse
import warnings
import numpy as np
import hdlike
from cobaya.yaml import yaml_load_file
from cobaya.run import run


parser = argparse.ArgumentParser(description='Test the initialization of HDlike with Cobaya, and evaluate the likelihood at the fiducial cosmology.')

parser.add_argument('--lensed',  action='store_true', help='Use lensed CMB data (delensed data is used by default).')
parser.add_argument('--desi', action='store_true', help='Include DESI BAO data (it is not included by default).')
parser.add_argument('--feedback', action='store_true', help='Use a model with baryonic feedback effects (a CDM-only model is used by default).')
parser.add_argument('--class', dest='use_class', action='store_true', help='Calculate the theory with CLASS instead of CAMB, using the settings in `example_cmbhd_class_lcdm_nrun_nnu_mnu.yaml`, and compare it to the CMB-HD bandpowers that were calculated with CLASS. Requires `--lensed`, and cannot be used with `--feedback`.')


args = parser.parse_args()

# tell the user what they asked for
cmb_type = 'lensed' if args.lensed else 'delensed'
print(f'Using {cmb_type} CMB-HD data.')
if args.desi:
    print('Using DESI BAO data.')
if args.feedback:
    print('Using a model with baryonic feedback.')
if args.use_class:
    if not args.lensed:
        raise ValueError("Cannot compute delensed power spectra using CLASS; pass `--lensed` with `--class`.")
    if args.feedback:
        raise ValueError("Cannot use the HMCode2020 baryonic feedback model with CLASS.")
    print('Calculating the theory with CLASS.')
print('-----\n')

# get the path to the hdlike files
data_dir = os.path.join(os.path.dirname(hdlike.__file__), 'data/')
data_path = lambda x: os.path.join(data_dir, x)

# read in the YAML file
yaml_fname = 'example_cmbhd_class_lcdm_nrun_nnu_mnu.yaml' if args.use_class else 'example_cmbhd.yaml'
info = yaml_load_file(yaml_fname)

if args.use_class:
    # the CLASS example file has the `/PATH/TO/` placeholder in the file names
    #  for the BBN table, the DESI BAO data, and the proposal matrix;
    #  replace it with the absolute path to the `hdlike/data` directory
    def replace_paths(x):
        if isinstance(x, dict):
            return {key: replace_paths(value) for key, value in x.items()}
        elif isinstance(x, str) and ('/PATH/TO/hdlike/hdlike/data/' in x):
            return x.replace('/PATH/TO/hdlike/hdlike/data/', data_dir)
        else:
            return x
    info = replace_paths(info)
    # the DESI BAO data is included in the example file; only use it if requested
    if 'bao.generic' in info['likelihood'].keys():
        info['likelihood'].pop('bao.generic')

# set what kind of data we're using
info['likelihood']['hdlike.hdlike.HDLike'] = {'delensed': not args.lensed,
        'baryonic_feedback': args.feedback}
if args.use_class:
    # the CMB-HD bandpowers calculated with CLASS are only available for v1.2 of the data
    info['likelihood']['hdlike.hdlike.HDLike']['use_class'] = True
    info['likelihood']['hdlike.hdlike.HDLike']['hd_data_version'] = 'v1.2'
if args.desi:
    desi_data_file = data_path('mock_desi_bao_rs_over_DV_data.txt')
    desi_cov_file = data_path('mock_desi_bao_rs_over_DV_cov.txt')
    info['likelihood']['bao.generic'] = {'measurements_file': desi_data_file, 'cov_file': desi_cov_file}
if args.feedback:
    info['theory']['camb']['extra_args']['halofit_version'] = 'mead2020_feedback'

# use the `evaluate` sampler at a fixed cosmology:
if args.use_class:
    # the same cosmology as below, in terms of the CLASS parameters that are
    #  sampled in the example file. NOTE that H0 is used in place of theta_MC,
    #  and that there are three massive neutrinos with N_eff = 3.044, to match
    #  the CMB-HD bandpowers that were calculated with CLASS.
    fiducial_params = {'ln_A_s_1e10': 3.044,
                       'n_s': 0.9649,
                       'H0': 67.36,
                       'omega_b': 0.02237,
                       'omega_cdm': 0.1200,
                       'mnu': 0.06,
                       'N_eff': 3.044,
                       'tau_reio': 0.0544,
                       'alpha_s': 0.0}
else:
    fiducial_params = {'logA': 3.044,
                       'ns': 0.9649,
                       'theta_MC_100': 1.04092,
                       'ombh2': 0.02237,
                       'omch2': 0.1200,
                       'mnu': 0.06,
                       'nnu': 3.046,
                       'tau': 0.0544}
info['sampler'] = {'evaluate': {'override': fiducial_params.copy()}}

# don't save the output
info['output'] = None


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

# expected results
expected_chi2_desi_values = {'camb': 8.38174e-06,
        'class': 7.36069e-04}
expected_chi2_hd_values = {'lensed': 6.21009,
        'delensed': 3.41060,
        'lensed_feedback': 6.17882,
        'delensed_feedback': 3.39521,
        'lensed_class': 2.40332e-03}
if args.use_class:
    expected_chi2_hd = expected_chi2_hd_values[f'{cmb_type}_class']
    expected_chi2_desi = expected_chi2_desi_values['class']
elif args.feedback:
    expected_chi2_hd = expected_chi2_hd_values[f'{cmb_type}_feedback']
    expected_chi2_desi = expected_chi2_desi_values['camb']
else:
    expected_chi2_hd = expected_chi2_hd_values[cmb_type]
    expected_chi2_desi = expected_chi2_desi_values['camb']

print('\n-----')
expected_chi2_tot = expected_chi2_hd
if args.desi:
    expected_chi2_tot += expected_chi2_desi

chi2_hd_diff = expected_chi2_hd - chi2_hd
if args.desi:
    chi2_desi_diff = expected_chi2_desi - chi2_desi
chi2_tot_diff = expected_chi2_tot - chi2_tot

print(f'chi^2 for CMB-HD = {chi2_hd}; expected chi^2 = {expected_chi2_hd} (difference = {chi2_hd_diff})')
if args.desi:
    print(f'chi^2 for DESI = {chi2_desi}; expected chi^2 = {expected_chi2_desi} (difference = {chi2_desi_diff})')
    print(f'Total chi^2 for CMB-HD + DESI = {chi2_tot}; expected chi^2 = {expected_chi2_tot} (difference = {chi2_tot_diff})')

if np.isclose(chi2_tot, expected_chi2_tot, atol=1e-3):
    print('Success! The chi^2 value from Cobaya matches the expected value.')
else:
    print('Test failed: the chi^2 value from Cobaya does not match the expected value.')
    print(f'Note: this may occur if you have modified the `{yaml_fname}` file.')




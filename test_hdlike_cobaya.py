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


args = parser.parse_args()

# tell the user what they asked for
cmb_type = 'lensed' if args.lensed else 'delensed'
print(f'Using {cmb_type} CMB-HD data.')
if args.desi:
    print('Using DESI BAO data.')
if args.feedback:
    print('Using a model with baryonic feedback.')
print('-----\n')

# read in the YAML file
info = yaml_load_file('example_cmbhd.yaml')

if 'classy' in info['theory'].keys():
    if not args.lensed:
        raise ValueError("Cannot compute delensed power spectra using CLASS.")
    if args.feedback:
        raise ValueError("Cannot use the HMCode2020 baryonic feedback model with CLASS.")

# get the path to the hdlike files
data_dir = os.path.join(os.path.dirname(hdlike.__file__), 'data/')
data_path = lambda x: os.path.join(data_dir, x)

# set what kind of data we're using
info['likelihood']['hdlike.hdlike.HDLike'] = {'delensed': not args.lensed,
        'baryonic_feedback': args.feedback}
if args.desi:
    desi_data_file = data_path('mock_desi_bao_rs_over_DV_data.txt')
    desi_cov_file = data_path('mock_desi_bao_rs_over_DV_cov.txt')
    info['likelihood']['bao.generic'] = {'measurements_file': desi_data_file, 'cov_file': desi_cov_file}
if args.feedback:
    info['likelihood']['hdlike.hdlike.HDLike']['data_file'] = data_path(f'hd_binnedTheorySpectra_lmin30_lmax20k_Lmin30_Lmax20k_{cmb_type}_feedback.txt')
    info['theory']['camb']['extra_args']['halofit_version'] = 'mead2020_feedback'

# use the `evaluate` sampler at a fixed cosmology:
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
expected_chi2_desi = 8.38174e-06
expected_chi2_hd_values = {'lensed': 3.60341,
        'delensed': 3.01047,
        'lensed_feedback': 3.58806,
        'delensed_feedback': 2.99985}
if args.feedback:
    expected_chi2_hd = expected_chi2_hd_values[f'{cmb_type}_feedback']
else:
    expected_chi2_hd = expected_chi2_hd_values[cmb_type]
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

if np.isclose(chi2_tot, expected_chi2_tot, atol=1e-3):
    print('Success! The chi^2 value from Cobaya matches the expected value.')
else:
    print('Test failed: the chi^2 value from Cobaya does not match the expected value.')
    print('Note: this may occur if you have modified the `example_cmbhd.yaml` file.')




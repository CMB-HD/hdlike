"""Automatically generate a Cobaya YAML file for the CMB-HD likelihood."""
import os
import argparse
import warnings
import numpy as np
import hdlike
from cobaya.yaml import yaml_load_file, yaml_dump_file

parser = argparse.ArgumentParser(description="Generate an input file for Cobaya (i.e., a `.yaml` file) for the CMB-HD likelihood.")
parser.add_argument('input_settings_filename', help='Pass a file, `input_settings_filename`, containing the input settings for the CMB-HD likelihood. A default file named `hdlike_settings.yaml` has been provided. NOTE that this is not the same as the input YAML file that is used to run Cobaya.')
parser.add_argument('cobaya_settings_filename', help='Pass a file, `cobaya_settings_filename`, containing the settings for Cobaya and CAMB. A default file named `camb_cobaya_settings.yaml` has been provided.')
args = parser.parse_args()

# default settings from `hdlike_settings.yaml`
defaults = {'delensed': True,
            'baryonic_feedback': False,
            'lmax': 20100, 'Lmax': 20100,
            'use_cmb_power_spectra': True, 
            'use_cmb_lensing_spectrum': True,
            'desi_bao': False,
            'use_fisher_proposal_widths': False,
            'output_dir': os.getcwd(), 'output_root': None,
            'yaml_dir': os.getcwd(), 'yaml_file_name': None}

# read in the settings passed by the user and set default values
hdlike_settings = yaml_load_file(args.input_settings_filename)
for name, value in hdlike_settings.items():
    if value is None:
        hdlike_settings[name] = defaults[name]
delensed = hdlike_settings['delensed']
baryonic_feedback = hdlike_settings['baryonic_feedback']
lmax = hdlike_settings['lmax']
Lmax = hdlike_settings['Lmax']
use_cmb_power_spectra = hdlike_settings['use_cmb_power_spectra']
use_cmb_lensing_spectrum = hdlike_settings['use_cmb_lensing_spectrum']
desi_bao = hdlike_settings['desi_bao']
use_fisher_proposal_widths = hdlike_settings['use_fisher_proposal_widths']
output_dir = hdlike_settings['output_dir'] 
output_root = hdlike_settings['output_root']
yaml_dir = hdlike_settings['yaml_dir'] 
yaml_file_name = hdlike_settings['yaml_file_name']

# get the path to the hdlike files and load the template Cobaya YAML file
data_dir = os.path.join(os.path.dirname(hdlike.__file__), 'data/')
data_path = lambda x: os.path.join(data_dir, x)

# settings for Cobaya
info = yaml_load_file(args.cobaya_settings_filename)

if 'classy' in info['theory'].keys():
    if delensed:
        raise ValueError("Cannot compute delensed power spectra using CLASS.")
    if baryonic_feedback:
        raise ValueError("Cannot use the HMCode2020 baryonic feedback model with CLASS.")


# update the `info` dict with the loaded settings
cmb_type = 'delensed' if delensed else 'lensed'
if baryonic_feedback:
    data_file = data_path(f'hd_binnedTheorySpectra_lmin30_lmax20k_Lmin30_Lmax20k_{cmb_type}_feedback.txt')
    hmcode_version = 'mead2020_feedback'
else:
    data_file = data_path(f'hd_binnedTheorySpectra_lmin30_lmax20k_Lmin30_Lmax20k_{cmb_type}.txt')
    hmcode_version = 'mead2016'
info['theory']['camb']['extra_args']['halofit_version'] = hmcode_version

hdlike_info = {'delensed': delensed,
               'data_file': data_file, # covmat is set within likelihood
               'baryonic_feedback': baryonic_feedback,
               'lmax': lmax, 'Lmax': Lmax,
               'use_cmb_power_spectra': use_cmb_power_spectra, 'use_cmb_lensing_spectrum': use_cmb_lensing_spectrum}
if 'likelihood' not in info.keys():
    info['likelihood'] = {'hdlike.hdlike.HDLike': hdlike_info}
else:
    info['likelihood']['hdlike.hdlike.HDLike'] = hdlike_info

if desi_bao:
    desi_data_file = data_path('mock_desi_bao_rs_over_DV_data.txt')
    desi_cov_file = data_path('mock_desi_bao_rs_over_DV_cov.txt')
    info['likelihood']['bao.generic'] = {'measurements_file': desi_data_file, 'cov_file': desi_cov_file}

# if using model with baryonic feedback and not sampling the feedback parameter, make sure it's fixed.
if baryonic_feedback:
    if 'HMCode_logT_AGN' not in info['params'].keys():
        info['params']['HMCode_logT_AGN'] = 7.8

# where to save the Cobaya output and the YAML file
root = f'cmbhd_{cmb_type}'
if (lmax < defaults['lmax']) or (Lmax < defaults['Lmax']):
    root = f'{root}_lmax{lmax}Lmax{Lmax}'
if not use_cmb_lensing_spectrum:
    root = f'{root}_tt_te_ee_bb'
elif not use_cmb_power_spectra:
    root = f'{root}_kk'
if desi_bao:
    root = f'{root}_desi'
if baryonic_feedback:
    root = f'{root}_feedback'

if output_root is None:
    output_root = root
if yaml_file_name is None:
    yaml_file_name = f'{root}.yaml'

info['output'] = os.path.join(output_dir, output_root)
print(f"The Cobaya output will be saved in the directory `{output_dir}`, with file names beginning with `{output_root}`.\n")

if use_fisher_proposal_widths:
    # get proposal widths from inverse of diagonal elements of Fisher matrix
    fisher_proposal_path = lambda x: os.path.join(data_path('proposal_cov/from_fisher/'), x)
    
    fisher_pcov_fnames = {'delensed': fisher_proposal_path('hd_fsky0pt6_lmin30_lmax20k_Lmin30_Lmax20k_withFG_delensed_lcdm_neff_mnu_theta_pcov.txt'),
                          'delensed_feedback': fisher_proposal_path('hd_fsky0pt6_lmin30_lmax20k_Lmin30_Lmax20k_withFG_delensed_lcdm_neff_mnu_feedback_pcov.txt'),
                          'delensed_bao': fisher_proposal_path('hd_fsky0pt6_lmin30_lmax20k_Lmin30_Lmax20k_withFG_delensed_desi_bao_lcdm_neff_mnu_pcov.txt'),
                          'delensed_bao_feedback': fisher_proposal_path('hd_fsky0pt6_lmin30_lmax20k_Lmin30_Lmax20k_withFG_delensed_desi_bao_lcdm_neff_mnu_feedback_pcov.txt'),
                          'lensed': fisher_proposal_path('hd_fsky0pt6_lmin30_lmax20k_Lmin30_Lmax20k_withFG_lensed_lcdm_neff_mnu_theta_pcov.txt'),
                          'lensed_feedback': fisher_proposal_path('hd_fsky0pt6_lmin30_lmax20k_Lmin30_Lmax20k_withFG_lensed_lcdm_neff_mnu_feedback_pcov.txt'),
                          'lensed_bao': fisher_proposal_path('hd_fsky0pt6_lmin30_lmax20k_Lmin30_Lmax20k_withFG_lensed_desi_bao_lcdm_neff_mnu_pcov.txt'),
                          'lensed_bao_feedback': fisher_proposal_path('hd_fsky0pt6_lmin30_lmax20k_Lmin30_Lmax20k_withFG_lensed_desi_bao_lcdm_neff_mnu_feedback_pcov.txt')}
    
    fisher_pcov_key = cmb_type
    if desi_bao:
        fisher_pcov_key = f'{fisher_pcov_key}_bao'
    if baryonic_feedback:
        fisher_pcov_key = f'{fisher_pcov_key}_feedback'
    fisher_pcov_fname = fisher_pcov_fnames[fisher_pcov_key]
    # load inverse Fisher matrix
    fisher_pcov = np.loadtxt(fisher_pcov_fname)
    # get parameter names from file header
    with open(fisher_pcov_fname, 'r') as f:
        fisher_pcov_header = f.readline()
    fisher_params = fisher_pcov_header.strip('# \n').split()
    # get the fisher matrix
    fisher_matrix = np.linalg.inv(fisher_pcov)
    # get proposal widths
    proposal_widths = {}
    for i, param in enumerate(fisher_params):
        proposal_widths[param] = 1 / np.sqrt(fisher_matrix[i,i])
    
    
    # or, use a proposal matrix from an MCMC run, if it exists
    use_proposal_matrix = True if (delensed and desi_bao) else False
    
    if 'mcmc' in info['sampler'].keys():
        if use_proposal_matrix:
            if baryonic_feedback:
                info['sampler']['mcmc']['covmat'] = data_path('proposal_cov/from_chains/hd_delensed_lcdm_nnu_mnu_bao_hmcode2020_feedback.txt')
            else:
                info['sampler']['mcmc']['covmat'] = data_path('proposal_cov/from_chains/hd_delensed_desi_bao_lcdm_nnu_mnu.txt')
            print(f"Set the proposal matrix to {info['sampler']['mcmc']['covmat']}")
        else:
            warnings.warn('No proposal matrix available; setting the `proposal` width for each parameter from the corresponding Fisher matrix.')
            for param in proposal_widths.keys():
                if type(info['params'][param]) == dict:
                    info['params'][param]['proposal'] = proposal_widths[param]

elif 'mcmc' in info['sampler'].keys():
    if 'covmat' not in info['sampler']['mcmc'].keys():
        proposal_cov_dir = os.path.join(data_path('proposal_cov'), 'from_chains')
        print(f'You have not set a proposal matrix. See `example_cmbhd.yaml` for how to include one. We have provided a few in the following directory: {proposal_cov_dir}')


# save the YAML file
yaml_fname = os.path.join(yaml_dir, yaml_file_name)
yaml_dump_file(yaml_fname, info)
print(f"\nSaved the Cobaya YAML file as {yaml_fname}.")
print(f"It is recommended that you test this by running the following command:")
print(f"\tcobaya-run {yaml_fname} --test")



"""tests the initialization of the CMB-HD likelihood and evaluate it at the fiducial cosmology"""
import os
import argparse
import numpy as np
import camb
import hdlike


# cosmological parameters for theory
cosmo_params = {'logA': 3.044,
                'ns': 0.9649,
                'theta_MC_100': 1.04092,
                'ombh2': 0.02237,
                'omch2': 0.1200,
                'mnu': 0.06,
                'nnu': 3.046,
                'tau': 0.0544,
                'HMCode_logT_AGN': 7.8, # only used with feedback
                }

parser = argparse.ArgumentParser(description='Test the initialization of the CMB-HD likelihood, and evaluate the likelihood at the fiducial cosmology.')

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

# get the path to the hdlike files
data_dir = os.path.join(os.path.dirname(hdlike.__file__), 'data/')
data_path = lambda x: os.path.join(data_dir, x)

# initialize the HD likelihood
hd_likelihood = hdlike.HDData(delensed=(not args.lensed), baryonic_feedback=args.feedback, use_desi_bao=args.desi)

# calculate the theory
lmax = 20100
cmb_spectra = ['tt', 'ee', 'bb', 'te']
hmcode_version = 'mead2020_feedback' if args.feedback else 'mead2016'
pars = camb.CAMBparams()
pars.set_cosmology(ombh2=cosmo_params['ombh2'], 
                   omch2=cosmo_params['omch2'], 
                   cosmomc_theta=cosmo_params['theta_MC_100']/100, 
                   tau=cosmo_params['tau'], 
                   num_massive_neutrinos=1, 
                   mnu=cosmo_params['mnu'], 
                   nnu=cosmo_params['nnu'])
pars.InitPower.set_params(As=np.exp(cosmo_params['logA']) * 1e-10, 
                          ns=cosmo_params['ns'])
pars.set_for_lmax(int(lmax)+500, lens_potential_accuracy=30, lens_margin=2050)
pars.set_accuracy(AccuracyBoost=1.1, lSampleBoost=3.0, lAccuracyBoost=3.0, DoLateRadTruncation=False)
pars.NonLinear = camb.model.NonLinear_both
pars.NonLinearModel.set_params(hmcode_version, 
                               HMCode_logT_AGN=cosmo_params['HMCode_logT_AGN'])
print('Calculating theory ...')
theo = {}
results = camb.get_results(pars)
results.calc_power_spectra()
powers = results.get_cmb_power_spectra(pars, CMB_unit='muK', lmax=lmax, raw_cl=True, spectra=['total'])
for i, s in enumerate(cmb_spectra):
    theo[s] = powers['total'][:lmax+1,i].copy()
    theo[s][:2] = 0
theo['pp'] = results.get_lens_potential_cls(lmax=lmax, raw_cl=True)[:,0] 
theo['pp'][:2] = 0
loglike_hd = hd_likelihood.log_likelihood(theo, camb_results=results)
chi2_hd = -2 * loglike_hd
if args.desi:
    z = hd_likelihood.get_desi_redshifts()
    rs_dv = results.get_BAO(z, pars)[:,0]
    loglike_desi = hd_likelihood.log_likelihood_desi(rs_dv)
    chi2_desi = -2 * loglike_desi
else:
    chi2_desi = 0
chi2_tot = chi2_hd + chi2_desi

# expected results
expected_chi2_desi = 2.60556e-9
expected_chi2_hd_values = {'lensed': 3.91772e-10,
        'delensed': 3.05169e-10,
        'lensed_feedback': 1.5083e-22,
        'delensed_feedback': 1.20353e-19}
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

print(f'\nchi^2 for CMB-HD = {chi2_hd}; expected chi^2 = {expected_chi2_hd} (difference = {chi2_hd_diff})')
if args.desi:
    print(f'chi^2 for DESI = {chi2_desi}; expected chi^2 = {expected_chi2_desi} (difference = {chi2_desi_diff})')
    print(f'Total chi^2 for CMB-HD + DESI = {chi2_tot}; expected chi^2 = {expected_chi2_tot} (difference = {chi2_tot_diff})')

if np.isclose(chi2_tot, expected_chi2_tot, atol=1e-3):
    print('Success! The chi^2 value matches the expected value.')
else:
    print('Test failed: the chi^2 value does not match the expected value.')
    print('Note: this may occur if you have modified the theory calculation within this file in any way.')




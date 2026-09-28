"""Test the initialization of the CMB-HD likelihood, and evaluate it at
the fiducial cosmology.

"""

import os
import argparse
import numpy as np
from hd_mock_data import hd_data
import hdlike


# command-line args:
formatter_class = argparse.ArgumentDefaultsHelpFormatter
description = ("Test the initialization of the CMB-HD likelihood, and "
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

# print out info about the input arguments:
cmb_type = 'lensed' if args.lensed else 'delensed'
print(f'Using {cmb_type} CMB-HD data.')
if args.desi:
    print('Using DESI BAO data.')
if args.feedback:
    print('Using a model with baryonic feedback.')
if args.use_class:
    print("Using CLASS.")

# initialize the HD likelihood:
hd_likelihood = hdlike.HDData(delensed=(not args.lensed),
                              baryonic_feedback=args.feedback,
                              use_desi_bao=args.desi,
                              use_class=args.use_class)

# calculate the theory:
print('Calculating theory ...')
hd_datalib = hd_data.HDMockData()
lmax = hd_datalib.lmax
cmb_spectra = ['tt', 'ee', 'bb', 'te']
if args.desi:
    z = hd_likelihood.get_desi_redshifts()
theo = {}
if args.use_class:
    results = None # passed to likelihood
    from classy import Class
    params = hd_datalib.class_settings(baryonic_feedback=args.feedback)
    class_results = Class()
    class_results.empty()
    class_results.set(params)
    class_results.compute()
    # lensed CMB power spectra:
    class_cls = class_results.lensed_cl(lmax)
    TCMB = class_params.get('T_cmb', 2.7255) * 1e6
    for s in cmb_spectra:
        theo[s] = class_cls[s] * TCMB**2
    # lensing potential power spectrum:
    theo['pp'] = class_results.raw_cl(lmax)['pp']
    theo['pp'][:2] = 0
    if args.desi:
        rs = class_results.rs_drag() # Mpc
        d_a = class_results.angular_distance(z) # Mpc (physical D_A)
        hubble = class_results.Hubble(z) # 1/Mpc, i.e. H(z)/c
        d_m = (1.0 + z) * d_a # Mpc (comoving angular-diameter distance)
        d_v = (d_m**2 * z / hubble)**(1.0 / 3.0) # Mpc
        rs_dv = rs / d_v
else:
    import camb
    params = hd_datalib.camb_settings(baryonic_feedback=args.feedback)
    pars = camb.set_params(**params)
    results = camb.get_results(pars)
    # CMB power spectra:
    results.calc_power_spectra()
    powers = results.get_cmb_power_spectra(pars, CMB_unit='muK', lmax=lmax,
                                           raw_cl=True, spectra=['total'])
    for i, s in enumerate(cmb_spectra):
        theo[s] = powers['total'][:lmax+1,i].copy()
        theo[s][:2] = 0
    # lensing potential power spectrum:
    theo['pp'] = results.get_lens_potential_cls(lmax=lmax, raw_cl=True)[:,0]
    theo['pp'][:2] = 0
    if args.desi:
        rs_dv = results.get_BAO(z, pars)[:,0]

# calculate the log-likelihood:
loglike_hd = hd_likelihood.log_likelihood(theo, camb_results=results)
chi2_hd = -2 * loglike_hd
if args.desi:
    loglike_desi = hd_likelihood.log_likelihood_desi(rs_dv)
    chi2_desi = -2 * loglike_desi
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
    # TODO: update these
    expected_chi2_desi = 4.380975e-3
    expected_chi2_hd_values = {'lensed': 2.627772,
                               'delensed': 3.016197,
                               'lensed_feedback': 2.558651,
                               'delensed_feedback': 2.976135}
chi2key = f'{cmb_type}_feedback' if args.feedback else cmb_type
expected_chi2_hd = expected_chi2_hd_values[chi2key]
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

if np.isclose(chi2_tot, expected_chi2_tot, atol=0.05):
    print('Success! The chi^2 value matches the expected value.')
else:
    print('Test failed: the chi^2 value does not match the expected value.')
    print('Note: this may occur if you have modified the theory calculation within this file in any way.')


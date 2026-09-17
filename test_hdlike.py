"""tests the initialization of the CMB-HD likelihood and evaluate it at the fiducial cosmology"""
import os
import argparse
import numpy as np
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

# the same cosmology in terms of the CLASS parameters, used with `--class`.
# NOTE that H0 is used in place of theta_MC, and that there are three massive
#  neutrinos with N_eff = 3.044, to match the CMB-HD bandpowers that were
#  calculated with CLASS (see the README).
class_cosmo_params = {'ln_A_s_1e10': 3.044,
                      'n_s': 0.9649,
                      'H0': 67.36,
                      'omega_b': 0.02237,
                      'omega_cdm': 0.1200,
                      'mnu': 0.06,
                      'N_eff': 3.044,
                      'tau_reio': 0.0544,
                      }

parser = argparse.ArgumentParser(description='Test the initialization of the CMB-HD likelihood, and evaluate the likelihood at the fiducial cosmology.')

parser.add_argument('--lensed',  action='store_true', help='Use lensed CMB data (delensed data is used by default).')
parser.add_argument('--desi', action='store_true', help='Include DESI BAO data (it is not included by default).')
parser.add_argument('--feedback', action='store_true', help='Use a model with baryonic feedback effects (a CDM-only model is used by default).')
parser.add_argument('--class', dest='use_class', action='store_true', help='Calculate the theory with CLASS instead of CAMB, and compare it to the CMB-HD bandpowers that were calculated with CLASS. Requires `--lensed`, and cannot be used with `--feedback`.')


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

# initialize the HD likelihood
# (the CMB-HD bandpowers calculated with CLASS are only available for v1.2 of the data)
hd_data_version = 'v1.2' if args.use_class else 'latest'
hd_likelihood = hdlike.HDData(delensed=(not args.lensed), baryonic_feedback=args.feedback, use_desi_bao=args.desi, hd_data_version=hd_data_version, use_class=args.use_class)

# calculate the theory
lmax = 20100
cmb_spectra = ['tt', 'ee', 'bb', 'te']
theo = {}
if args.use_class:
    import yaml
    from classy import Class
    # use the same CLASS settings as `class_cobaya_settings.yaml`
    settings_fname = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'class_cobaya_settings.yaml')
    with open(settings_fname, 'r') as f:
        class_settings = yaml.safe_load(f)['theory']['classy']['extra_args']
    # CLASS reads the BBN table itself, so it needs the absolute path
    data_dir = os.path.join(os.path.dirname(hdlike.__file__), 'data/')
    class_settings['sBBN file'] = os.path.join(data_dir, 'PRIMAT21_class_format.dat')
    # split the sum of the neutrino masses evenly over the three massive
    #  neutrinos to get CLASS's `m_ncdm`, and convert N_eff to CLASS's `N_ur`
    #  (see `class_cobaya_settings.yaml`)
    n_ncdm = class_settings['N_ncdm']
    m_ncdm = class_cosmo_params['mnu'] / n_ncdm
    class_params = {'ln_A_s_1e10': class_cosmo_params['ln_A_s_1e10'],
                    'n_s': class_cosmo_params['n_s'],
                    'H0': class_cosmo_params['H0'],
                    'omega_b': class_cosmo_params['omega_b'],
                    'omega_cdm': class_cosmo_params['omega_cdm'],
                    'tau_reio': class_cosmo_params['tau_reio'],
                    'm_ncdm': ','.join([str(m_ncdm)] * n_ncdm),
                    'N_ur': class_cosmo_params['N_eff'] - n_ncdm * 1.0131966}
    print('Calculating theory ...')
    results = Class()
    results.set({**class_settings, **class_params})
    results.compute()
    # CLASS returns dimensionless spectra, so convert them to uK^2
    T_cmb_muK = class_settings['T_cmb'] * 1e6
    lensed_cls = results.lensed_cl(lmax)
    for s in cmb_spectra:
        theo[s] = lensed_cls[s][:lmax+1].copy() * T_cmb_muK**2
        theo[s][:2] = 0
    theo['pp'] = results.raw_cl(lmax)['pp'][:lmax+1].copy()
    theo['pp'][:2] = 0
    loglike_hd = hd_likelihood.log_likelihood(theo)
else:
    import camb
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
    if args.use_class:
        # calculate r_s / d_V(z) in the same way as CAMB's `get_BAO`
        rs = results.rs_drag()
        rs_dv = np.zeros(len(z))
        for i, zi in enumerate(z):
            d_m = (1 + zi) * results.angular_distance(zi)
            d_v = (d_m**2 * zi / results.Hubble(zi))**(1/3)
            rs_dv[i] = rs / d_v
    else:
        rs_dv = results.get_BAO(z, pars)[:,0]
    loglike_desi = hd_likelihood.log_likelihood_desi(rs_dv)
    chi2_desi = -2 * loglike_desi
else:
    chi2_desi = 0
chi2_tot = chi2_hd + chi2_desi

# expected results
expected_chi2_desi_values = {'camb': 2.60556e-9,
        'class': 7.36069e-04}
expected_chi2_hd_values = {'lensed': 3.91772e-10,
        'delensed': 3.05169e-10,
        'lensed_feedback': 1.5083e-22,
        'delensed_feedback': 1.20353e-19,
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




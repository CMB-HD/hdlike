"""Save a YAML file to run Cobaya using the CMB-HD defaults."""

import os
import argparse
import warnings
from cobaya.yaml import yaml_load_file, yaml_dump_file
from hd_mock_data import hd_data
from hdlike import hdlike


# ---------- command line args ----------

# text printed out when `--help` is passed:
formatter_class = argparse.ArgumentDefaultsHelpFormatter
description = ("Generate an input YAML file to use the CMB-HD likelihood with "
               "Cobaya. By default, a YAML file will be generated to use the "
               "CMB-HD likelihood, including the latest available CMB-HD mock "
               "data with delensed CMB power spectra the CMB lensing "
               "convergence power spectrum, sample the six LCDM parameters with "
               "the `mcmc` sampler, and use CAMB with the CMB-HD accuracy "
               "parameters for the theory calculation. You may change these "
               "defaults by passing the arguments described below.")
file_help = ("Path to an existing Cobaya YAML file (optional). If provided, "
             "a copy of this file will be saved with the CMB-HD likelihood "
             "(and, if `--bao` is passed, the BAO likelihood with the mock "
             "DESI data) added to (or updated in) the `likelihood` block, "
             "the `theory` block updated with the CMB-HD accuracy settings, "
             "and any other missing blocks added with the CMB-HD defaults. "
             "Note: if your YAML file has a `theory` block, we will only "
             "check if it contains a key for `camb` or `classy`; if neither "
             "are found, one of them will be added (`camb` by default). ")
output_help = ("The output prefix for the Cobaya run. If an existing YAML "
               "file is passed with the `output` defined, it will be replaced.")
fname_help = ("The name of the YAML file. By default, if `--output` was passed "
              "(or a YAML file was passed with the `output` defined), and it "
              "is not the path to a directory, e.g. if `--output=chains/root`, "
              "the YAML file will be `root.yaml`. Otherwise, a default file "
              "name, based on the information passed in, will be used..")
overwrite_help = ("Whether to overwrite an existing YAML file.")
lensed_help = ("Whether to use lensed CMB-HD mock data instead of delensed; "
               "the latter is used by default. If you passed an existing "
               "YAML file with `delensed=True`, passing `--lensed` will "
               "overwrite this.")
bao_help = ("Whether to include mock DESI BAO data provided with `hdlike` "
            "by adding the `bao.generic` likelihood to the `likelihood` block. "
            "An error will be raised if the `likelihood` block of an existing "
            "YAML file already exists with different data files.")
class_help = ("Whether to use CLASS instead of CAMB for the theory calculation. "
              "Ignored if you passed an existing YAML file with a `theory` "
              "block that already has a key for `camb` or `classy`. Note that: "
              "(1) CLASS does not calculate delensed power spectra, so you must "
              "pass `--lensed` if you are using CLASS; (2) The CMB-HD "
              "bandpowers and CLASS accuracy settings are only available for "
              "CMB-HD mock data versions 1.2 and higher; (3) CLASS must be "
              "modified in order to use the default CMB-HD accuracy settings, "
              "and to vary the effective number of relativistic species with "
              "three massive neutrinos; see the `hdlike` readme file and "
              "Appendix A of Cheslog et. al. (2026).")
feedback_help = ("Use the HMCode 2020 non-linear model with baryonic feedback "
                 "effects, instead of the default 2016 HMCode CDM-only model. "
                 "If `--feedback` is passed and `--feedback_fixed` is not passed, "
                 "the parameter characterizing the feedback strength will be "
                 "added to the `params` block as a sampled parameter. Note that "
                 "passing `--feedback` will overwrite the non-linear model in "
                 "the `theory` block of an existing YAML file.")
vary_feedback_help = ("Whether to fix the baryonic feedback parameter of the "
                      "HMCode 2020 non-linear model with baryonic feedback "
                      "effects, if that model is being used. By default, it "
                      "will be sampled.")
theo_help = ("If you passed an existing YAML file that already contains a "
             "`theory` block with a key for `camb` or `classy`, pass "
             "`--overwrite_theory` to overwrite any existing parameters "
             "defined in the `camb` or `classy` block with the fiducial "
             "values for CMB-HD. By default, any existing theory parameter "
             "values will not be changed, but any un-defined CMB-HD theory "
             "parameters will be added.")
H0_help = ("Whether to sample the Hubble parameter H0 instead of the "
           "angular size of the sound horizon at last scattering. "
           "By default, H0 is a derived parameter.")
Neff_help = ("Whether to vary the effective number of relativistic species. "
             "Ignored if this parameter already appears in the `params` block.")
mnu_help = ("Whether to vary the sum of the neutrino masses. Note that three "
            "massive neutrinos are assumed for CMB-HD mock data versions "
            ">= 1.2, and when using CLASS, the total mass is split evenly "
            "between them. Ignored if this parameter already appears in the "
            "`params` block.")
nrun_help = ("Whether to vary the running of the scalar spectral index. "
             "Ignored if this parameter already appears in the `params` block.")
cov_help = ("Whether to use a parameter covariance matrix provided by "
            "`hdlike`, from a previous run on the mock CMB-HD and DESI "
            "BAO data, as the proposal matrix (`covmat` in the `mcmc` "
            "sampler block). Note that we do not provide this matrix for "
            "every possible set of sampled parameters and data, but the "
            "provided files should be a good starting point. If you passed "
            "a YAML file that already has a `covmat` under `mcmc` in the "
            "`sampler` block, its value will be overwritten. Ignored if a "
            "different sampler is being used.")
hd_help = ("The version of the CMB-HD mock data to use. By default, the "
           "latest available data will be used. Overrides any existing "
           "`hd_data_version` specified in the input YAML file, if one "
           "was passed.")

# define and parse the arguments:
parser = argparse.ArgumentParser(formatter_class=formatter_class, description=description)
parser.add_argument('input_file', nargs='?', help=file_help)
parser.add_argument('-o', '--output', help=output_help)
parser.add_argument('-f', '--filename', help=fname_help)
parser.add_argument('--overwrite', action='store_true', help=overwrite_help)
parser.add_argument('--lensed', action='store_true', help=lensed_help)
parser.add_argument('--bao', action='store_true', help=bao_help)
parser.add_argument('--use_class', action='store_true', help=class_help)
parser.add_argument('--feedback', action='store_true', help=feedback_help)
parser.add_argument('--feedback_fixed', action='store_true', help=vary_feedback_help)
parser.add_argument('--overwrite_theory', action='store_true', help=theo_help)
parser.add_argument('--sampleH0', action='store_true', help=H0_help)
parser.add_argument('--Neff', action='store_true', help=Neff_help)
parser.add_argument('--mnu', action='store_true', help=Neff_help)
parser.add_argument('--nrun', action='store_true', help=nrun_help)
parser.add_argument('--use_param_covmat', action='store_true', help=cov_help)
parser.add_argument('--hd_data_version', default='latest', help=hd_help,
                    choices=[*hd_data.HDMockData.data_versions, 'latest'])
args = parser.parse_args()

# say this after warnings:
overwrite_warn_msg = ("Note: If you resolve this warning by re-running this "
                      "script again with the same output YAML file name, pass "
                      "`--overwrite` to overwrite the file; otherwise, pass a "
                      "new file name.")


# ---------- set up ----------

hd_datalib = hd_data.HDMockData(version=args.hd_data_version)

info = {} if (args.input_file is None) else yaml_load_file(args.input_file)

output = args.output if (args.output is not None) else info.get('output')
if output is None:
    msg = ("No `output` provided for the YAML file, so Cobaya will not save any"
           f" output files when run with this YAML file. {overwrite_warn_msg}")
    warnings.warn(msg)
else:
    info['output'] = output
    if 'resume' not in info:
        info['resume'] = True

# check if CLASS is being used:
use_class = args.use_class
if 'theory' in info:
    if 'classy' in info['theory']:
        use_class = True
# check if CLASS is an option for this version of HD mock data:
hd_datalib.get_compatible_version(hd_datalib.class_theo_versions, 
                                  'mock CMB-HD data for CLASS')

    
# ---------- likelihood block ----------

def hdlike_path(relative_path):
    hdlike_dir = os.path.split(hdlike.__file__)[0]
    return os.path.join(hdlike_dir, relative_path)

# load in the defaults for HDLike:
hdlike_defaults = yaml_load_file(hdlike_path('HDLike.yaml'))

# check if 'info' has a `likelihood` block, and, if so, 
# whether it already has a key for hdlike:
if 'likelihood' not in info:
    info['likelihood'] = {}
hdlike_info = info['likelihood'].get('hdlike.hdlike.HDLike', {})
hdlike_info = {**hdlike_info, 
               'delensed': (not args.lensed), 
               'baryonic_feedback': args.feedback}
if hd_datalib.version != hd_datalib.latest_version:
    hdlike_info['hd_data_version'] = hd_datalib.version

# make sure each key in `hdlike_info` has a corresponding key in `hdlike_defaults`:
unrecognized_keys = [key for key in hdlike_info if (key not in hdlike_defaults)]
if len(unrecognized_keys) > 0:
    raise ValueError("The following keyword arguments are not "
                     f"accepted by `HDLike`: {unrecognized_keys}")
# check if `delensed=True` and `use_class=True` were both passed:
if use_class:
    delensed = hdlike_info.get('delensed', not args.lensed)
    if delensed:
        raise ValueError("When using CLASS, you must either pass `--lensed`, "
                         "or set `delensed=False` in the `HDLike` likelihood "
                         "block of an existing YAML file. Note that the "
                         "default is `delensed=True`.")

# only save the non-default settings in the YAML file:
hdlike_settings = {}
for key, value in hdlike_info.items():
    # we are only comparing `str`, `int`, or `bool`, so this should be ok:
    if value != hdlike_defaults[key]:
        hdlike_settings[key] = value
# will want to know if `delensed` is `True` or `False` later:
delensed = hdlike_settings.get('delensed', hdlike_defaults['delensed'])

if len(hdlike_settings) > 0:
    info['likelihood']['hdlike.hdlike.HDLike'] = hdlike_settings
else:
    info['likelihood']['hdlike.hdlike.HDLike'] = None

# BAO:
if args.bao: # add paths to mock DESI files:
    desi_data_fname = 'mock_desi_bao_rs_over_DV_data.txt'
    desi_cov_fname = 'mock_desi_bao_rs_over_DV_cov.txt'
    # check if the `bao.generic` likelihood is already being used with different files:
    if 'bao.generic' in info['likelihood']:
        bao_block = info['likelihood']['bao.generic']
        bao_data_file = bao_block.get('measurements_file', desi_data_fname)
        bao_cov_file = bao_block.get('cov_file', desi_data_fname)
        if (desi_data_fname not in bao_data_file) or (desi_cov_fname not in bao_cov_file):
            raise ValueError("You passed `--bao`, but your input YAML file "
                             "already has a `bao.generic` block in the "
                             "`likelihood` block.")
    if use_class:
        desi_data_fname = f'class_{desi_data_fname}'
    bao_settings = {'measurements_file': hdlike_path(os.path.join('data', desi_data_fname)),
                    'cov_file': hdlike_path(os.path.join('data', desi_cov_fname))}
    info['likelihood']['bao.generic'] = bao_settings
    

# ---------- params block ----------
        
tau_prior = 0.007 if (hd_datalib.version in ['v1.0', 'v1.1']) else 0.005

H0_fixed = {'min': 20, 'max': 100, 'latex': r'H_0'}
H0_varied = {'prior': {'min': 20, 'max': 100}, 
             'ref': {'dist': 'norm', 'loc': 67.36, 'scale': 0.3}, 
             'proposal': 0.015, 'latex': r'H_0'}
H0_settings = H0_varied if args.sampleH0 else H0_fixed

# first, define settings for all available params (available here, not in general);
# then translate them to CLASS names; then only keep the sampled params in the `params` block
hdlike_camb_params_settings = {'ombh2': {'prior': {'min': 0.005, 'max': 0.1},
                                         'ref': {'dist': 'norm', 'loc': 0.0224, 'scale': 5e-05},
                                         'proposal': 9.7e-06,
                                         'latex': r'\Omega_\mathrm{b} h^2'},
                               'omch2': {'prior': {'min': 0.001, 'max': 0.99},
                                         'ref': {'dist': 'norm', 'loc': 0.12, 'scale': 0.0005},
                                         'proposal': 2.1e-05,
                                         'latex': r'\Omega_\mathrm{c} h^2'},
                               'logA': {'prior': {'min': 2, 'max': 4},
                                        'ref': {'dist': 'norm', 'loc': 3.05, 'scale': 0.001},
                                        'proposal': 0.0003,
                                        'latex': r'\log(10^{10} A_\mathrm{s})'},
                               'As': {'value': 'lambda logA: 1e-10*np.exp(logA)',
                                      'latex': r'A_\mathrm{s}',
                                      'derived': True},
                               'ns': {'prior': {'min': 0.8, 'max': 1.2},
                                      'ref': {'dist': 'norm', 'loc': 0.965, 'scale': 0.002},
                                      'proposal': 0.0002,
                                      'latex': r'n_\mathrm{s}'},
                               'tau': {'prior': {'dist': 'norm', 'loc': 0.054, 'scale': tau_prior},
                                       'ref': {'dist': 'norm', 'loc': 0.054, 'scale': 0.006},
                                       'proposal': 0.0002,
                                       'latex': r'\tau'},
                               'theta_MC_100': {'prior': {'min': 0.5, 'max': 10},
                                                'ref': {'dist': 'norm', 'loc': 1.04071, 'scale': 0.0004},
                                                'proposal': 5e-05,
                                                'latex': r'100\theta_\mathrm{MC}',
                                                'drop': True,
                                                'renames': 'theta'},
                               'cosmomc_theta': {'value': 'lambda theta_MC_100: 1.e-2*theta_MC_100',
                                                 'derived': False},
                               'H0': H0_settings,
                               'nnu': {'prior': {'min': 0.05, 'max': 10},
                                       'ref': {'dist': 'norm', 'loc': 3.044, 'scale': 0.005},
                                       'proposal': 0.0013,
                                       'latex': r'N_\mathrm{eff}'},
                               'mnu': {'prior': {'min': 0, 'max': 5},
                                       'ref': {'dist': 'norm', 'loc': 0.06, 'scale': 0.01},
                                       'proposal': 0.0014,
                                       'latex': r'\sum m_\nu'},
                               'nrun': {'prior': {'min': -0.2, 'max': 0.2},
                                        'ref': {'dist': 'norm', 'loc': 0.0, 'scale': 0.005},
                                        'proposal': 0.001,
                                        'latex': r'\alpha_s'},
                               'HMCode_logT_AGN': {'prior': {'dist': 'norm', 'loc': 7.8, 'scale': 0.0006 * 7.8},
                                                   # !! TODO !!
                                                   #'prior': {'min': 7.6, 'max': 8.0}, 
                                                   'ref': {'dist': 'norm', 'loc': 7.8, 'scale': 0.03}, 
                                                   'proposal': 0.006, 
                                                   'latex': r'\log_{10}(T_\mathrm{AGN}/\mathrm{K})'},
                               'sigma8': {'latex': r'\sigma_8'},
                              }
# for CLASS, first add parameters that can't just
# be copied from the CAMB params dict and renamed:
hdlike_class_params_settings = {'A_s': {'value': 'lambda ln_A_s_1e10: 1e-10*np.exp(ln_A_s_1e10)',
                                        'latex': r'A_\mathrm{s}',
                                        'derived': True},
                               'theta_s_100': {'prior': hdlike_camb_params_settings['theta_MC_100']['prior'],
                                                'ref': {**hdlike_camb_params_settings['theta_MC_100']['ref'], 
                                                        'loc': 1.04165},
                                                'proposal': hdlike_camb_params_settings['theta_MC_100']['proposal'],
                                                'latex': r'100\theta_\mathrm{s}'},
                                'mnu': {**hdlike_camb_params_settings['mnu'], 'drop': True},
                                'm_ncdm': {'value': 'lambda mnu: str(mnu/3.) + "," + str(mnu/3.) + "," + str(mnu/3.)',
                                           'latex': 'm_\\nu', 'derived': False},
                                'N_eff': {**hdlike_camb_params_settings['nnu'], 'drop': True},
                                'N_ur': {'value': 'lambda N_eff: N_eff - 3*1.0131966',
                                         'latex': 'N_\\mathrm{ur}', 'derived': False},
                               }
# for the remaining parameters, the only difference between
# CAMB and CLASS in the params block is the parameter name:
camb2class = {'logA': 'ln_A_s_1e10', 'As': 'A_s', 'ns': 'n_s', 'tau': 'tau_reio',
              'ombh2': 'omega_b', 'omch2': 'omega_cdm', 'H0': 'H0',
              'nrun': 'alpha_s',  'HMCode_logT_AGN': 'log10T_heat_hmcode', 
              'sigma8': 'sigma8'}
for camb_name, class_name in camb2class.items():
    if class_name not in hdlike_class_params_settings:
        hdlike_class_params_settings[class_name] = hdlike_camb_params_settings[camb_name]

# make the params block dict:
params = ['logA', 'As', 'ns', 'tau', 'ombh2', 'omch2', 'H0']
if use_class:
    params = [camb2class[p] for p in params]
if not args.sampleH0:
    theta_names = ['theta_s_100'] if use_class else ['theta_MC_100', 'cosmomc_theta']
    params = [*params, *theta_names]
if args.Neff:
    nnu_names = ['N_eff', 'N_ur'] if use_class else ['nnu']
    params = [*params, *nnu_names]
if args.mnu:
    mnu_names = ['mnu', 'm_ncdm'] if use_class else ['mnu']
    params = [*params, *mnu_names]
if args.nrun:
    nrun = 'alpha_s' if use_class else 'nrun'
    params.append(nrun)
if args.feedback and (not args.feedback_fixed):
    feedback_param = 'log10T_heat_hmcode' if use_class else 'HMCode_logT_AGN'
    params.append(feedback_param)
params.append('sigma8')
if use_class:
    params_block = {p: hdlike_class_params_settings[p] for p in params}
else:
    params_block = {p: hdlike_camb_params_settings[p] for p in params}

if 'params' not in info:
    info['params'] = params_block
else: # don't overwrite, just add additional params
    info['params'] = {**params_block, **info['params']} 


# ---------- theory block ----------
        
theo_name = 'classy' if use_class else 'camb'
if 'theory' not in info:
    info['theory'] = {}
if theo_name not in info['theory']:
    info['theory'][theo_name] = {}
if 'extra_args' not in info['theory'][theo_name]:
    info['theory'][theo_name]['extra_args'] = {}
theo_settings = info['theory'][theo_name]['extra_args']

# load in all of the HD settings for CAMB or CLASS:
if use_class:
    all_theo_settings = hd_datalib.class_settings(baryonic_feedback=args.feedback)
else:
    # check if user has already defined the CAMB parameter named 
    # `lens_margin` or `lens_output_margin`, depending on the CAMB version:
    v2_camb_names = None
    if 'lens_margin' in theo_settings:
        v2_camb_names = False
    elif 'lens_output_margin' in theo_settings:
        v2_camb_names = True
    all_theo_settings = hd_datalib.camb_settings(baryonic_feedback=args.feedback,
                                                 v2_camb_names=v2_camb_names)
    all_theo_settings['theta_H0_range'] = [20, 100]
    # if it wasn't already defined, warn user about this CAMB parameter name:
    if v2_camb_names is None:
        use_v2_camb_names = hd_datalib._use_v2camb
        # `use_v2_camb_names` is `True` or `False` based on `camb.__version__`
        # if CAMB can be imported, or it is `None` if `import camb` fails,
        # in which case we default to the name used in the latest version:
        name = 'lens_margin' if (use_v2_camb_names is False) else 'lens_output_margin'
        if use_v2_camb_names is not None:
            import camb # get the exact version
            camb_version = camb.__version__
            camb_info = (f"You are using CAMB version {camb_version} "
                         "(determined by importing `camb`), so the YAML "
                         f"file will contain `{name}`.")
        else:
            camb_info ("CAMB could not be imported with `import camb`, so the "
                       "YAML file will contain `lens_output_margin`.")
        msg = ("The CAMB accuracy parameters used for CMB-HD includes the "
               "parameter named `lens_output_margin` in CAMB versions >= 2.0.0, "
               "which is named `output_margin` in previous CAMB versions. "
               f"{camb_info} If you are use a CAMB version that is incompatible "
               "with this name when you run Cobaya, you will need to edit the "
               "YAML file to update the name of this parameter.")
        warnings.warn(msg)
    
# update the theory block:
if args.overwrite_theory:
    theo_settings = {**theo_settings, **all_theo_settings}
else:
    theo_settings = {**all_theo_settings, **theo_settings}
# make sure the correct non-linear model is used if `--feedback` was passed:
if args.feedback:
    if use_class:
        theo_settings['hmcode_version'] = '2020_baryonic_feedback'
    else:
        theo_settings['halofit_version'] = 'mead2020_feedback'
# remove any parameters that appear in the `params` block:
theo_settings = {k: v for (k, v) in theo_settings.items() if (k not in info['params'])}
info['theory'][theo_name]['extra_args'] = theo_settings


# ---------- sampler block ----------

mcmc_settings = {'Rminus1_stop': 0.01, 'Rminus1_cl_stop': 0.2}
if use_class:
    mcmc_settings['max_tries'] = '400d'

if 'sampler' not in info:
    info['sampler'] = {'mcmc': mcmc_settings}
elif 'mcmc' in info['sampler']:
    info['sampler']['mcmc'] = {**mcmc_settings, **info['sampler']['mcmc']}

if args.use_param_covmat and ('mcmc' in info['sampler']):
    rel_pcov_path = os.path.sep.join(['data', 'proposal_cov', 'from_chains'])
    pcov_path = lambda x: os.path.join(hdlike_path(rel_pcov_path), x)
    
    if delensed:
        if args.feedback:
            pcov_fname = pcov_path('hd_delensed_lcdm_nnu_mnu_bao_hmcode2020_feedback.txt')
        else:
            pcov_fname = pcov_path('hd_delensed_desi_bao_lcdm_nnu_mnu.txt')
    elif use_class:
        if args.Neff or args.mnu or args.nrun:
            pcov_fname = pcov_path('hd_lensed_desi_bao_class_lcdm_nrun_nnu_mnu.txt')
        else:
            pcov_fname = pcov_path('hd_lensed_desi_bao_class_lcdm.txt')
    else:
        if args.Neff or args.mnu or args.nrun:
            pcov_fname = pcov_path('hd_lensed_desi_bao_camb_lcdm_nrun_nnu_mnu.txt')
        else:
            pcov_fname = pcov_path('hd_lensed_desi_bao_camb_lcdm.txt')
    info['sampler']['mcmc']['covmat'] = pcov_fname

    
# ---------- save the file ----------

# default YAML file name if none is provided:
cmb_type = 'delensed' if delensed else 'lensed'
file_root_info = ['hdlike']
if hd_datalib.version != hd_datalib.latest_version:
    file_root_info.append(f'{hd_datalib.version}HDdata')
has_cmb = hdlike_settings.get('has_cmb_power_spectra', hdlike_defaults['has_cmb_power_spectra'])
has_clkk = hdlike_settings.get('has_cmb_lensing_spectrum', hdlike_defaults['has_cmb_lensing_spectrum'])
if not has_clkk:
    file_root_info.append(f'{cmb_type}_tt_te_ee_bb')
elif not has_cmb:
    file_root_info.append('kk')
else:
    file_root_info.append(cmb_type)
if args.bao:
    file_root_info.append('mockDESI')
if args.feedback:
    file_root_info.append('feedback')
lcdm_info = 'lcdmH0' if args.sampleH0 else 'lcdm' 
file_root_info.append(lcdm_info)
if args.Neff:
    file_root_info.append('Neff')
if args.mnu:
    file_root_info.append('mnu')
if args.nrun:
    file_root_info.append('nrun')
if use_class:
    file_root_info.append('CLASS')
file_root = '_'.join(file_root_info)

# get the file name from info passed in, or use the default,
# and check if the file already exists:
if args.filename is not None:
    yaml_file = args.filename
elif (output is not None) and (not output.endswith(os.path.sep)):
    root = output.split(os.path.sep)[-1]
    yaml_file = f'{root}.yaml'
else:
    yaml_file = f'{file_root}.yaml'
if os.path.exists(yaml_file) and (not args.overwrite):
    errmsg = (f"The YAML file {yaml_file} already exists. You must either "
              "pass a new file name, or pass `--overwrite` to overwrite "
              "an existing file.")
    raise FileExistsError(errmsg)


yaml_dump_file(yaml_file, info, error_if_exists=(not args.overwrite))
print(f"\nSaved the Cobaya YAML file to {yaml_file}.")
print(f"It is recommended that you test this by running the following command:")
print(f"    cobaya-run {yaml_file} --test")


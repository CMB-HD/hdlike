# Readme for HDLike

This is a mock CMB-HD likelihood including lensed and delensed $TT/TE/EE/BB$ CMB + lensing $\kappa\kappa$ spectra from multipoles 30 to 20,000.  We also include a mock DESI BAO likelihood. The likelihood can be used with Cobaya.  Please cite [MacInnis, Sehgal, and Rothermel (2023)](https://arxiv.org/abs/2309.03021). The option to calculate the theory with CLASS instead of CAMB was added for [Cheslog et. al. (2026)](https://arxiv.org/abs/XXXX.XXXXX); please also cite that work if you use it.


## Installation of likelihood


### Installation requirements

To use the CMB-HD likelihood, you must install Python version >= 3, [NumPy](https://numpy.org/), and [hdMockData](https://github.com/CMB-HD/hdMockData).

To test the likelihood by running `test_hdlike.py`, you must also install [CAMB](https://camb.readthedocs.io/en/latest/) (or CLASS, to run the test with the `--class` flag).

To use the likelihood with Cobaya, you must have Python version >= 3.8 and [Cobaya](https://cobaya.readthedocs.io/en/latest/index.html).

To calculate the theory with [CLASS](https://github.com/lesgourg/class_public) instead of CAMB, you must also install CLASS and its Python wrapper `classy`, and modify CLASS as described in the [Using CLASS instead of CAMB](#using-class-instead-of-camb) section below.

(The likelihood has been tested with Cobaya version 3.3.2 and CAMB version 1.5.0. The CLASS option was tested with CLASS version 3.3.0.)


### Installation instructions

Simply clone this repository and install the code with `pip`. Navigate to the directory in which you would like to place the `hdlike` directory (i.e., the directory where this `README.md` file is located), and then run the following commands:

```
git clone https://github.com/CMB-HD/hdlike.git
cd hdlike
pip install . --user
```

To uninstall the code, use the command `pip uninstall hdlike`. (Note that you may have to navigate away from the `hdlike` directory before running this command).


## Testing the likelihood

To test the likelihood, run the following command:

```
python test_hdlike.py
```

The test code above assumes the input mock CMB-HD data contains delensed CMB spectra with no baryonic feedback and no DESI BAO.  You can also add a flag to test the lensed CMB spectra, spectra including baryonic feedback effects, and/or include DESI BAO by running the following commands:

```
python test_hdlike.py --lensed
python test_hdlike.py --feedback
python test_hdlike.py --desi
```

You can use more than one flag, e.g. `python test_hdlike.py --desi --feedback`.

This will test that the code is calculating the correct likelihood by matching the output value to a precomputed likelihood value.

To test the likelihood with the theory calculated by CLASS instead of CAMB (see the [Using CLASS instead of CAMB](#using-class-instead-of-camb) section below), add the `--class` flag. This requires the `--lensed` flag, and cannot be used with `--feedback`:

```
python test_hdlike.py --lensed --class
python test_hdlike.py --lensed --class --desi
```

This uses the CLASS settings in `class_cobaya_settings.yaml`, and compares the theory to the CMB-HD bandpowers that were calculated with CLASS. Note that CLASS must be modified as described in the section below, or the test will fail.


### Testing the likelihood with Cobaya

To test the interface between hdlike and Cobaya, run the following command:

```
python test_hdlike_cobaya.py
```

The same flag options apply as above for `test_hdlike.py`.  This will also test that the code is outputting the correct likelihood by matching the output value to a precomputed likelihood.  In addition, it will test the Cobaya initialization.  With the `--class` flag, the test instead loads the settings from `example_cmbhd_class_lcdm_nrun_nnu_mnu.yaml` (see below), and checks that Cobaya can initialize the likelihood with `classy`.


## Using the likelihood with Cobaya 

We provide an example file that can be input into Cobaya and run, called `example_cmbhd.yaml`.  We also provide the four input files used for the CMB-HD MCMC chains in Cheslog et. al. (2026), which use lensed CMB-HD spectra with mock DESI BAO, for a $\Lambda$CDM model (`lcdm`) and a $\Lambda$CDM + $\alpha_\mathrm{s}$ + $N_\mathrm{eff}$ + $\sum m_\nu$ model (`lcdm_nrun_nnu_mnu`), with the theory calculated by either CAMB or CLASS:

```
example_cmbhd_camb_lcdm.yaml
example_cmbhd_camb_lcdm_nrun_nnu_mnu.yaml
example_cmbhd_class_lcdm.yaml
example_cmbhd_class_lcdm_nrun_nnu_mnu.yaml
```

In these four files, you must replace `/PATH/TO/` with the absolute path to your `hdlike` directory (the directory in which this README is located). Each one points at a proposal matrix from its own converged chain, in the `hdlike/data/proposal_cov/from_chains` directory.

To run the built-in Cobaya test of the initialization, use the following command:

```
cobaya-run example_cmbhd.yaml --test
```

Note that `test_hdlike_cobaya.py` loads information from the file `example_cmbhd.yaml` (or from `example_cmbhd_class_lcdm_nrun_nnu_mnu.yaml`, with the `--class` flag); if you’d like to modify the example files, please create your own copy.

MCMC chains can be run with the following command:

```
cobaya-run your_cmbhd.yaml
```

To run multiple chains with MPI you should follow the instructions for your machine.  For NERSC, we provide a job script template, named `nersc_perlmutter_job_template.sb` (see the NERSC section below) that can be run from this `hdlike` directory with the following command:

```
sbatch nersc_perlmutter_job_template.sb
```

(You should include the name of the CMB project to which you are charging your hours somewhere in the `job-name`; e.g., `#SBATCH --job-name=hdlike_CMBEXP`. Please see below for more information about using Cobaya on NERSC.)

To analyze the MCMC chains, we provide an example Jupyter notebook named `hdlike_cobaya_results.ipynb`. We also provide a set of pre-computed chains that can be used to test the notebook.

## Using CLASS instead of CAMB

The likelihood can also be used with the theory calculated by CLASS instead of CAMB, by using `classy` under the `theory` block of your Cobaya `.yaml` file and setting the option `use_class: True` under HDLike in the `likelihood` block:

```
likelihood:
	hdlike.hdlike.HDLike:
		delensed: False
		use_class: True
```

When `use_class` is `True`, the likelihood compares the theory to CMB-HD bandpowers that were calculated with CLASS (using the settings of Cheslog et. al. (2026)), rather than with CAMB, so that the mock data and the theory are calculated consistently. These bandpowers are provided with `hdlike` in the `hdlike/data` directory, and are only available for the latest version of the CMB-HD data (`hd_data_version: v1.2`), with lensed spectra and no baryonic feedback. If `baryonic_feedback` is `True`, the likelihood uses the bandpowers calculated with CAMB instead, and warns you that the feedback model in CLASS may not match. The likelihood will raise an error if `use_class` does not match the theory code Cobaya is using.

__Note__ the following before using CLASS:
- __The CLASS source code must be modified__ as described in Appendix A of Cheslog et. al. (2026), or the calculation will fail. In `source/lensing.c`, the variables `num_mu` and `index_mu` must be changed from `int` to `long long` (and `icount` to `unsigned long long`), so that the lensing calculation does not overflow at the multipoles used for CMB-HD. In `source/input.c`, the check that rejects a negative `N_ur` must be commented out, so that $N_\mathrm{eff}$ can be varied below its standard value with three massive neutrinos. These are at lines 124 and 2470 of CLASS version 3.3.4; the chains in Cheslog et. al. (2026) used version 3.3.0. A warning is issued as a reminder whenever `use_class` is `True`.
- Delensed spectra cannot be calculated with CLASS, so you must set `delensed: False`. The likelihood will raise an error otherwise.
- CLASS reads the BBN table named by `sBBN file` in the `classy` block itself, so it must be an absolute path. We provide the table used in Cheslog et. al. (2026), `hdlike/data/PRIMAT21_class_format.dat`; replace `/PATH/TO/` in the example files with the path to your `hdlike` directory.
- The parameter names in the `params` block are the CLASS names (e.g. `omega_b`, `omega_cdm`, `tau_reio`, `ln_A_s_1e10`, `n_s`), and $H_0$ is sampled in place of $\theta_\mathrm{MC}$. The sum of the neutrino masses and $N_\mathrm{eff}$ are sampled as `mnu` and `N_eff`, and converted to CLASS's `m_ncdm` and `N_ur` in the `params` block; see `class_cobaya_settings.yaml`.

The example files `example_cmbhd_class_lcdm.yaml` and `example_cmbhd_class_lcdm_nrun_nnu_mnu.yaml` (see above) contain the CLASS settings used in Cheslog et. al. (2026), and `class_cobaya_settings.yaml` holds the same settings for use with `generate_cobaya_input_file.py` (see below).

## Running MCMC chains on NERSC

A NERSC job script template, named `nersc_perlmutter_job_template.sb`, is also included. 

Some NERSC/Perlmutter-specific notes:
- Save your chains in your `SCRATCH` directory. This will involve modifying the `output` path in `your_cmbhd.yaml`. To find your `SCRATCH` directory, use the command `echo $SCRATCH`. However, be aware that this is not permanent storage: see the [NERSC documentation](https://docs.nersc.gov/filesystems/perlmutter-scratch/) for more details.
- You'll need to follow the instructions in the [NERSC documentation](https://docs.nersc.gov/development/languages/python/using-python-perlmutter/#mpi4py-on-perlmutter) to use MPI with `mpi4py`.
- You must provide the name of the allocation account being used on the line `#SBATCH --account=` (after the `=` sign; e.g. `#SBATCH --account=mp107` for CMB).
- You _should_ include the name of the CMB project you are charging your hours to somewhere in the `job-name` (e.g., `#SBATCH --job-name=hdlike_CMBEXP`). 
- If you're activating a conda environment within the job script, you need to use `source activate` instead of the `conda activate` used on the login nodes.

---

## Additional helpful information

Below we have provided some additional information (but you do not need to read further in order to use the CMB-HD likelihood as described above).


### Generating a new input `.yaml` file for Cobaya

We also provide a script that the user can run to generate their own `.yaml` file to input into Cobaya. To use this script, adjust the parameters of the `hdlike_settings.yaml` file. (One can also adjust the settings for CAMB and Cobaya in `camb_cobaya_settings.yaml`, for example to change which parameters are fixed/varied, such as baryonic feedback. One can also specify a proposal matrix in this file; we strongly recommend using a proposal matrix to speed up the convergence of the MCMC chains.) 

To make the new input file, run the following command:

```
python generate_cobaya_input_file.py hdlike_settings.yaml camb_cobaya_settings.yaml
```

To calculate the theory with CLASS instead, pass `class_cobaya_settings.yaml` in place of `camb_cobaya_settings.yaml`. This sets `use_class: True` for the likelihood, and requires `delensed: False` in `hdlike_settings.yaml`. (Note that you must still replace `/PATH/TO/` in the `sBBN file` entry of `class_cobaya_settings.yaml` with the path to your `hdlike` directory.)

### Note on using DESI BAO with Cobaya

You may see the following error message when trying to use the provided DESI BAO likelihood with Cobaya:

```
cobaya.component.ComponentNotInstalledError: The data for this likelihood has not been correctly installed. To install it, run `cobaya-install bao.generic`
```

If this occurs, running the `cobaya-install bao.generic` command will clone the Cobaya `bao_data` [repository](https://github.com/CobayaSampler/bao_data) and fix the issue. See the [Cobaya docs](https://cobaya.readthedocs.io/en/latest/installation_cosmo.html) for more information about the `cobaya-install` command.


### Using your own Cobaya configuration file

If you have your own Cobaya configuration `.yaml` file already, you can use the CMB-HD likelihood following the directions below.


#### To add the CMB-HD likelihood

To use the CMB-HD delensed $TT/TE/EE/BB$ CMB + $\kappa\kappa$ mock spectra and covariance matrix, you only need to mention HDLike in your Cobaya configuration `.yaml` file or python dictionary, under the `likelihood` block. For example, if running Cobaya from a `.yaml` file, your `likelihood` block would be

```
likelihood:
	hdlike.hdlike.HDLike:
```

To instead use CMB-HD lensed $TT/TE/EE/BB$ CMB + $\kappa\kappa$ mock spectra and covariance matrix, set the option `delensed: False` under HDLike in the `likelihood` block:

```
likelihood:
	hdlike.hdlike.HDLike:
		delensed: False
```

Additionally, you can change the maximum multipole for the CMB or lensing mock spectra and covariance matrix by setting the option `lmax` or `Lmax`, respectively, to an integer below the default `lmax = 20100` or `Lmax = 20100`. The minimum multipole for all spectra is fixed to 30. You can exclude the lensing power spectrum data by setting the option `use_cmb_lensing_spectrum: False`, or exclude the CMB $TT/TE/EE/BB$ data by setting the option `use_cmb_power_spectra: False`.

A basic example YAML file is provided, named `example_cmbhd.yaml`. The proposal widths of parameters in the `parameters` block were obtained from a CMB-HD Fisher matrix, located in the `hdlike/data/proposal_cov/from_fisher` directory. We also provide a few proposal matrices from MCMC runs in the `hdlike/data/proposal_cov/from_chains` directory, including those from the CMB-HD chains in Cheslog et. al. (2026), for both CAMB and CLASS.

- __Note__ that you should provide your own `output` path at the top of the file; the chain files that are output may be large. See the [Cobaya docs](https://cobaya.readthedocs.io/en/latest/output.html#output-shell) for details.


#### To add the DESI BAO likelihood

To combine CMB-HD with mock DESI BAO data, simply add the following lines to your `likelihood` block:

```
likelihood:
	hdlike.hdlike.HDLike:
		delensed: True
	bao.generic:
		measurements_file: /PATH/TO/hdlike/hdlike/data/mock_desi_bao_rs_over_DV_data.txt
		cov_file: /PATH/TO/hdlike/hdlike/data/mock_desi_bao_rs_over_DV_cov.txt
```

The example YAML file `example_cmbhd.yaml` contains a sample `bao.generic` block, but you must replace `/PATH/TO/` with the absolute path to your `hdlike` directory (the directory in which this README is located). See the [Cobaya docs](https://cobaya.readthedocs.io/en/latest/likelihood_bao.html) for more information about the generic BAO likelihood.

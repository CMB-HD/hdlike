# `HDLike`: Likelihood for CMB-HD

This is a mock CMB-HD likelihood including lensed and delensed $TT/TE/EE/BB$ CMB + lensing $\kappa\kappa$ spectra from multipoles 30 to 20,000.  We also include a mock DESI BAO likelihood. The likelihood can be used with Cobaya.  If you use this code, please cite:
- [MacInnis, Sehgal, and Rothermel (2023)](https://arxiv.org/abs/2309.03021)
- If you use CLASS, please also cite [Cheslog, Finson, MacInnis, Sehgal, Afshordi, Nerval, and Hložek (2026)](https://arxiv.org/abs/2609.35964)
- If you use the CMB-HD mock data, please also cite the appropriate references given in the [`hdMockData` repository](https://github.com/CMB-HD/hdMockData#forecasting-data-for-cmb-hd) for the mock data version used (the latest version by default)


## Installation of likelihood


### Installation requirements

To use the CMB-HD likelihood, you must install Python version >= 3, [NumPy](https://numpy.org/), and [hdMockData](https://github.com/CMB-HD/hdMockData).

To test the likelihood by running `test_hdlike.py`, you must also install [CAMB](https://camb.readthedocs.io/en/latest/) or [CLASS](https://github.com/lesgourg/class_public). **If you are using CLASS, see [this note](#using-class)** below.

To use the likelihood with Cobaya, you must have Python version >= 3.8 and [Cobaya](https://cobaya.readthedocs.io/en/latest/index.html).

(The likelihood has been tested with Cobaya version 3.5.4, CAMB version 1.5.4, and CLASS version 3.3.4)


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

The command above will use the latest mock CMB-HD data including delensed CMB power spectra and the lensing convergence power spectrum, without mock DESI BAO data, and use CAMB for the theory calculation, for a non-linear model without baryonic feedback.  You can also add a flag to test the lensed CMB spectra, spectra including baryonic feedback effects, include mock DESI BAO data, and/or use CLASS instead of CAMB by running the following commands:

```
python test_hdlike.py --lensed
python test_hdlike.py --feedback
python test_hdlike.py --desi
python test_hdlike.py --use_class --lensed
```

You can use more than one flag, e.g. `python test_hdlike.py --desi --feedback`, but you **must** always pass `--lensed` when you are using CLASS, since CLASS does not calculate delensed CMB power spectra.

This will test that the code is calculating the correct likelihood by matching the output value to a precomputed likelihood value.


### Testing the likelihood with Cobaya

To test the interface between hdlike and Cobaya, run the following command:

```
python test_hdlike_cobaya.py
```

The same flag options apply as above for `test_hdlike.py`.  This will also test that the code is outputting the correct likelihood by matching the output value to a precomputed likelihood.  In addition, it will test the Cobaya initialization.


## Using the likelihood with Cobaya 

We provide two example files that can be input into Cobaya and run, called `example_cmbhd_CAMB.yaml` to use with CAMB, or `example_cmbhd_CLASS.yaml` to use with CLASS.  

To run the built-in Cobaya test of the initialization with CAMB or CLASS, use the following command(s):

```
cobaya-run example_cmbhd_CAMB.yaml --test
cobaya-run example_cmbhd_CLASS.yaml --test
```

Note that `test_hdlike_cobaya.py` loads information from these example YAML files; if you’d like to modify the example files, please create your own copies.

MCMC chains can be run with the following command:

```
cobaya-run your_cmbhd.yaml
```

To run multiple chains with MPI you should follow the instructions for your machine  (see, e.g., the [NERSC](#running-mcmc-chains-on-nersc) section below).  

To analyze the MCMC chains, we provide an example Jupyter notebook named `hdlike_cobaya_results.ipynb`. We also provide a set of pre-computed chains that can be used to test the notebook.


---

## Using CLASS

Note that CLASS does not calculate delensed CMB power spectra, but the default in `hdlike` is to use the delensed CMB-HD mock data, *and* to calculate the delensed theory prediction for a given sample of cosmological parameters. **If you are using CLASS, you must always use lensed CMB power spectra** by setting `delensed: False` under `hdlike.hdlike.HDLike` in the `likelihood` block of your YAML file.

### CLASS modifications required to use the fiducial CMB-HD CLASS parameters

**CLASS must be modified** in order to use the CMB-HD CLASS accuracy settings in `example_cmbhd_CLASS.yaml`, and in order to vary the effective number of relativistic species. These modifications must be made *before* compiling and installing CLASS (e.g., via `make` or `pip install`). 

The instructions are given in Appendix A of Cheslog et. al. (2026) and repeated below. We will use `$CLASS_DIR` for the path to the CLASS repository directory (named `class_public` by default); i.e., the directory that contains the CLASS "readme" and `explanatory.ini` files.

- In `$CLASS_DIR/source/lensing.c` ([line 124](https://github.com/lesgourg/class_public/blob/v3.3.4/source/lensing.c#L124) in version 3.3.4, or [lines 130-131](https://github.com/lesgourg/class_public/blob/v3.4.0/source/lensing.c#L130) in version 3.4.0), you must update the type of `num_mu` and `index_mu` to be `long long`, and update the type of `icount` to be `unsigned long long`. After these modifications, the relevant lines in the file should be:
    ```c
    
    unsigned long long icount;
    long long num_mu , index_mu;
    
    
    ```
  This is necessary in order to use `l_max_scalars` above approximately 14,000 when `accurate_lensing=1`.

- In `$CLASS_DIR/source/input.c`, you must comment out the following line ([line 2470](https://github.com/lesgourg/class_public/blob/v3.3.4/source/input.c#L2470) in version 3.3.4 or [line 2548](https://github.com/lesgourg/class_public/blob/v3.4.0/source/input.c#L2548) in version 3.4.0):
    ```c
    
    /*
    class_test(pba->Omega0_ur<0,errmsg,"You cannot set the density of ultra-relativistic relics (dark radiation/neutrinos) to negative values. You might have input a total Neff smaller than what your massive neutrinos require minimally (around 1.02 * N_ncdm * deg_ncdm).");
    */
    
    
    ```
  This is necessary in order to vary `Neff` below approximately 3.0396 (or, equivalently, the `N_ur` parameter below zero) with three massive neutrinos (`N_ncdm=3`). 
  
- In order to use the `sBBN file` in `example_cmbhd_CLASS.yaml` (available [here](https://github.com/CMB-HD/hdMockData/blob/main/hd_mock_data/data/theory/PRIMAT_Yp_DH_ErrorMC_2021_CLASS.dat), which is a copy of the corresponding file provided by [CAMB](https://github.com/cmbant/CAMB/blob/master/camb/PRIMAT_Yp_DH_ErrorMC_2021.dat) that has been formatted for CLASS), you must place a copy of that file in the `$CLASS_DIR/external/bbn/` directory. You may copy it over yourself, or follow these instructions:

  First, make sure that `hdMockData` has been correctly installed by running the command `python -c "from hd_mock_data import hd_data` ; if nothing happens, you may proceed; otherwise, follow the [instructions](https://github.com/CMB-HD/hdMockData/tree/main#installation) to install `hd_mock_data`. 
  
  Then, run the following two commands from within your `$CLASS_DIR`:
  
  ```bash
  CLASS_SBBN_FILE=$(python -c "from hd_mock_data import hd_data; print(hd_data.class_sbbn_file())")
  
  ```
  
  ```bash
  cp $CLASS_SBBN_FILE external/bbn/
  ```


## Note about CAMB versions

The fiducial set of CAMB parameters used by CMB-HD includes a parameter named `lens_output_margin` starting in CAMB version 2.0.0, or named `lens_margin` in previous versions. When you test the likelihood or generate new YAML files using the python scripts provided here, we attempt to import CAMB in order to determine which CAMB version you are using, and use the appropriate parameter name for that version. If CAMB cannot be imported, we default to the newer `lens_output_margin` name; this is the name used in the `example_cmbhd_CAMB.yaml` file.



---

## Additional helpful information

Below we have provided some additional information (but you do not need to read further in order to use the CMB-HD likelihood as described above).


### Generating a new input `.yaml` file for Cobaya

We also provide a Python script, `generate_cobaya_input_file.py`, that you can run to generate a new `.yaml` file to input into Cobaya. 

By default, running `python generate_cobaya_input_file.py` will generate a Cobaya YAML file that uses the CMB-HD likelihood with the latest delensed CMB and CMB lensing mock data, uses CAMB for the theory calculation with the fiducial CMB-HD accuracy settings, and varies the six $\Lambda\mathrm{CDM}$ parameters $\left(\Omega_\mathrm{b} h^2,\, \Omega_\mathrm{c} h^2,\, 100\theta_\mathrm{MC}\, \ln\left(10^{10} A_\mathrm{s}\right),\, n_\mathrm{s},\, \tau \right)$ with a prior of $\sigma(\tau) = 0.005$ using the [`mcmc` sampler](https://cobaya.readthedocs.io/en/latest/sampler_mcmc.html).

These defaults may be changed by passing different options to `generate_cobaya_input_file.py`; run `python generate_cobaya_input_file.py --help` to see the full list of options. You may also pass an existing YAML file to modify (a copy will be saved). 

For example, to update the `example_cmbhd_CAMB.yaml` file provided here to use lensed CMB spectra, include mock DESI BAO data (with the correct paths to the data files), and use a non-linear model with baryonic feedback effects, and use a new `output` path, you would use the following command:

```bash
python generate_cobaya_input_file.py example_cmbhd_CAMB.yaml -o hdlike_chains/cmbhd_lensed_bao_feedback --lensed --bao --feedback
```

Note that you should always check and test the saved YAML file, and modify it if necessary, before running MCMC chains.

You may also use one of the pre-computed proposal matrices for the parameter covariance (the [`covmat` option](https://cobaya.readthedocs.io/en/latest/sampler_mcmc.html#covariance-matrix-of-the-proposal-pdf) of the `mcmc` sampler) provided with `hdlike` by passing `--use_param_covmat`. These matrices are derived from chains run with `hdlike`. Note that there is not a proposal matrix available for each combination of data, theory model, and set of varied parameters, but these should be a good starting point. 


### Note on using DESI BAO with Cobaya

You may see the following error message when trying to use the provided DESI BAO likelihood with Cobaya:

```
cobaya.component.ComponentNotInstalledError: The data for this likelihood has not been correctly installed. To install it, run `cobaya-install bao.generic`
```

If this occurs, running the `cobaya-install bao.generic` command will clone the Cobaya `bao_data` [repository](https://github.com/CobayaSampler/bao_data) and fix the issue. See the [Cobaya documentation](https://cobaya.readthedocs.io/en/latest/installation_cosmo.html) for more information about the `cobaya-install` command.


### Using your own Cobaya configuration file

If you have your own Cobaya configuration `.yaml` file already, you can use the CMB-HD likelihood following the directions below.


#### To add the CMB-HD likelihood

To use the CMB-HD delensed $TT/TE/EE/BB$ CMB + $\kappa\kappa$ mock spectra and covariance matrix, you only need to mention HDLike in your Cobaya configuration `.yaml` file or python dictionary, under the `likelihood` block. For example, if running Cobaya from a `.yaml` file, your `likelihood` block would be

```yaml
likelihood:
	hdlike.hdlike.HDLike:
```

To instead use CMB-HD lensed $TT/TE/EE/BB$ CMB + $\kappa\kappa$ mock spectra and covariance matrix, set the option `delensed: False` under HDLike in the `likelihood` block; **note** that this is *required* if you are using CLASS:

```yaml
likelihood:
	hdlike.hdlike.HDLike:
		delensed: False
```

Additionally, you can change the maximum multipole for the CMB or lensing mock spectra and covariance matrix by setting the option `lmax` or `Lmax`, respectively, to an integer below the default `lmax = 20100` or `Lmax = 20100`. The minimum multipole for all spectra is fixed to 30. You can exclude the lensing power spectrum data by setting the option `use_cmb_lensing_spectrum: False`, or exclude the CMB $TT/TE/EE/BB$ data by setting the option `use_cmb_power_spectra: False`.

Two basic example YAML files are provided, named `example_cmbhd_CAMB.yaml` and `example_cmbhd_CLASS.yaml`. The proposal widths of parameters in the `parameters` block were obtained from a CMB-HD Fisher matrix, located in the `hdlike/data/proposal_cov/from_fisher` directory. We also provide a few proposal matrices from MCMC runs in the `hdlike/data/proposal_cov/from_chains` directory.

- __Note__ that you should provide your own `output` path at the top of the file; the chain files that are output may be large. See the [Cobaya documentation](https://cobaya.readthedocs.io/en/latest/output.html#output-shell) for details.


#### To add the DESI BAO likelihood

To combine CMB-HD with mock DESI BAO data, simply add the following lines to your `likelihood` block:

```yaml
likelihood:
	hdlike.hdlike.HDLike:
		delensed: True
	bao.generic:
		measurements_file: /PATH/TO/hdlike/hdlike/data/mock_desi_bao_rs_over_DV_data.txt
		cov_file: /PATH/TO/hdlike/hdlike/data/mock_desi_bao_rs_over_DV_cov.txt
```

The example YAML files `example_cmbhd_CAMB.yaml` and `example_cmbhd_CLASS.yaml` contain a sample `bao.generic` block, but you must replace `/PATH/TO/` with the absolute path to your `hdlike` directory (the directory in which this README is located). See the [Cobaya documentation](https://cobaya.readthedocs.io/en/latest/likelihood_bao.html) for more information about the generic BAO likelihood.


### Running MCMC chains on NERSC

(Note: this section was last updated in 2023; refer to the NERSC documentation for more recent information.)

A NERSC job script template, named `nersc_perlmutter_job_template.sb`, is also included. 

Some NERSC/Perlmutter-specific notes:
- Save your chains in your `SCRATCH` directory. This will involve modifying the `output` path in `your_cmbhd.yaml`. To find your `SCRATCH` directory, use the command `echo $SCRATCH`. However, be aware that this is not permanent storage: see the [NERSC documentation](https://docs.nersc.gov/filesystems/perlmutter-scratch/) for more details.
- You'll need to follow the instructions in the [NERSC documentation](https://docs.nersc.gov/development/languages/python/using-python-perlmutter/#mpi4py-on-perlmutter) to use MPI with `mpi4py`.
- You must provide the name of the allocation account being used on the line `#SBATCH --account=` (after the `=` sign; e.g. `#SBATCH --account=mp107` for CMB).
- You _should_ include the name of the CMB project you are charging your hours to somewhere in the `job-name` (e.g., `#SBATCH --job-name=hdlike_CMBEXP`). 
- If you're activating a conda environment within the job script, you need to use `source activate` instead of the `conda activate` used on the login nodes.



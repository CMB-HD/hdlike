# Readme for HDLike

This is a mock CMB-HD likelihood including lensed and delensed $TT/TE/EE/BB$ CMB + lensing $\kappa\kappa$ spectra from multipoles 30 to 20,000.  We also include a mock DESI BAO likelihood. The likelihood can be used with Cobaya.  Please cite MacInnis, Sehgal, and Rothermel (2023).


## Installation of likelihood


### Installation requirements

To use the CMB-HD likelihood, you must install Python version >= 3 and [NumPy](https://numpy.org/).

To test the likelihood by running `test_hdlike.py`, you must also install [CAMB](https://camb.readthedocs.io/en/latest/).

To use the likelihood with Cobaya, you must have Python version >= 3.8 and [Cobaya](https://cobaya.readthedocs.io/en/latest/index.html).

(The likelihood has been tested with Cobaya version 3.3.2 and CAMB version 1.5.0.)


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


### Testing the likelihood with Cobaya

To test the interface between hdlike and Cobaya, run the following command:

```
python test_hdlike_cobaya.py
```

The same flag options apply as above for `test_hdlike.py`.  This will also test that the code is outputting the correct likelihood by matching the output value to a precomputed likelihood.  In addition, it will test the Cobaya initialization.


## Using the likelihood with Cobaya 

We provide an example file that can be input into Cobaya and run, called `example_cmbhd.yaml`.  

To run the built-in Cobaya test of the initialization, use the following command:

```
cobaya-run example_cmbhd.yaml --test
```

Note that `test_hdlike_cobaya.py` loads information from the file `example_cmbhd.yaml`; if you’d like to modify the example file, please create your own copy.

MCMC chains can be run with the following command:

```
cobaya-run example_cmbhd.yaml
```

To run multiple chains with MPI you should follow the instructions for your machine.  For NERSC, we provide a job script template, named `nersc_perlmutter_job_template.sb` (see the NERSC section below) that can be run from this `hdlike` directory with the following command:

```
sbatch nersc_perlmutter_job_template.sb
```

(You should include the name of the CMB project to which you are charging your hours somewhere in the `job-name`; e.g., `#SBATCH --job-name=hdlike_CMBEXP`. Please see below for more information about using Cobaya on NERSC.)


---

## Additional helpful information

Below we have provided some additional information (but you do not need to read further in order to use the CMB-HD likelihood as described above).


### Generating a new input `.yaml` file for Cobaya

We also provide a script that the user can run to generate their own `.yaml` file to input into Cobaya. To use this script, adjust the parameters of the `hdlike_settings.yaml` file. (One can also adjust the settings for CAMB and Cobaya in `camb_cobaya_settings.yaml`, for example to change which parameters are fixed/varied, such as baryonic feedback. One can also specify a proposal matrix in this file.) 

To make the new input file, run the following command:

```
python generate_cobaya_input_file.py hdlike_settings.yaml camb_cobaya_settings.yaml
```

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

A basic example YAML file is provided, named `example_cmbhd.yaml`. The proposal widths of parameters in the `parameters` block were obtained from a CMB-HD Fisher matrix, located in the `hdlike/data/proposal_cov/from_fisher` directory. We also provide a few proposal matrices from MCMC runs in the `hdlike/data/proposal_cov/from_chains` directory.

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


### Running MCMC chains on NERSC with Cobaya

A NERSC job script template, named `nersc_perlmutter_job_template.sb`, is also included.

Some NERSC/Perlmutter-specific notes:
- You must provide the name of the allocation account being used on the line `#SBATCH --account=` (after the `=` sign; e.g. `#SBATCH --account=mp107` for CMB).
- You _should_ include the name of the CMB project you are charging your hours to somewhere in the `job-name` (e.g., `#SBATCH --job-name=hdlike_CMBEXP`). 
- If you're activating a conda environment within the job script, you need to use `source activate` instead of the `conda activate` used on the login nodes.
- You'll need to follow the instructions in the [NERSC documentation](https://docs.nersc.gov/development/languages/python/using-python-perlmutter/#mpi4py-on-perlmutter) to use MPI with `mpi4py`.
- You'll likely want to save your chains in your `SCRATCH` directory, but be aware that this is not permanent storage: see the [NERSC documentation](https://docs.nersc.gov/filesystems/perlmutter-scratch/) for more details.

# __The latest version of the CMB-HD mock data is located in the `hdMockData` [repository](https://github.com/CMB-HD/hdMockData)__.

By default, the latest CMB-HD mock data is used in the likelihood.

The `hdlike/data` directory previously contained  the original version `'v1.0'` of the CMB-HD data used in MacInnis et. al. (2023). Those files, along with the latest versions, are now located in the `hdMockData` repository [here](https://github.com/CMB-HD/hdMockData). We also provide Python functions with `hdMockData` that can be used to easily access the files.

Below, we describe the data files used in the likelihood, and the methods within the the `HDMockData` class of the `hdMockData` repository that can be used to access the data. See the `hdMockData` [repository](https://github.com/CMB-HD/hdMockData) for more information. 

---

The data needed for the likelihood includes the mock CMB-HD lensed/delensed theory spectra, lensed/delensed covariance matrices, lensing power spectrum noise (for delensed spectra at each step of the chain), and a binning file, along with the mock DESI BAO data.


- combined 90 and 150 GHz,
- lensing spectra & covmat blocks from `lmin=30` to `Lmax=20100`, 
- CMB spectra & covmat blocks binned from `lmin=30` to `lmax=20100`, 
- including residual extragalactic foregrounds in TT noise,
- using Advanced Simons Observatory noise for BB below $\ell = 1000$,
- lensing reconstruction noise is a MV combination of Hu & Okamoto TT, TE, EE, EB, TB calculated iteratively (Hotinli et. al., arXiv:2111.15036) below L ~ 5000, and the same Hu & Okamoto TE, EE, EB, TB + simulation-based HDV TT estimator above L ~ 5000 (from Han et. al. 2022, arXiv:2112.02109).

CMB-HD data descriptions and format of the data files: 
1. Binning file: Contains the bin edges used to bin the data. See the `bin_edges` method of the `HDMockData` class.
2. Binned theory spectra: The binned lensed or delensed CMB TT, TE, EE, BB and lensing convergence kappakappa power spectra. See the `mcmc_bandpowers` method of the `HDMockData` class.
  - Each file contains a single column of data with length `5 * nbin` holding the binned theory spectra in the order `TT`, `TE`, `EE`, `BB`, `kk` (lensing power), where `nbin` is the number of bins for a _single_ binned spectrum (e.g., TT). The first `nbin` elements are the binned TT spectra, the next `nbin` elements are the binned TE spectra, etc.
  - The CMB spectra are in units of uK^2 in $C_l$'s (i.e., no $l(l+1)/2\pi$ factors). The lensing convergence power spectrum is in the form $C_L^{\kappa\kappa} = [L (L + 1)]^2 C_L^{\phi\phi} / 4$, (i.e., we do not follow the CAMB convention with $2\pi$ in the denominator) where $C_L^{\phi\phi}$ is the lensing potential power spectrum. 
3. Binned, block covariance matrices: The full covariance matrix for the binned spectra, containing 25 blocks for the auto- and cross-covariance between the different spectra. It has a shape `(5 * nbin, 5 * nbin)`, with the blocks in the same order and units as the spectra described above. See the `block_covmat` method of the `HDMockData` class.
4. Lensing noise: Only used for __delensed__ spectra. This file contains two columns; the second is the _unbinned_ lensing reconstruction noise for the lensing convergence data, and the first is the corresponding multipoles. The noise follows the same convention as the lensing convergence spectrum, i.e. $N_L^{\kappa\kappa} = [L (L + 1)]^2 N_L^{\phi\phi} / 4$. See the `lensing_noise_spectrum` method of the `HDMockData` class.


Note: the mock CMB and lensing data was made with higher CAMB accuracy than the default - see the example Cobaya `.yaml` file for the accuracy settings used.


Mock DESI BAO data files:
- `mock_desi_bao_rs_over_DV_data.txt` and `mock_desi_bao_rs_over_DV_cov.txt` : The first file contains the mock DESI BAO data; the first column gives the redshift, and the second gives the quantity $r_s/d_V(z)$, where $d_V(z) \equiv \left[(1+z)d_A(z)\right]^{2/3} \left[cz / H(z)\right]^{1/3}$, $d_A(z)$ is the angular diameter distance, and $H(z)$ is the expansion rate. The second file contains the (diagonal) covariance matrix for this data, derived from Tables 2.3 and 2.5 in arXiv:1611.00036.


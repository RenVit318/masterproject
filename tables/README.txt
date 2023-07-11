The csv files in this directory contain the results of the work done by R. Kievit for his Master Thesis in coorporation with Dr. J. de Bruijne (ESA) and Dr. A.G.A. Brown (Leiden University) titled "Connecting Mother and Daughter: Creating an Ideal Hipparcos - Gaia Crossmatch". For the scientific context behind this work a copy of the thesis can be requested through correspondence with the main author at rens@renskievit.com. 

All files are ASCII encoded, each row is an entry and each column is separated by a tab. Strings at the top of each column give the corresponding name. This file contains a brief overview of each table and the contents therin.


 > Hipparcos_Mix.csv

Linear combination of the original Hipparcos reduction (ESA, 1997) and a later reduction (van Leeuwen, 2007) with respective weights 0.4 and 0.6 as prescribed by Brandt (2021). All columns are named like the van Leeuwen, (2007) version of the catalogue except 'ra_dec_corr' which is computed from the u_ij elements in their catalogue using the method described in Appendix B of Michalik et al. (2014). Additionaly 0.6 mas and 0.25 mas/yr are added to the position and proper motion uncertainties respectively following Kovaevsky et al. (1997). 


 > GaiaBaseCat.csv

Data from Gaia DR3 (Gaia Collaboration et al., 2023) with numerous adjustments, which are (in order of application):
 - Only all sources with G < 14 are selected.
 - The reference frame correction by Cantat-Gaudin & Brandt (2021) is applied.
 - The zero-point parallax bias is corrected following the recipe by Lindegren et al. (2021).
 - All astrometric uncertainties are mutliplied by a factor 1.37 (Brandt, 2021).
 - The astrometric position and uncertainties are backpropagated from the Gaia epoch (J2016.0) to the Hipparcos epoch (J1991.25) using the EPOCH_PROP and EPOCH_PROP_UNCERTAINTIY functions on the Gaia archive. 

The table contains both the data in the original- and the Hipparcos epochs. All column names follow the Gaia DR3 notation, but the backpropagated data is indicated with a '_prop' suffix. 


 > GMM_params.csv

Weights, means and covariance of the sklearn Gaussian Mixture Model fitted on the five features given in the xm_results.csv table with K=9 and random_state=42. Following the definition of the GMM likelihood given in the thesis, the column names represent:
 - alpha_k   : Weight of component k
 - mu_ki     : The mean of component k in feature/dimension i
 - Sigma_kij : The covariance of component k between features i and j


 > xm_results.csv

The ultimate crossmatch result of this work. This table contains all Hipparcos - Gaia pairs (present in Hipparcos_Mix and GaiaBaseCat) that are within 10" of each other in J1991.25. For each of the pairs the five features are given (defined in Table 3.2), and also
 - P_ki : The probability that this Hipparcos - Gaia pair belongs to GMM component i
 - ncomp: The number of Gaia objects assigned to a single Hipparcos objects based on a photometric analysis
 - comp_idx: If ncomp > 1 this indicates which Gaia source is the primary, which the secondary etc. This classification is mostly based on magnitude with the brightest source usually being the primary.
 - xm_flag: Describes the corresponding neighbour pair based on:
    0. Not a believable match
    1. Believable crossmatch based on the GMM results
    2. Secondary component based on photometry (\autoref{sec:binaries})
    3. Believable crossmatch based on high $\sigma_{\mu}/\mu$ in component 3 (\autoref{sec:green})
    4. Believable crossmatch based on forward propagation in component 7 (\autoref{sec:pink})
    5. Crossmatch determined based on individual analysis

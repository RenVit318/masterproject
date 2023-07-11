############################
#
# This file contains all the functions that are used to create
# the final data products of this thesis which are stored in
#   masterproject/tables
#
############################

import numpy as np
from astropy.io import fits
from astropy.table import Table
import pickle

TPATH = '../tables/'
DPATH = '../../data/'
RPATH = '../results/'


def hipp_mix():
    hdul = fits.open(DPATH+'Hipparcos_mix.fits', memmap=True) # Necessary for large files
    table = Table(hdul[1].data)

    # Non-position astrometry is given to two digits in Hip2
    # Photometric accuracies given here are the same as in the Archive
    round_to_prec = ['pm_ra', 'pm_de', 'plx', 'e_pm_ra', 'e_pm_de', 'e_plx', 'ra_dec_corr']
    for name in round_to_prec:
        # np very slightly less accurate, but built-in function only works on single floats    
        # Hit in accuracy is not a problem at the level we are looking at
        table[name] = np.round(table[name].data, 2) 

    table.write(TPATH+'Hipparcos_Mix.csv', delimiter='\t', format='ascii', overwrite=True)
    

def gaia_base_cat():
    hdul = fits.open(DPATH+'GaiaBaseCat_Final.fits')#, memmap=True) # memmap for large files
    table = Table(hdul[1].data)

    table.write(TPATH+'GaiaBaseCat.csv', delimiter='\t', format='ascii', overwrite=True)


def gmm_weights():
    # Open GMM 
    with open(RPATH+'gmm_M5_K9_noadj.pickle', 'rb') as f:
        gmm = pickle.load(f)

    K = gmm.n_components # Number of components
    M = len(gmm.means_[0]) # Number of features

    # Create the csv file column by column 
    table_columns = []

    # Could either be 0 - 8 or 1 - 9
    table_columns.append(fits.Column(name='k', format='I', array=np.arange(K)))

    # Component Weights
    table_columns.append(fits.Column(name='alpha_k', format='D', array=gmm.weights_))

    # Component Means
    for i in range(M):
        means_k = gmm.means_[:, i]
        table_columns.append(fits.Column(name=f'mu_k{i}', format='D', array=means_k))

    # Component Covariances
    for i in range(M):
        for j in range(i, M):
            covar_k = gmm.covariances_[:, i, j]
            table_columns.append(fits.Column(name=f'Sigma_k{i}{j}', format='D', array=covar_k))
            
    # Write to file
    table_hdu = fits.BinTableHDU.from_columns(table_columns)

    table = Table(table_hdu.data)
    table.write(TPATH+'GMM_params.csv', delimiter='\t', format='ascii',  overwrite=True)


def xm_results():
    # Open table containing IDs and Feature data
    ml_data = fits.open(RPATH+'all_neighbours_ml_data_final.fits')[1].data
    probs = np.loadtxt(RPATH+'gmm_probs.txt')
    binaries = np.loadtxt(RPATH+'binary_labels.txt')
    flag = np.loadtxt(RPATH+'xm_flag.txt')

    # Create the csv file column by column
    table_columns = []
    
    # Hipparcos and Gaia IDs for crossmatching purpose
    table_columns.append(fits.Column(name=f'Hip_ID', format='K', array=ml_data['hip_id']))
    table_columns.append(fits.Column(name=f'Gaia_ID', format='K', array=ml_data['gaia source_id']))
    
    # Five Features
    features = ['G_Gpred', 'D', 'delta_pm_tot', 'delta_pm_angle', 'delta_plx']
    for feat in features:
        table_columns.append(fits.Column(name=feat, format='D', array=ml_data[feat]))

    # Nine base component probabilities
    for i in range(probs.shape[1]):
        table_columns.append(fits.Column(name=f'P_k{i}', format='D', array=probs[:,i]))

    # Binary indexes
    table_columns.append(fits.Column(name='ncomp', format='I', array=binaries[0]))
    table_columns.append(fits.Column(name='comp_idx', format='I', array=binaries[1]))

    # XM Flag
    table_columns.append(fits.Column(name='xm_flag', format='I', array=flag))

    # write to file
    table_hdu = fits.BinTableHDU.from_columns(table_columns)
    table = Table(table_hdu.data)
    table.write(TPATH+'xm_results.csv', delimiter='\t', format='ascii', overwrite=True)


def main():
    #hipp_mix()
    gaia_base_cat()
    #gmm_weights()
    #xm_results()

if __name__ == '__main__':
    main()

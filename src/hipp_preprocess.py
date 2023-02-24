##################
#
# Python scripts to preprocess the Hipparcos catalogue for the Gaia-Hipparcos crossmatch
# In this file:
#   - Covariance extraction code
#       - Transforms the U matrix as given by van Leeuwen+07 in Hipparcos 2 into a usable covariance matrix form
#   - Hipparcos mixing
#       - Mixes Hipparcos 1 and 2 with the 60/40 split prescribed by Brandt+21
#
##################

import numpy as np
from astropy.io import fits

def mix_hipparcos():
    pass

def extract_cov(tab):
    """Transform the u-values given in Hipparcos-2 (van Leeuwen, 2007) into a complete covariance matrix
       This code is based on a Python 2 implementation of the same procedure by A.G.A. Brown"""

    num = tab.shape[0]
    print(f'Converting U -> \Sigma for {num} Hip2 objects')

    # Setup arrays and parameters
    U = np.zeros((num,5,5))
    sigma = np.zeros((num,5))
    nobs = tab['Ntr']
    gof = tab['f2']
    nu = nobs-5
    Q = nu * (np.sqrt(2./(9.*nu)) * gof + 1 - (2./(9.*nu)))**3.
    u = np.sqrt(Q/nu) # confusing naming. Change?
    u_sq = u**2.

    # Populate sigma
    sigma_names = ['e_ra_rad', 'e_de_rad', 'e_plx', 'e_pm_ra', 'e_pm_de']
    for i, field in enumerate(sigma_names):
       sigma[:,i] = tab[field]

    # Populate U
    u_index = 1
    for j in range(5):
        for i in range(j+1):
            U[:, i, j] = tab[f'u{u_index}']
            if i == j:
                U[:, i, j] *= (u/sigma[:,i])
            u_index += 1

    U_transpose = np.transpose(U, axes=(0, 2, 1))
    # In case of a 3D array, a matmul treats it as a stack of matrices and treats
    # the first index as the stack index. 
    # N = Normal matrix (notation from Michalik et al. (2014))
    N = np.matmul(U_transpose, U)
    N_inverse = np.linalg.inv(N) 
    covar = np.reshape(u_sq[:, None, None] * N_inverse, N_inverse.shape)

    
    print(len(covar[:,0,1][np.abs(covar[:,0,1]) > 1e2]))
    idx = np.argmax(covar[:, 0, 1][covar[:, 0, 1] < 1e10])
    print(tab['hip'][idx])
    print(sigma[idx])
    print(gof[idx], nobs[idx])
    print(Q[idx], u[idx])
    print(covar[idx])
    print(U_transpose[idx])
    print(N[idx])
    print(N_inverse[idx])
    
    return covar  



def main():
    hipp_hdul = fits.open('../../data/Hipparcos2.fits')
    covar = extract_cov(hipp_hdul[1].data)

if __name__ == '__main__':
    main()

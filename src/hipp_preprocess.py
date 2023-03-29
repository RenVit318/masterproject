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
    """Combine the Hipparcos-1 and Hipparcos-2 catalogues based on the mixing factor from Brandt+21"""
    hip1 = fits.open('../../data/Hipparcos1.fits')[1].data
    hip2 = fits.open('../../data/Hipparcos2.fits')[1].data
    mix_factor = 0.6  # Brandt (2021)

    mix_columns = []
    hip1_mix = ['ra', 'de', 'pmra', 'pmde', 'plx', 'e_radeg', 'e_dedeg', 'e_pmra', 'e_pmde', 'e_plx', 'dera']
    hip2_mix = ['ra', 'dec', 'pm_ra', 'pm_de', 'plx', 'e_ra_rad', 'e_de_rad', 'e_pm_ra', 'e_pm_de', 'e_plx', 'ra_dec_corr']
    extra_data = ['hp_mag', 'b_v', 'v_i']

    print('Starting Hip1 Indexing')
    # Match indices on hip because hip1 has 263 more objects than hip2
    hip1_idx = np.full(hip1.shape[0], True)
    num_unmatched = 0
    for i in range(hip1.shape[0]):
        if hip1['hip'][i] not in hip2['hip']:
            hip1_idx[i] = False
            num_unmatched += 1
            print(f'No match found for hip1={hip1["hip"][i]}')
    print(f'Could not find a match for {num_unmatched} Hipparcos-1 objects')
    hip2_covar = extract_cov(hip2)[:, 0, 1]

    # Add in hip-id
    mix_columns.append(fits.Column(name='hip', format='K', array=hip2['hip']))  # Format K = 64bit int

    print('Starting Astrometric Mixing')
    # Mix the astrometric parameters
    for i in range(5):
        print(hip1_mix[i], hip2_mix[i])
        mixed_data = mix_factor * hip2[hip2_mix[i]] + (1 - mix_factor) * hip1[hip1_mix[i]][hip1_idx]
        mix_columns.append(fits.Column(name=hip2_mix[i], format='D', array=mixed_data))  # Format D = double precision float

    print('next step')
    # Mix uncertainties on astrometric parameters
    for i in range(5, 10):
        print(hip1_mix[i], hip2_mix[i])
        mixed_data = np.sqrt(mix_factor*hip2[hip2_mix[i]]**2 + (1-mix_factor) * hip1[hip1_mix[i]][hip1_idx]**2)

        if hip2_mix[i] == 'e_ra_rad':
            mixed_data = np.sqrt(mixed_data**2 + 0.6**2) # Reference frame correction
            e_ra_mix = mixed_data
        if hip2_mix[i] == 'e_de_rad':            
            mixed_data = np.sqrt(mixed_data**2 + 0.6**2)
            e_de_mix = mixed_data
        if hip2_mix[i] == 'e_pm_ra' or hip2_mix[i] == 'e_pm_de':
            mixed_data = np.sqrt(mixed_data**2 + 0.25**2)

        mix_columns.append(fits.Column(name=hip2_mix[i], format='D', array=mixed_data))


    # Mix the correlations via the covariance matrix
    cov_radec_hip1 = hip1['dera'] * hip1['e_radeg'] * hip1['e_dedeg']
    #cov_radec_hip2 = hip2['ra_dec_corr'] * hip2['e_ra_rad'] * hip2['e_de_rad']
    cov_radec_mix = mix_factor * hip2_covar + (1 - mix_factor) * cov_radec_hip1[hip1_idx]
    rho_radec_mix = cov_radec_mix / (e_ra_mix * e_de_mix)

    mix_columns.append(fits.Column(name='ra_dec_corr', format='D', array=rho_radec_mix))

    print('Starting Photometric Extraction')
    # Add in photometric parameters. Choose here to pick them from vL
    for i in range(len(extra_data)):
        mix_columns.append(fits.Column(name=extra_data[i], format='D', array=hip2[extra_data[i]]))

    print('Starting HDUList Creation and Saving')
    table_hdu = fits.BinTableHDU.from_columns(mix_columns)

    table_hdu.writeto('../../data/Hipparcos_mix.fits', overwrite=True)
    print('Process Completed')


def extract_cov(tab):
    """Transform the u-values given in Hipparcos-2 (van Leeuwen, 2007) into a complete covariance matrix
       This code is based on a Python 2 implementation of the same procedure by A.G.A. Brown"""

    num = tab.shape[0]
    print(f'Converting U -> \Sigma for {num} Hip2 objects')

    # Setup arrays and parameters
    U = np.zeros((num, 5, 5))
    sigma = np.zeros((num, 5))
    nobs = tab['Ntr']
    gof = tab['f2']
    nu = nobs - 5 # There are some stars with n > 5 but this is a good approximation
    Q = nu * (np.sqrt(2. / (9. * nu)) * gof + 1 - (2. / (9. * nu))) ** 3.
    u = np.sqrt(Q / nu)  # confusing naming. Change?
    u_sq = u ** 2.

    # Populate sigma
    sigma_names = ['e_ra_rad', 'e_de_rad', 'e_plx', 'e_pm_ra', 'e_pm_de']
    for i, field in enumerate(sigma_names):
        sigma[:, i] = tab[field]

    # Populate U
    u_index = 1
    for j in range(5):
        for i in range(j + 1):
            U[:, i, j] = tab[f'u{u_index}']
            if i == j:
                U[:, i, j] *= (u / sigma[:, i])
            u_index += 1

    U_transpose = np.transpose(U, axes=(0, 2, 1))
    # In case of a 3D array, a matmul treats it as a stack of matrices and treats
    # the first index as the stack index.
    # N = Normal matrix (notation from Michalik et al. (2014))
    N = np.matmul(U_transpose, U)
    N_inverse = np.linalg.inv(N)
    covar = np.reshape(u_sq[:, None, None] * N_inverse, N_inverse.shape)

    return covar



def main():
    #hipp_hdul = fits.open('../../data/Hipparcos2.fits')
    #covar = extract_cov(hipp_hdul[1].data)
    mix_hipparcos()

if __name__ == '__main__':
    main()

############################
#
# Code and functions to perform cone-search based
# cross-matches between two object catalogues
#
############################

import numpy as np
import time
from numba import njit
from astropy.io import fits
import copy


@njit  # numba wrapper to improve processing speed by factor ~5
def conesearch_noerr(ra_s, de_s, ra_b, de_b,
                     dec_cutoff=3. / 3600.,
                     conesearch_radius=1. / 3600.):
    """Performs a simple error-less cone search around coordinates of the 'small' dataset to find objects in the 'big'
    dataset within a projected circle with conesearch_radius.

    Inputs: (ra, dec)_s [degrees]: coordinates of dataset to be crossmatched
            (ra, dec)_b [degrees]: coordinates of dataset to crossmatch to
            dec_cutoff  [degrees]: maximum difference between object in s and b to compute angular distance for
            conesearch_radius [degrees]: angular size of the cone

    Outputs: xm_table [ndarray of size k x 3]: array containing the indices of the cross-matched objects in s and b,
                                               and the angular size in degrees between them."""

    num_objects = ra_s.shape[0]
    xm_table = np.zeros((int(num_objects * 2), 3))
    k = 0
    for i in range(num_objects):
        dc_idxs = np.where(
            np.abs(de_b - de_s[i]) < dec_cutoff)  # Reduce computational load by calculating fewer distances?
        theta_ar = np.sqrt(
            ((ra_s[i] - ra_b[dc_idxs]) * np.cos(np.radians(de_b[dc_idxs]))) ** 2. + (de_s[i] - de_b[dc_idxs]) ** 2.)
        match = np.where(theta_ar < conesearch_radius)
        if len(match[0]) > 0:
            for j in match[0]:
                xm_table[k, :] = [i, dc_idxs[0][j], theta_ar[j] * 3600.]
                k += 1

    return xm_table[:k, :]


def read_gaia_hipp_data(gaia_path, hipp_path,
                        num_hipp='all',
                        **kwargs):
    """Read Gaia and Hipparcos data for cross-match purposes. """
    tab_g = fits.open(gaia_path)
    tab_h = fits.open(hipp_path)

    if not num_hipp == 'all':
        # NOTE: When cropping FITS Table like this, some header vals might be invalid
        tab_h[1].data = tab_h[1].data[:num_hipp]

    return tab_h, tab_g


def convert_to_ids(xm_table, tab_g, tab_h):
    """Takes an array with [hip_array_number, gaia_array_number, distance]
    and converts it to     [hip.hip, gaia.source_id, distance]"""
    tab = copy.deepcopy(xm_table)
    for i in range(xm_table.shape[0]):
        tab[i][0] = tab_h[1].data['hip'][int(xm_table[i][0])]
        tab[i][1] = tab_g[1].data['source_id'][int(xm_table[i][1])]

    return tab


def test_conesearch():
    # Mock data, replace with Gaia query later
    num_samples = 100

    # Real Data
    dpath = '../../data/'
    gaia_path = dpath + 'gaia_stars_sel12_noerr.fits'
    hipp_path = dpath + 'hipp_stars_noerr.fits'

    use_real_data = False

    # Conesearch args.
    dec_cutoff = 3. / 3600.  # deg, compute distance only for these objects
    conesearch_radius = 1. / 3600.  # deg
    ####

    if use_real_data:
        tab_h, tab_g = read_gaia_hipp_data(gaia_path, hipp_path, num_hipp=num_samples)
        ra_g = np.array(tab_g[1].data['ra_prop'], dtype=np.float64)
        de_g = np.array(tab_g[1].data['dec_prop'], dtype=np.float64)
        ra_h = np.array(tab_h[1].data['ra'], dtype=np.float64)
        de_h = np.array(tab_h[1].data['dec'], dtype=np.float64)
        print(f"Data Read.. ({ra_h.shape[0]})")


    print(f"""Performing Cross Match..\n- Simple {conesearch_radius * 3600} arcsecond cone search""")
    t0 = time.time()
    xm_table = conesearch_noerr(ra_h, de_h, ra_g, de_g)
    # xm_table = convert_to_ids(xm_table, tab_g, tab_h)
    t_e = time.time() - t0
    print(f"Conesearch completed, time elapsed: {t_e:.2f} s ({(t_e / ra_h.shape[0] * 1e3):2f} ms / object)")
    print(f"Number of cross match candidates:   {xm_table.shape[0]}")

    savename = f'conesearch_{int(conesearch_radius * 3600)}as_gaia_sel12_2'
    #np.save(f'../results/' + savename, xm_table)
    print(f"Table succesfully saved as {savename}.npy")


def main():
    test_conesearch()


if __name__ == '__main__':
    main()

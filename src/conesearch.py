import numpy as np
from astroquery.gaia import Gaia
import time
from numba import njit, jit, vectorize
from astropy.io import fits


@njit
def conesearch(ra_s, de_s, ra_b, de_b):
    """ Cross-Match in a loop-based way. Attempt to optimize using numba
    PROBLEM: Loop works relatively quickly, but data extraction has proven to be difficult"""

    num_objects = ra_s.shape[0]
    xm_table = np.zeros((int(num_objects * 1.2), 3))
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
        #else:
            #print(f'No 1" XM found for star {i}')  # with Hipparcos ID {tab_h[1].data["hip"][i]}')


    return xm_table[:k, :]


def create_mock_data(num_samples_h, num_samples_g=int(3e6)):
    """Make mock Hipparcos- and Gaia data.
    3e6 Gaia stars and num_samples Hipparcos objects"""
    # Hipparcos
    ra_h = np.random.rand(num_samples_h) * 360.  # deg [0, 360]
    de_h = np.random.rand(num_samples_h) * 180. - 90  # deg [-90, +90]

    ra_g = np.random.rand(num_samples_g) * 360.  # deg [0, 360]
    de_g = np.random.rand(num_samples_g) * 180. - 90.  # deg, [-90, +90]

    return ra_h, de_h, ra_g, de_g


def read_gaia_hipp_data(gaia_path, hipp_path, num_hipp='all'):
    tab_g = fits.open(gaia_path)
    tab_h = fits.open(hipp_path)

    if not num_hipp == 'all':
        # NOTE: When cropping FITS Table like this, some header vals might be invalid
        tab_h[1].data = tab_h[1].data[:num_hipp]

    return tab_g, tab_h


def convert_to_ids(xm_table, tab_g, tab_h):
    """Takes an array with [hip_array_number, gaia_array_number, distance]
    and converts it to     [hip.hip, gaia.source_id, distance]"""
    for i in range(xm_table.shape[1]):
        xm_table[i][0] = tab_h[1].data['hip'][int(xm_table[i][0])]
        xm_table[i][1] = tab_g[1].data['source_id'][int(xm_table[i][1])]

    return xm_table


# Mock data, replace with Gaia query later
num_samples = 'all'

# Real Data
dpath = '../data/'
gaia_path = dpath + 'gaia_stars_sel12_noerr.fits'
hipp_path = dpath + 'hipp_stars_noerr.fits'

use_real_data = True

# Conesearch args.
dec_cutoff = 3. / 3600.  # deg, compute distance only for these objects
conesearch_radius = 1. / 3600.  # deg
####

if use_real_data:
    tab_g, tab_h = read_gaia_hipp_data(gaia_path, hipp_path, num_hipp=num_samples)
    ra_g = np.array(tab_g[1].data['ra_prop'], dtype=np.float64)
    de_g = np.array(tab_g[1].data['dec_prop'], dtype=np.float64)
    ra_h = np.array(tab_h[1].data['ra'], dtype=np.float64)
    de_h = np.array(tab_h[1].data['dec'], dtype=np.float64)
    print(f"Data Read.. ({ra_h.shape[0]})")
else:
    ra_h, de_h, ra_g, de_g = create_mock_data(num_samples)
    print(f"Mock Data Created.. ({ra_h.shape[0]} objects)")

print(f"""Performing Cross Match..\n- Simple {conesearch_radius * 3600} arcsecond cone search""")
t0 = time.time()
xm_table = conesearch(ra_h, de_h, ra_g, de_g)
#xm_table = convert_to_ids(xm_table, tab_g, tab_h)
t_e = time.time() - t0
print(f"Conesearch completed, time elapsed: {t_e:.2f} s ({(t_e/ra_h.shape[0])*1e3:2f} ms / object)")
print(f"Number of cross match candidates:   {xm_table.shape[0]}")

savename = f'conesearch_{int(conesearch_radius * 3600)}as_gaia_sel12'
np.save(f'../results/'+savename, xm_table)
print(f"Table succesfully saved as {savename}")

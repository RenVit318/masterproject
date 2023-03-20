############################
#
# Functions for the second step in the crossmatch pipeline
# Multiple methods of describing the best neighbour for a
# Hipparcos star with multiple possible matches
#
############################
from typing import List

import numpy as np
from astropy.io import fits
from data_queries import query_extra_data
from astrometry_equations import predict_G_bv, compute_error_normalized_distance


def nearest(matches, distances, extra_gaia_data, extra_hipp_data):
    """Determines the best match based only on distance
    Requires the input to be: (ID1, ID2, distance)"""
    best_match_idx = np.argmin(distances)
    return matches[best_match_idx]


def brightest(matches, distances, extra_gaia_data, extra_hipp_data):
    """Determines the best match based only on the brightest object"""
    mag = extra_gaia_data
    best_match_idx = np.argmin(mag)
    return matches[best_match_idx]


def likeliest_magnitude(matches, distances, extra_gaia_data, extra_hipp_data):
    """Predict the Gaia magnitude using Hp_mag and B-V colour and return the smallest difference between G and G_pred"""
    G = extra_gaia_data
    hp_mag = extra_hipp_data[:, 0]
    b_v = extra_hipp_data[:, 1]

    G_pred, _ = predict_G_bv(hp_mag, b_v)
    best_match_idx = np.argmin(np.abs(G_pred - G))
    return matches[best_match_idx]


def likeliest_position(matches, distances, extra_gaia_data, extra_hipp_data):
    """Determines the best match based on distance and errors using the error normalized computation code"""
    pos_gaia = np.array([extra_gaia_data[:, 0], extra_gaia_data[:, 1]]) * 3.6e6 # Convert deg -> mas
    pos_hipp = np.array([extra_hipp_data[:, 0], extra_hipp_data[:, 1]]) * 3.6e6
    unc_gaia = np.array([extra_gaia_data[:, 2]**2, extra_gaia_data[:, 3]**2, extra_gaia_data[:, 4]])
    unc_hipp = np.array([extra_hipp_data[:, 2]**2, extra_hipp_data[:, 3]**2, extra_hipp_data[:, 4]])

    D = compute_error_normalized_distance(pos_gaia, pos_hipp, unc_gaia, unc_hipp, method='full')
    best_match_idx = np.argmin(D)
    return matches[best_match_idx]


def check_hip_sorted(array):
    """If called, performs a check whether or not the Hipparcos array is sorted low to high
    This is because the best neighbour selection function works on the assumption that they are"""
    for i in range(len(array) - 1):
        if array[i + 1] < array[i]:
            print("This array is not sorted")
            return False
    print("This array is sorted")
    return True


def merge_extra_data(tab_xm, extra_data_names, cat, cat_type=None, dpath='../../data/'):
    """Get data from one of the GaiaCats"""
    full_cat = fits.open(dpath+cat)[1].data
    sel_cat = np.zeros((tab_xm.shape[0], len(extra_data_names)))

    if cat_type == 'gaia':
        indexes = tab_xm[:, 2]
    elif cat_type == 'hipp':
        indexes = np.zeros(tab_xm.shape[0], dtype=np.int64)
        for i in range(tab_xm.shape[0]):
            indexes[i] = int(np.where(full_cat['hip'] == tab_xm[i][0])[0][0])

    for i, name in enumerate(extra_data_names):
        sel_cat[:, i] = full_cat[name][indexes[i]]

    return tab_xm, sel_cat


def select_best_neighbour(tab_xm, distances, best_match_selection, GaiaCat=None, HippCat=None, **kwargs):
    """From a table containing all possible cross-matches, selects the best neighbour
       for each hipparcos object with >1 reported neighbour"""

    # Check which selection type is inputted, and appended what we need to do to get extra data
    # Options are:
    #   - None (can be anything) : No extra data is required for the selection
    #   - merge : Merge xm-table with a GaiaCat in Python
    #   - query : Merge xm-table with the Gaia table in the archive.
    #             Slowest by far. Use only if required data not in the GaiaCat

    if best_match_selection == 'nearest':
        select_func = nearest
        get_extra_gaia_data_type = None
        get_extra_hipp_data_type = None

    elif best_match_selection == 'brightest':
        select_func = brightest
        get_extra_gaia_data_type = 'query'  # CHANGE. For testing purposes
        extra_gaia_data_names = ['phot_g_mean_mag']

    elif best_match_selection == 'likeliest_magnitude':
        select_func = likeliest_magnitude
        # DON'T CHANGE INDEXING HERE AS IT WILL MESS UP THE SELECTION CODE
        get_extra_gaia_data_type = 'query'
        extra_gaia_data_names = ['phot_g_mean_mag']
        get_extra_hipp_data_type = 'query'
        extra_hipp_data_names = ['hp_mag', 'b_v']

    elif best_match_selection == 'likeliest_position':
        select_func = likeliest_position
        # DON'T CHANGE INDEXING HERE AS IT WILL MESS UP THE SELECTION CODE
        get_extra_gaia_data_type = 'merge'
        extra_gaia_data_names = ['ra_prop', 'dec_prop', 'e_ra_prop', 'e_de_prop', 'ra_dec_prop']  # Also need Hipparcos errors here?
        get_extra_hipp_data_type = 'merge'
        extra_hipp_data_names = ['ra', 'dec', 'e_ra_rad', 'e_de_rad', 'ra_dec_corr']

    # We cannot easily sort the array due to RAM constraints. So we have to assume Hip array is sorted
    if not check_hip_sorted(tab_xm[:, 0]):
        raise NotImplementedError(
            "The provided Hipparcos ID's are not sorted. Currently this function assumes they are. Fix this")

    # Gather ancillary data from the Gaia Archive. 
    if get_extra_gaia_data_type == 'merge':
        _, extra_gaia_data_table = merge_extra_data(tab_xm, extra_gaia_data_names, cat=GaiaCat, cat_type='gaia')
    elif get_extra_gaia_data_type == 'query':
        _, extra_gaia_data_table = query_extra_data(tab_xm, extra_gaia_data_names)
    else:
        print(f"Extra data parameter {get_extra_gaia_data_type} unknown. Continuing without gathering extra data")

    # Gather ancillary data from the Hipparcos catalogues
    if get_extra_hipp_data_type == 'merge':
        _, extra_hipp_data_table = merge_extra_data(tab_xm, extra_hipp_data_names, cat=HippCat, cat_type='hipp')
    elif get_extra_hipp_data_type == 'query':
        _, extra_hipp_data_table = query_extra_data(tab_xm, extra_hipp_data_names, cat='hipp')
    else:
        print(f"Extra data parameter {get_extra_hipp_data_type} unknown. Continuing without gathering extra data")

    vals, idx_start, count = np.unique(tab_xm[:, 0], return_counts=True, return_index=True)

    best_matches_array = np.zeros((len(vals), 3), dtype=np.int64)  # save the best matches in here

    # Call selection function
    for i in range(len(vals)):
        if count[i] > 1:
            matches = tab_xm[idx_start[i]:idx_start[i] + count[i]]
            extra_gaia_data = extra_gaia_data_table[idx_start[i]:idx_start[i] + count[i]]
            extra_hipp_data = extra_hipp_data_table[idx_start[i]:idx_start[i] + count[i]]
            best_match = select_func(matches, distances, extra_gaia_data, extra_hipp_data)

        else:
            best_match = tab_xm[idx_start[i]]
        best_matches_array[i] = best_match

    return best_matches_array


def main():
    r_path = '../results/'
    xm_tab_name = 'crossmatch+PMC+EIB_1as_all_neighbours.npy'
    best_match_selection = 'brightest'

    tab_xm = np.load(r_path + xm_tab_name)
    print(tab_xm)
    addition_index = xm_tab_name.find('_all')
    savename = xm_tab_name[0:addition_index]

    print(f'Starting single selection run w/ {best_match_selection}')
    print(f"Saving into '{savename}'...")
    best_matches_array = select_best_neighbour(tab_xm, None, best_match_selection)
    print(best_matches_array)
    #_, res = query_extra_data(best_matches_array, ['phot_g_mean_mag'])


if __name__ == '__main__':
    main()

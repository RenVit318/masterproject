############################
#
# Functions for the second step in the crossmatch pipeline
# Multiple methods of describing the best neighbour for a
# Hipparcos star with multiple possible matches
#
############################
import numpy as np
from data_queries import query_extra_data


def nearest(matches, distances):
    """Determines the best match based only on distance
    Requires the input to be: (ID1, ID2, distance)"""
    best_match_idx = np.argmin(distances)
    return matches[best_match_idx]


def brightest(matches, mag):
    """Determines the best match based only on the brightest object"""
    best_match_idx = np.argmin(mag)
    return matches[best_match_idx]


def likeliest_position():
    """Determines the best match based on distance and errors"""
    pass


def check_hip_sorted(array):
    """If called, performs a check whether or not the Hipparcos array is sorted low to high
    This is because the best neighbour selection function works on the assumption that they are"""
    for i in range(len(array) - 1):
        if array[i + 1] < array[i]:
            print("This array is not sorted")
            return False
    print("This array is sorted")
    return True


def merge_extra_data():
    """Get data from one of the GaiaCats. Can do this with pandas"""
    tab_xm, extra_data_table = None, None
    return tab_xm, extra_data_table


def select_best_neighbour(tab_xm, distances, best_match_selection, GaiaCat=None):
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
        get_extra_data = None  # Defines if we need to query for extra data
        extra_data_table = distances
    elif best_match_selection == 'brightest':
        select_func = brightest
        get_extra_data = 'query'  # CHANGE. For testing purposes
        extra_data_names = ['phot_g_mean_mag']
    elif best_match_selection == 'likeliest_position':
        select_func = likeliest_position
        get_extra_data = 'merge'
        extra_data_names = ['e_ra_prop, e_de_prop, ra_dec_prop']  # Also need Hipparcos errors here?

    # We cannot easily sort the array due to RAM constraints. So we have to assume Hip array is sorted
    if not check_hip_sorted(tab_xm[:, 0]):
        raise NotImplementedError(
            "The provided Hipparcos ID's are not sorted. Currently this function assumes they are. Fix this")

    # Gather ancillary data. 
    if get_extra_data == 'merge':
        tab_xm, extra_data_table = merge_extra_data(tab_xm, extra_data_names)
    elif get_extra_data == 'query':
        tab_xm, extra_data_table = query_extra_data(tab_xm, extra_data_names)
    else:
        print(f"Extra data parameter {get_extra_data} unknown. Continuing without gathering extra data")

    vals, idx_start, count = np.unique(tab_xm[:, 0], return_counts=True, return_index=True)

    best_matches_array = np.zeros((len(vals), 3), dtype=np.int64)  # save the best matches in here

    # Call selection function
    for i in range(len(vals)):
        if count[i] > 1:
            matches = tab_xm[idx_start[i]:idx_start[i] + count[i]]
            extra_data = extra_data_table[idx_start[i]:idx_start[i] + count[i]]
            best_match = select_func(matches, extra_data)

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

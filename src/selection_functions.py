############################
#
# Functions for the second step in the crossmatch pipeline
# Multiple methods of describing the best neighbour for a
# Hipparcos star with multiple possible matches
#
############################
import numpy as np
from data_queries import query_extra_data

def nearest(matches):
    """Determines the best match based only on distance
    Requires the input to be: (ID1, ID2, distance)"""
    best_match_idx = np.argmin(matches[:,2])
    return matches[best_match_idx, 0:2]

def brightest(matches, mag):
    """Determines the best match based only on the brightest object"""
    best_match_idx = np.argmin(mag)
    return matches[best_match_idx, 0:2]


def likeliest_position():
    """Determines the best match based on distance and errors"""
    pass


def check_hip_sorted(array):
    """If called, performs a check whether or not the Hipparcos array is sorted low to high
    This is because the best neighbour selection function works on the assumption that they are"""
    for i in range(len(array)-1):
        if array[i+1] < array[i]:
            print("This array is not sorted")
            return False
    print("This array is sorted")
    return True
            

def merge_extra_data(tab_xm, extra_data_names, GaiaCat):
    """Merge tab_xm and selected attributes from GaiaCat together based on source_id"""
    


def select_best_neighbour(tab_xm, best_match_selection, savename, GaiaCat=None):
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
        get_extra_data = None # Defines if we need to query for extra data
    elif best_match_selection == 'brightest':
        select_func = brightest
        get_extra_data = 'query' # CHANGE. For testing purposes
        extra_data_names = ['phot_g_mean_mag']
    elif best_match_selection == 'likeliest_position':
        select_func = likeliest_position
        get_extra_data = 'merge'
        extra_data_names = ['e_ra_prop, e_de_prop, ra_dec_prop']  # Also need Hipparcos errors here?

    # Gather ancillary data. 
    if get_extra_data == 'merge':
        tab_xm, extra_data_table = merge_extra_data(tab_xm, extra_data_names, GaiaCat)
    elif get_extra_data == 'query':
        tab_xm, extra_data_table = query_extra_data(tab_xm[:100], extra_data_names)
    else:
        print(f"Extra data parameter {get_extra_data} unkown. Continuing without querying gathering data")
    

    # Check all objects with >1 candidate. Based on
    # https://stackoverflow.com/questions/30003068/how-to-get-a-list-of-all-indices-of-repeated-elements-in-a-numpy-array
    # The sorting step is probably not necessary because the crossmatch algorithm automatically sorts on hip. But we do do it
    #hip = tab_xm[:, 0]
    #hip_sort_idxs = np.argsort(hip)  
    #tab_xm_sorted = tab_xm[:, hip_sort_idxs]
    # But, we cannot easily sort the array due to RAM constraints. So we have to assume Hip array is sorted
    if not check_hip_sorted(tab_xm[:,0]):
        raise NotImplementedError("The provided Hipparcos ID's are not sorted. Currently this function assumes they are. Fix this")

    vals, idx_start, count = np.unique(tab_xm[:,0], return_counts=True, return_index=True)
    
    best_matches_array = np.zeros((len(vals), 2)) # save the best matches in here
    
    # Call selection function
    for i in range(len(vals)):
        if count[i] > 1:
            matches = tab_xm[idx_start[i]:idx_start[i]+count[i]]
            if get_extra_data is not None:
                extra_data = extra_data_table[idx_start[i]:idx_start[i]+count[i]]
                best_match = select_func(matches, extra_data)
            else:
                best_match = select_func(matches)
        else:
            best_match = tab_xm[idx_start[i], 0:2] # get the Hipparcos and Gaia ID
        best_matches_array[i] = best_match       

    return best_matches_array



def main():
    r_path = '../results/'
    xm_tab_name = 'crossmatch+PMC+EIB_1as_all_neighbours.npy'
    best_match_selection = 'brightest'

    tab_xm = np.load(r_path + xm_tab_name)
    addition_index = xm_tab_name.find('_all')
    savename = xm_tab_name[0:addition_index]

    print(f'Starting single selection run w/ {best_match_selection}')
    print(f"Saving into '{savename}'...")
    best_matches_array = select_best_neighbour(tab_xm, best_match_selection, savename)




if __name__ == '__main__':
    main()

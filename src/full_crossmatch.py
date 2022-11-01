############################
#
# Complete cross-match function capable of running all cross-match related code
# Currently can read or query data (+ all necessary preprocessing)
# and perform an n arcsecond conesearch selecting all neighbours
#
# PLANNED: Best neighbour selection functions, automatic analytics creation
#
############################

import time
from conesearch import conesearch_noerr, save_xm_results
from data_queries import read_gaia_hipp_data, gaia_login
from table_functions import extract_sky_positions, convert_to_ids
from gaia_preprocess import full_preprocess
from selection_functions import select_best_neighbour

import numpy as np #Remove

def full_run_crossmatch(data_type, data_kwargs, conesearch_params, 
                        save_file=False, savename='crossmatch'):
    """Execute a complete iteration of the cross-match algorithm with either a locally stored table or a queried table from
       the Gaia Archive. The latter is not really feasible for mag_lim>~12 as processing takes up to 8 hours, store those locally
       Procedure:
        - Import data locally or with full_preprocess() from gaia_preprocess.py"""

    t0 = time.time()
    # Gather cross-match data
    if data_type == 'local':
        tab_h, tab_g = read_gaia_hipp_data(**data_kwargs)
    elif data_type == 'query_gaia':
        gaia_login()
        tab_g = full_preprocess(return_cat=True, save_cat=False, **data_kwargs)
        tab_h, _ = read_gaia_hipp_data(**data_kwargs)

    ra_h, de_h = extract_sky_positions(tab_h)
    ra_g, de_g = extract_sky_positions(tab_g, ra_id='ra_prop', de_id='dec_prop')
    

    print(f"Data Imported. Import runtime: {time.time() - t0:.2f} s\nNo. Gaia objects: {ra_g.shape[0]}\nNo. Hipparcos objects: {ra_h.shape[0]}")

    print(f"Starting {int(conesearch_params['conesearch_radius']*3600)} arcsecond Conesearch")
    t1 = time.time()
    tab_xm = conesearch_noerr(ra_h, de_h, ra_g, de_g, conesearch_params['conesearch_radius'])
    unique_candidates, c = np.unique(tab_xm, axis=0,return_counts=True) 

    tab_xm_ids = convert_to_ids(tab_xm, tab_g, tab_h)
    u,c = np.unique(tab_xm_ids, axis=0, return_counts=True)

    print(f"Conesearch Done. Conesearch runtime: {time.time() - t1} s")

    # Save functions
    if save_file:
        save_xm_results(tab_xm_ids, savename+'_all_neighbours')
        print(f"Cross-Match results saved as {savename}_all_neighbours.npy")

    selection_type = conesearch_params['best_match_selection']
    if selection_type is not None:
        one_xm_tab = select_best_neighbour(tab_xm_ids, **conesearch_params)
        
        if save_file:
            save_xm_results(one_xm_tab, f'{savename}_best_neighbour_{selection_type}')
            print(f"Cross-Match results saved as {savename}_best_neighbour.npy")
    

    print(f"One Cross-Match Completed.\nTotal runtime: {time.time() - t0} s")

def make_crossmatch_savename(cat, conesearch_radius, basename='crossmatch'):
    """"""
    cat_plus_idx = cat.find('+')
    conesearch_as = int(conesearch_radius * 3600) # 'breaks' for non-integer (in as) radii 
    extensions = '' 
    if cat_plus_idx != -1: # If it is, there is no '+' in the string
        extensions = cat[cat_plus_idx:cat.find('.')] # cat name ends in .fits

    savename = basename + extensions + f'_{conesearch_as}as'
    return savename
    

def experiment_conesearch():
    """"""
    ######### USER PARAMETERS #########
    # Data
    dpath = '../../data/'
    data_type = 'local'  # local or query_gaia
    gaia_cat = 'GaiaBaseCat+PMC+EIB.fits'
    data_kwargs = {
        'gaia_path': dpath + gaia_cat,                        # Data Path variables
        'hipp_path': dpath + 'hipp_stars_noerr.fits',
        'min_g_mag': 14,                                      # Preprocess variables
        'apply_pm_corr': False,                               # Usually aren't going to need these, but just in case
        'error_inflation_type': 'Brandt21',
        'batch_size': int(1e5),
        'num_hipp': 'all'
    }

    # Conesearch
    conesearch_params = {
        'conesearch_radius': 1./3600.,
        'best_match_selection': None
    }
    ###################################

    savename = make_crossmatch_savename(gaia_cat, conesearch_params['conesearch_radius'])
    print(f"Starting full crossmatch run. Saving into basename: {savename}")
    full_run_crossmatch(data_type, data_kwargs, conesearch_params, save_file=True,
                        savename=savename)

    ## Analytical functions here ##


def main():
    experiment_conesearch()


if __name__ == '__main__':
    main()

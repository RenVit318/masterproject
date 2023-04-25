############################
#
# All of the functions to be called when we want to run a
# large amount of one code in one time/loop over
# functions with different settings
#
# Here we can easily work in the functions without breaking anything
#
############################

import numpy as np
from data_queries import delete_unlabeled_jobs
from gaia_preprocess import full_preprocess
from hipp_preprocess import mix_hipparcos
from full_crossmatch import full_run_crossmatch, make_crossmatch_savename
from selection_functions import select_best_neighbour
from conesearch import save_xm_results

def make_all_catalogues(mag_lim=14, gaia_epoch=2016.0, hipp_epoch=1991.25, batch_size=int(1e5),
                        read_local=False, data_path=None): # maybe save this file somewhere if ever necessary
    """Given a certain set of standard values, iterates over all tunable variables"""

    for apply_pm_corr in [True, False]:
        for error_inflation_type in [None, 'Brandt21']:
            # Skip Selection
            
            # Make name
            savename = "GaiaBaseCat"
            if apply_pm_corr:
                savename += "+PMC"  # Proper Motion Correction
            if error_inflation_type == 'Brandt21':
                savename += "+EIB"  # Error Inflation Brandt

            print(f"Starting Gaia Preprocessing Into: {savename}")

            full_preprocess(mag_lim, gaia_epoch, hipp_epoch, batch_size=batch_size, read_local=read_local, data_path=data_path,
                            apply_pm_corr=apply_pm_corr, error_inflation_type=error_inflation_type,
                            savename=savename)
    delete_unlabeled_jobs(no_check=True)


def make_xm_catalogues():

    dpath = '../../data/'
    data_type = 'local'  # local or query_gaia
    catalogues = ['GaiaBaseCat.fits', 'GaiaBaseCat+PMC.fits', 'GaiaBaseCat+EIB.fits', 'GaiaBaseCat+PMC+EIB.fits']
    conesearch_radii_as = [1, 2, 3, 5, 10]

    for gaia_cat in catalogues:
        for radius in conesearch_radii_as:

            # Crossmatch Parameters
            data_kwargs = {
                'gaia_path': dpath + gaia_cat,  # Data Path variables
                'hipp_path': dpath + 'hipp_stars_noerr.fits',
            }
            conesearch_params = {
                'conesearch_radius': radius / 3600.,
                'best_match_selection': 'all'
            }

            savename = make_crossmatch_savename(gaia_cat, conesearch_params['conesearch_radius'])
            print(f"Starting Crossmatch Iteration: {savename}")
            full_run_crossmatch(data_type, data_kwargs, conesearch_params,
                                save_file=True, savename=savename)


def make_marrese_xm_catalogues():
    """Crossmatch with only the Marrese objects to study statistics on 'good' stars"""
    dpath = '../../data/'
    data_type = 'local'
    catalogues = ['GaiaBaseCat+EIB.fits', 'GaiaBaseCat+PMC+EIB.fits']
    conesearch_radii_as = [1, 2, 5]

    for gaia_cat in catalogues:
        for radius in conesearch_radii_as:
            # Crossmatch Parameters
            data_kwargs = {
                'gaia_path': dpath + gaia_cat,  # Data Path variables
                'hipp_path': dpath + 'marrese_hipp_table.fits',
            }
            conesearch_params = {
                'conesearch_radius': radius / 3600.,
                'best_match_selection': 'all'
            }

            savename = make_crossmatch_savename(gaia_cat, conesearch_params['conesearch_radius'], basename='marrese_crossmatch')
            print(f"Starting Crossmatch Iteration: {savename}")
            full_run_crossmatch(data_type, data_kwargs, conesearch_params,
                                save_file=True, savename=savename)


def full_crossmatch_result():
    """Crossmatch with all Hipparcos objects to make our own xm table"""
    dpath = '../../data/'
    data_type = 'local'  # local or query_gaia
    catalogues = ['GaiaBaseCat+PMC+EIB.fits']
    hipp_cat = 'Hipparcos_mix.fits'
    conesearch_radii_as = [2, 3, 5, 10]

    for gaia_cat in catalogues:
        for radius in conesearch_radii_as:

            # Crossmatch Parameters
            data_kwargs = {
                'gaia_path': dpath + gaia_cat,  # Data Path variables
                'hipp_path': dpath + hipp_cat,
                'GaiaCat': gaia_cat,
                'HippCat': hipp_cat
            }
            conesearch_params = {
                'conesearch_radius': radius / 3600.,
                'best_match_selection': 'all'
            }

            savename = make_crossmatch_savename(gaia_cat, conesearch_params['conesearch_radius'], basename='final_crossmatch')
            print(f"Starting Crossmatch Iteration: {savename}")
            full_run_crossmatch(data_type, data_kwargs, conesearch_params,
                                save_file=True, savename=savename)    

def select_best_matches():
    """Use already existing 'all-neighbour' crossmatches to select the best matches with
    new selection functions."""
    rpath = '../results/'
    xm_basename = 'final_crossmatch+PMC+EIB_'
    GaiaCat = 'GaiaBaseCat+PMC+EIB.fits'
    HippCat = 'Hipparcos_mix.fits'
    conesearch_radii_as = [2, 3, 5, 10]
    selection_methods = ['likeliest_magnitude', 'likeliest_position']
    
    for cs_radius in conesearch_radii_as:
        cs_basename = f'{xm_basename}{cs_radius}as_'
        tab_xm = np.load(f'{rpath}{cs_basename}all_neighbours.npy')
        
        for method in selection_methods:
            best_matches = select_best_neighbour(tab_xm, None, method, GaiaCat=GaiaCat, HippCat=HippCat)
            savename = f'{cs_basename}_best_neighbour_{method}.npy'
            save_xm_results(best_matches, savename)
    

def make_final_cats():
    """Make the final Gaia Base Cat including everything from Brandt21 and Cantat-Gaudin+21,
    and the final Hipparcos mix catalogue, where we include all known uncertainties."""
    # 1. GaiaCat
    full_preprocess(mag_lim=14, gaia_epoch=2016., hipp_epoch=1991.25, batch_size=int(1e5), read_local=True,
                    data_path='../../data/gaia_maglim14.fits', apply_pm_corr=True, error_inflation_type='Brandt21', apply_color_correction=True,
                    apply_zp_corr=True, return_cat=False, save_cat=True, savename='GaiaBaseCat_all')

    #2. Hipparcos_mix
#    mix_hipparcos() # oops sligthly less well documented :(
    


def main():
    #make_marrese_xm_catalogues()
    #full_crossmatch_result()
    #select_best_matches()
    make_final_cats()


if __name__ == '__main__':
    main()


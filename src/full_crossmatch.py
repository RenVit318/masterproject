import time
from astroquery.gaia import Gaia
from conesearch import conesearch_noerr, save_xm_results
from data_queries import query_gaia_simple_conesearch, read_gaia_hipp_data, gaia_login
from table_functions import extract_sky_positions, convert_to_ids


def full_run_crossmatch(data_type, data_kwargs, conesearch_params, save_file=False):
    """"""

    t0 = time.time()
    # Gather cross-match data
    if data_type == 'local':
        tab_h, tab_g = read_gaia_hipp_data(**data_kwargs)
    elif data_type == 'query_gaia':
        tab_g = query_gaia_simple_conesearch(**data_kwargs)
        tab_h, _ = read_gaia_hipp_data(**data_kwargs)

    ra_h, de_h = extract_sky_positions(tab_h)
    ra_g, de_g = extract_sky_positions(tab_g)

    print(f"Data Imported. Import runtime: {time.time() - t0:.2f} s\nNo. Gaia objects: {ra_g.shape[0]}\nNo. Hipparcos objects: {ra_h.shape[0]}")

    print(f"Starting {int(conesearch_params['conesearch_radius']*3600)} arcsecond Conesearch")
    t1 = time.time()
    tab_xm = conesearch_noerr(ra_h, de_h, ra_g, de_g, **conesearch_params)
    tab_xm_ids = convert_to_ids(tab_xm, tab_g, tab_h)
    print(f"Conesearch Done. Conesearch runtime: {time.time() - t1} s")

    # Save functions
    savename = 'conesearch_1as_13p7'
    save_xm_results(tab_xm_ids, savename)
    print(f"Cross-Match results saved as {savename}")

    print(f"One Cross-Match Completed.\nTotal runtime: {time.time() - t0} s")


def experiment_conesearch():
    """"""
    ######### USER PARAMETERS #########
    # Data
    dpath = '../../data/'
    data_type = 'query_gaia'  # local or query_gaia
    data_kwargs = {
        'gaia_path': dpath+'gaia_stars_sel12_noerr.fits',   # Data Path variables
        'hipp_path': dpath + 'hipp_stars_noerr.fits',
        'min_p_mag': 13.7,                                    # Query variables
        'back_prop_gaia': True,
        'apply_frame_rot': False

    }

    # Conesearch
    conesearch_params = {
        'conesearch_radius': 1./3600.
    }
    ###################################
    gaia_login()
    full_run_crossmatch(data_type, data_kwargs, conesearch_params, save_file=True)

    ## Analytical functions here ##


def main():
    experiment_conesearch()


if __name__ == '__main__':
    main()

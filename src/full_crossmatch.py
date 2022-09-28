from astroquery.gaia import Gaia
from conesearch import conesearch_noerr, read_gaia_hipp_data
from data_queries import query_simple_conesearch

def full_run_crossmatch():






def experiment_conesearch():
    """"""
    ######### USER PARAMETERS #########
    # Data
    dpath = '../../data/'
    data_type = 'local'  # local or query
    data_kwargs = {
        'gaia_path': dpath+'gaia_stars_sel12_noerr.fits',   # Data Path variables
        'hipp_path': dpath + 'hipp_stars_noerr.fits',
        'min_p_mag': 12,                                    # Query variables
        'apply_frame_rot': False
    }

    # Conesearch
    conesearch_params = {
        'dec_cutoff': 3./3600.,
        'conesearch_radius': 1./3600.
    }
    ###################################

    # Gather cross-match data
    if data_type == 'local':
        tab_h, tab_g = read_gaia_hipp_data(**data_kwargs)
    elif data_type == 'query':
        tab_h, tab_g = query_simple_conesearch(**data_kwargs)

    # Question on implementation: How do we extract ra and dec information from queries tables
    # ideally we would do this in a way that is consistent for both .fits files and for queried tables

    tab_xm = conesearch_noerr(ra_h, de_h, ra_g, de_g, **conesearch_params)


def main():
    experiment_conesearch()


if __name__ == '__main__':
    main()

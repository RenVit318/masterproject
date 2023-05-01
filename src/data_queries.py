############################
#
# !!NOTE: THIS FILE CONTAINS AUTOMATED GAIA LOGIN, PLEASE CHANGE USER WHEN RUNNING THIS CODE!!
#
# All functions requiring communication with the Gaia Archive.
# Data Queries, Job Removal, etc..
#
############################

from astroquery.gaia import Gaia
from astropy.io import fits
from astropy.table import Table
import time


def gaia_login(user='rkievit', password="Gaia3-Hipp2"):
    global username
    username = user
    Gaia.login(user=user, password=password)


def launch_job(query, output_format='votable'):
    job = Gaia.launch_job_async(query=query , output_format=output_format)
    return job


def get_data(job):
    return job.get_results()


def delete_unlabeled_jobs(no_check=False):
    """ Deletes all jobs in user directory!!!
    BELOW IS CURRENTLY UNTRUE. HOW TO FIX THIS??
    Delete all jobs that weren't given a specific name. Takes on two assumptions:
    1. The job name string is exactly 14 characters long
    2. The job name ends with the letter 'O'
    Note that: therefore, if any given job name has these characters, it is deleted."""

    if no_check:  # Extra safety point
        pass
    else:
        delete_check = input(
            "Running this code means all Gaia jobs of logged in user will be deleted, continue?.. (y/n)")
        if delete_check == 'y':
            pass
        else:
            print("Understood. Not deleting jobs")
            return

    gaia_login()  # Just make sure we're logged in
    jobs = Gaia.list_async_jobs()
    jobs_to_remove = []

    for job in jobs:
        print(job)
        if job.jobid.endswith('O') and len(job.jobid) == 14:  # True for all jobs
            jobs_to_remove.append(job.jobid)
            # Gaia.remove_jobs(job.jobid)
    #    else:
    #        print(f"Job not removed: {job.jobid}")
    Gaia.remove_jobs(jobs_to_remove)


def read_gaia_hipp_data(gaia_path, hipp_path,
                        num_hipp='all',
                        **kwargs):
    """Read Gaia and Hipparcos data for cross-match purposes. """
    tab_g = fits.open(gaia_path)
    tab_h = fits.open(hipp_path)

    if not num_hipp == 'all':
        # NOTE: When cropping FITS Table like this, some header vals might be invalid
        tab_h[1].data = tab_h[1].data[:num_hipp]

    return tab_h[1].data, tab_g[1].data


def query_gaia_simple_conesearch(mag_lim=None,
                                 back_prop_gaia=False,
                                 epoch_1=1991.25, epoch_2=2016.0,
                                 apply_frame_rot=False,
                                 **kwargs):
    """Queries Gaia IDs, RA and Dec (back-propagated) for a simple cone search cross-match
    NOTE: propagation takes a long time ~90m for 3e6 sources, after preprocess can do everything locally."""
    query = ""
    query = query + f"SELECT source_id, ra, dec, parallax, pmra, pmra_error, pmdec, pmdec_error, phot_g_mean_mag "
    if back_prop_gaia:
        query = query + ",array_element(a0, 1) as ra_prop, \
                          array_element(a0, 2) as dec_prop, \
                          array_element(a0, 3) as parallax_prop, \
                          array_element(a0, 4) as pmra_prop, \
                          array_element(a0, 5) as pmdec_prop, \
                          array_element(a0, 6) as rv_prop "

        query = query + f"FROM \
                         ( \
                         SELECT *, EPOCH_PROP(ra, dec, parallax, pmra, pmdec, radial_velocity, {epoch_2}, {epoch_1}) as a0  "

    query = query + "FROM gaiadr3.gaia_source AS gaia "

    if mag_lim is not None:
        query = query + f"WHERE phot_g_mean_mag < {mag_lim} "

    if back_prop_gaia:
        query = query + ') as p'

    job = launch_job(query)
    return get_data(job)


def query_gaia_preprocess(mag_lim):
    """"Query all 6-parameter astrometry with the associated covariance components of Gaia data with M < mag_lim """
    query = f""" SELECT source_id, ra, ra_error, dec, dec_error, parallax, parallax_error, pmra, pmra_error, pmdec, 
                       pmdec_error, radial_velocity, radial_velocity_error, ra_dec_corr, ra_parallax_corr, ra_pmra_corr,
                       ra_pmdec_corr, dec_parallax_corr, dec_pmra_corr, dec_pmdec_corr,	parallax_pmra_corr, 
                       parallax_pmdec_corr, pmra_pmdec_corr, phot_g_mean_mag, bp_rp,
                       nu_eff_used_in_astrometry, pseudocolour, ecl_lat, astrometric_params_solved
			    FROM gaiadr3.gaia_source 
			    WHERE phot_g_mean_mag < {mag_lim} """
    job = launch_job(query)
    return get_data(job)


def query_gaia_zpcorr(gaia_ids, table_name='source_id_zp'):
    """Query the 5 parameters required to compute the zero point correction from Lindegren et al. (2021)"""
    import numpy as np
    gaia_login()

    table = Table([gaia_ids], names=['source_id']) 
    Gaia.upload_table(upload_resource=table, table_name=table_name)
    
    try:
        query = f"""SELECT source_id, parallax, phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolour, ecl_lat, astrometric_params_solved
                    FROM gaiadr3.gaia_source as gaia
                    JOIN user_{username}.{table_name} USING (source_id)"""
        job = launch_job(query, output_format='fits')
        data = get_data(job)
        Gaia.delete_user_table(table_name)            
    
        # Check if we are not shuffled
        if not np.array_equal(data['source_id'], gaia_ids):
            raise ValueError("Crossmatches array was shuffled during transfer with the Archive.")
    
        return data
    except: # remove user table if smth breaks
        Gaia.delete_user_table(table_name)
        raise

def propagate_batches_error(batches, epoch1, epoch2, num_runs=2):
    """Multithreading solution to propagating Gaia positions and errors. Full runtime of 16e6 entries takes >2h
    But by batching we reduce batch runtime to only ~few minutes
    TODO: increase archive communiation time? If slow remove meta batches because we only have a few jobs running
    at the same time, so we can just remove the user table in the individual function"""
    import threading
    # How many threads are passively active ? Use this to check when all jobs are done
    num_threads_passive = threading.active_count()
    print(f"\tNumber of passive cores: {num_threads_passive}")

    # Split the table up into two parts, because it is still just too big to all be saved as user tables.
    # Also adds better usability if we ever want to scale up the preprocessing
    # Meta batches might not be the best name, but it's a batch of batches.
    if num_runs == 1:
        meta_batches = batches
    elif num_runs == 2:
        middle_idx = int(len(batches) // 2)
        meta_batches = [batches[:middle_idx], batches[middle_idx:]]
    elif num_runs > 2:
        raise NotImplementedError("Currently cannot cut a list of tables into more than 2 meta batches.")

    counter = 0
    for batches in meta_batches:
        i = counter
        for batch in batches:
            batch_name = f'gaia_astrometric_batch_{i}'
            Gaia.upload_table(upload_resource=batch,
                              table_name=batch_name)  # rkievit user space is increased to 2GB, full dataset is ~1.9GB
            # The threading module won't wait untill the archive job is completed, but instead will send all batches
            # at the same time through propagate_error_one, which saves all tables in ../results/
            threading.Thread(target=propagate_error_one, args=(batch_name, epoch1, epoch2,)).start()
            i += 1  # track batch index through all meta batches

        all_threads_done = False
        while not all_threads_done:
            time.sleep(60)  # Expected processing time @ 1e5 objects ~10m. Set 1m sleep time to avoid overuse
            if threading.active_count() == num_threads_passive:
                print(f"Still processing, number of active cores: {threading.active_count()}", end='\r')
                all_threads_done = True  # Finish the while loop once all tables are saved
        i = counter
        for _ in range(len(batches)):
            batch_name = f'gaia_astrometric_batch_{i}'
            Gaia.delete_user_table(batch_name)
            i += 1
        counter = i  # make sure we continue counting such that there's no overlap


def propagate_error_one(table_name, epoch1, epoch2, save_path='../results'):
    """Query position and error propagation of user table, and save resulting table into a file"""
    query = f"""SELECT  source_id, ra, dec, pmra, pmdec, parallax,
                        array_element(a0, 1) as ra_prop,
                        array_element(a0, 2) as dec_prop,
                        array_element(a0, 3) as parallax_prop,
                        array_element(a0, 4) as pmra_prop,
                        array_element(a0, 5) as pmdec_prop,
                        array_element(a0, 6) as rv_prop,
                        array_element(a1, 1) as e_ra_prop,
                        array_element(a1, 2) as e_de_prop,
                        array_element(a1, 7)as ra_dec_prop,
                        phot_g_mean_mag, bp_rp
                FROM (
                    SELECT  *,
                                EPOCH_PROP(ra,dec,parallax,pmra,pmdec,radial_velocity,{epoch1},{epoch2}) as a0,
                            EPOCH_PROP_ERROR(
                             ra, dec, parallax, pmra, pmdec, radial_velocity,
                             ra_error, dec_error, parallax_error, pmra_error, pmdec_error, radial_velocity_error,
                             ra_dec_corr, ra_parallax_corr, ra_pmra_corr, ra_pmdec_corr,
                             dec_parallax_corr, dec_pmra_corr, dec_pmdec_corr,
                             parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr,
                             {epoch1}, {epoch2}) as a1
                    FROM user_{username}.{table_name}
                    ) as p"""
    job = launch_job(query)
    tab = get_data(job)
    tab.write('../../results/' + table_name + '.fits', format='fits', overwrite=True)


def query_extra_data(tab_xm, extra_data_names, cat='gaia', cat_dpath='../../data/',
                     return_strings=False):
    """Query ancillary data from the archive to be used for the best neighbour selection functions
    Gaia communication is via .fits format so we have to convert np array -> fits and back"""
    import numpy as np

    gaia_ids = np.array(tab_xm[:, 1], dtype=np.int64)
    table = Table([tab_xm[:, 0], gaia_ids],
                  names=('hip', 'source_id'))  # , dtype=(int, int, np.float64))
    table_name = 'xm_table'
    gaia_login()

    # Upload table with error catching failsafe
    i = 0
    while True:
        try:
            Gaia.upload_table(upload_resource=table, table_name=table_name)
            break
        except:
            print(f"Table {table_name} for {username} already exists. Trying again with _{i}..")
            if not table_name[-2] == '_':
                table_name += f'_{i}'
            else:
                table_name = table_name[0:-2] + f'_{i}'
            i += 1

    try: # Remove table if something goes wrong here
        # Query data
        if cat == 'gaia':
            cat_name = 'gaiadr3.gaia_source'
            match_on = 'source_id'
        elif cat == 'hipp':
            cat_name = 'public.hipparcos_newreduction'
            match_on = 'hip'
        elif cat == 'hipp1':
            cat_name = 'public.hipparcos'
            match_on = 'hip'
        else:  # Upload one of the user tables
            # NOTE: This does not work for a full GaiaBaseCat because it is too big.
            #data_table = fits.open(f'{cat_dpath}{cat}.fits')[1]
            data_table = Table.read(f'{cat_dpath}{cat}.fits', format='fits')
            data_table_name = 'GaiaCat'
            cat = 'tab'  # For neatness of the query
            Gaia.upload_table(upload_resource=data_table, table_name=data_table_name)
            cat_name = f'user_{username}.GaiaCat'
            match_on = 'source_id'

        query = "SELECT xm.*"
        for name in extra_data_names:
            query += f', {cat}.{name}'
        query += f"""
                 FROM {cat_name} as {cat}
                 JOIN user_{username}.{table_name} AS xm USING ({match_on})"""

        job = launch_job(query)
        res = get_data(job)
        if return_strings:
            extra_data = res
        else:
            extra_data = np.zeros((tab_xm.shape[0], len(extra_data_names)))#, dtype=np.float64)
            for i, name in enumerate(extra_data_names):
                extra_data[:, i] = res[name]

        # Just in case the Archive shuffled our matches, reassign all hip and source_id values to match the extra data
        # new_tab_xm = np.zeros((tab_xm.shape[0], 2), dtype=np.int64)
        # new_tab_xm[:, 0] = res['hip']
        # new_tab_xm[:, 1] = res['source_id']
        # But it doesn't, and tracking changes to GaiaCat_ID is not doable here, so just error check

        if not np.array_equal(res['hip'], tab_xm[:, 0]):
            raise ValueError("Crossmatches array was shuffled during transfer with the Archive.")
        else:
            print("Post data query array is not shuffled")

    except:
        Gaia.delete_user_table(table_name)
        raise

    # Remove tables
    Gaia.delete_user_table(table_name)
    if cat not in ['gaia', 'hipp', 'hipp1']:
        Gaia.delete_user_table(data_table_name)

    return tab_xm, extra_data


def main():
    # gaia_login()
    # table = query_gaia_simple_conesearch(back_prop_gaia=True, mag_lim=5)
    # return table
    delete_unlabeled_jobs()
    #import numpy as np
    #tab_xm = np.load('../results/crossmatch_1as_all_neighbours.npy')
    #print(query_extra_data(tab_xm, ['pmra', 'pmdec']))


if __name__ == '__main__':
    main()

############################
#
# Complete preprocessing code for the Gaia Catalogue
# Currently can apply magnitude limit, proper motion correction and error inflation
# Additionally currently contains function to create all catalgoues 'simultaneously' (move?)
#
# Callable within full_crossmatch for crossmatching, or seperately to save as .fits file 
#
############################

import numpy as np
from astropy.table import Table, vstack
from astropy.io import fits
from astrometry_equations import sind, cosd
from table_functions import extract_sky_positions, extract_proper_motions, batch_table, read_tables
from data_queries import gaia_login, query_gaia_preprocess, propagate_batches_error, delete_unlabeled_jobs, query_gaia_zpcorr
import time
from numba import njit

@njit
def numba_edr3ToICRF(pmra, pmde, ra, dec, G, table1, apply_color_correction=False):
    """
    Input: source position , coordinates ,
    and G magnitude from Gaia DR3.
    Output: corrected proper motion.
    """

    Gmin = table1[0]
    Gmax = table1[1]

    corrected_pmra = np.zeros(pmra.shape[0])
    corrected_pmde = np.zeros(pmde.shape[0])

    for i in range(pmra.shape[0]):
        if G[i] >= 13:
            corrected_pmra[i] = pmra[i]
            corrected_pmde[i] = pmde[i]
        else:

            # pick the appropriate omegaXYZ for the source’s magnitude:
            omegaX = table1[2][(Gmin <= G[i]) & (Gmax > G[i])][0]
            omegaY = table1[3][(Gmin <= G[i]) & (Gmax > G[i])][0]
            omegaZ = table1[4][(Gmin <= G[i]) & (Gmax > G[i])][0]

            pmraCorr = -1 * np.sin(np.radians(dec[i])) * np.cos(np.radians(ra[i])) * omegaX - np.sin(np.radians(dec[i])) * np.sin(np.radians(ra[i])) * omegaY + np.cos(
                np.radians(dec[i])) * omegaZ
            pmdecCorr = np.sin(np.radians(ra[i])) * omegaX - np.cos(np.radians(ra[i])) * omegaY

            corrected_pmra[i] = pmra[i] - pmraCorr / 1000.
            corrected_pmde[i] = pmde[i] - pmdecCorr / 1000.
        #print(f"PM Corr: {i+1}/{pmra.shape[0]}", end='\r')
    if apply_color_correction:
        corrected_pmra = np.sqrt(corrected_pmra**2 + 10e-3**2)
        corrected_pmde = np.sqrt(corrected_pmde**2 + 10e-3**2)        

    return corrected_pmra, corrected_pmde


def edr3ToICRF(pmra, pmde, ra, dec, G):
    """Prepares the PMC algorithm such that numba can work with it
    numba does not support np.fromstring, so we have to pass it an array instead"""

    table1 = """  0.0 9.0 18.4 33.8 -11.3
                  9.0 9.5 14.0 30.7 -19.4
                  9.5 10.0 12.8 31.4 -11.8
                  10.0 10.5 13.6 35.7 -10.5
                  10.5 11.0 16.2 50.0 2.1
                  11.0 11.5 19.4 59.9 0.2
                  11.5 11.75 21.8 64.2 1.0
                  11.75 12.0 17.7 65.6 -1.9
                  12.0 12.25 21.3 74.8 2.1
                  12.25 12.5 25.7 73.6 1.0
                  12.5 12.75 27.3 76.6 0.5
                  12.75 13.0 34.9 68.9 -2.9 """
    table1 = np.fromstring(table1, sep=' ').reshape((12, 5)).T

    return numba_edr3ToICRF(pmra, pmde, ra, dec, G, table1)

def error_inflation(table, inflated_errors_array=['ra_error', 'dec_error', 'parallax_error', 'pmra_error', 'pmdec_error'],
                    inflation_type='Brandt21'):
    """Inflate the Gaia errors according to methods prescribed in one of multiple papers"""

    if inflation_type == 'Brandt21':
        inflation_value = 1.37  # Brandt 2021
        for val in inflated_errors_array:
            table[val] = table[val] * inflation_value
    else:
        raise ValueError("This inflation method is not known or not yet implemented.")

    return table

def zero_point_corr(tab):
    """Applies the zero point correction to the parallax as discussed by Lindegren+21 and implemented
    by them in the gaiadr3-zeropoint code"""
    from zero_point import zpt
    zpt.load_tables()

    mask = np.isfinite(tab['parallax'] ) # only apply plx correction if we have plx information
    zp_corr = zpt.get_zpt(tab['phot_g_mean_mag'][mask],
                          tab['nu_eff_used_in_astrometry'][mask],
                          tab['pseudocolour'][mask],
                          tab['ecl_lat'][mask],
                          tab['astrometric_params_solved'][mask])
    tab['parallax'][mask] -= zp_corr
    return tab['parallax']

    # query correct data
    #N = len(gaia_ids)

    # 16 milion entries exceeds the 2GB limit on the archive
    #tab1 = query_gaia_zpcorr(gaia_ids[:N//2])
    #tab2 = query_gaia_zpcorr(gaia_ids[N//2:])

    # Append 
    #print(tab1['parallax'], tab2['parallax'])
    #print((tab1['parallax'], tab2['parallax']).flatten())

    #plx = np.append(tab1['parallax'], tab2['parallax'])
    #g_mag = np.append(tab1['phot_g_mean_mag'], tab2['phot_g_mean_mag'])
    #nu_eff = np.append(tab1['nu_eff_used_in_astrometry'], tab2['nu_eff_used_in_astrometry'])
    #pseudo_col = np.append(tab1['pseudocolour'], tab2['pseudocolour'])
    #ecl_lat = np.append(tab1['ecl_lat'], tab2['ecl_lat'])
    #params_solved = np.append(tab1['astrometric_params_solved'], tab2['astrometric_params_solved'])
    #print(plx)
    #mask = np.isfinite(plx) # only apply plx correction if we have plx information
    #zp_corr = zpt.get_zpt(g_mag[mask],
    #                      nu_eff[mask],
    #                      pseudo_col[mask],
    #                      ecl_lat[mask],
    #                      params_solved[mask])
    #plx[mask] -= zp_corr
    #return plx



def full_preprocess(mag_lim=14, gaia_epoch=2016., hipp_epoch=1991.25, batch_size=None, read_local=False, data_path=None,
                    apply_pm_corr=False, error_inflation_type=None, apply_color_correction=False, apply_zp_corr=False,
                    return_cat=False, save_cat=True, savename='GaiaCat_noname'):
    """Completely preprocess the Gaia data from the Archive, into something usable for 
    crossmatching with Hipparcos
    1. Extract all proper Gaia data with M<mag_lim
    2. Apply PM correction (Cantat-Gaudin & Brandt 2021)
    3. Apply zero-point parallax correction (Lindegren et al. 2021)
    4. Apply Error Inflation (Brandt 2021, ..)   
    5. Batch data for processing speed
    6. Apply backpropagation in the archive
    7. Create complete table"""

    gaia_login()
    t1 = time.time()
    print(f"Starting Gaia Data Preprocessing with M_lim = {mag_lim}..")

    # 1.
    if not read_local:
        all_gaia_maglim = query_gaia_preprocess(mag_lim)
    else:
        if 'vot' in data_path:
            all_gaia_maglim = Table.read(data_path)  # Reading full data takes ~1h, way too long just download it.
        elif 'fits' in data_path:
            all_gaia_maglim = fits.open(data_path)[1].data
    print(f"Gaia Archival Data imported. Time elapsed {time.time() - t1:.0f}s \n "
          f"Total Gaia objects: {int(len(all_gaia_maglim))}")

    # 2.
    t2 = time.time()
    if apply_pm_corr:
        ra, dec = extract_sky_positions(all_gaia_maglim)
        pmra, pmde = extract_proper_motions(all_gaia_maglim)
        g_mag = np.array(all_gaia_maglim['phot_g_mean_mag'], dtype=np.float64) # set dtype for numba

        pmra_corr, pmde_corr = edr3ToICRF(pmra, pmde, ra, dec, g_mag)

        all_gaia_maglim['pmra'] = pmra_corr
        all_gaia_maglim['pmdec'] = pmde_corr
        print(f"Proper Motion Correction Applied. Time elapsed {time.time() - t2:.0f}s")

    #3.
    t3 = time.time()
    print('Starting Zero Point Parallax Correction')
    if apply_zp_corr:
        all_gaia_maglim['parallax'] = zero_point_corr(all_gaia_maglim)
        print(f"Zero Point Parallax Correction Applied. Time elapsed {time.time() - t3:.0f}s")

    # 4.
    t4 = time.time()
    if error_inflation_type is not None:
        all_gaia_maglim = error_inflation(all_gaia_maglim, inflation_type=error_inflation_type)
        print(f"Error Inflation with Method {error_inflation_type} Applied. Time elapsed {time.time() - t4:.0f}s")

    # 5.
    gaia_batches = batch_table(all_gaia_maglim, batch_size=batch_size)

    print(f"Gaia data batched into {len(gaia_batches)} separate batches\n"
          f"Starting data propagation..")
    del all_gaia_maglim

    # 6.
    t5 = time.time()
    gaia_login()  # Do this as late as possible to make sure it does not expire
    propagate_batches_error(gaia_batches, gaia_epoch, hipp_epoch)
    print(f"All Gaia data propagated and saved. Propagation time taken: {time.time() - t5:.0f}s")
    del gaia_batches
    
    # 7.
    gaia_batches_prop = read_tables(table_path='../../results/gaia_astrometric_batch_*', multiple=True)
    gaia_tab_prop = vstack(gaia_batches_prop)  # astropy.table function. Should arrange everything automatically

    # Add in this error catcher in case the file already exists to ensure both no runtime is lost, and old
    # catalogues are not automatically deleted
    if save_cat:
        try:
            gaia_tab_prop.write('../../data/'+savename+'.fits', format='fits')
        except OSError:
            saved = False
            while not saved:
                print(f"Cannot write file {savename}.fits into ../../data/. Please insert different name:")
                savename = input()
                try:
                    gaia_tab_prop.write('../../data/' + savename + '.fits', format='fits')
                except OSError:
                    continue
                saved = True  # If no error it should call this
    print(f"Completed. Total runtime: {time.time() - t1:.0f}s")

    if return_cat:
        return gaia_tab_prop


def main():
    mag_lim = 5   # determined with H-G relations
    gaia_epoch = 2016.0
    hipp_epoch = 1991.25
    batch_size = int(135)  # Note: Code currently does not work if only one batch is made
    read_local = False  # UvL Vdesk can't properly download the .vot file. Fix this once on location.

    data_path = '../../data/gaia_process_maglim14.vot'

    #make_all_catalogues(mag_lim, gaia_epoch, hipp_epoch, batch_size, read_local, data_path)
    full_preprocess(mag_lim, gaia_epoch, hipp_epoch, batch_size, savename='GaiaBaseCat5')
    
if __name__ == '__main__':
    main()

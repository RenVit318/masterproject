import numpy as np
from astropy.table import Table, vstack
from astrometry_equations import sind, cosd
from table_functions import extract_sky_positions, extract_proper_motions, batch_table, read_tables
from data_queries import gaia_login, query_gaia_preprocess, propagate_batches_error
import time


def edr3ToICRF(pmra, pmde, ra, dec, G):
    """
    Input: source position , coordinates ,
    and G magnitude from Gaia EDR3.
    Output: corrected proper motion.
    TODO? : Numerical efficiency could be improved by a lot by using numpy instead of for loop!!
    """
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

            pmraCorr = -1 * sind(dec[i]) * cosd(ra[i]) * omegaX - sind(dec[i]) * sind(ra[i]) * omegaY + cosd(dec[i]) * omegaZ
            pmdecCorr = sind(ra[i]) * omegaX - cosd(ra[i]) * omegaY

            corrected_pmra[i] = pmra[i] - pmraCorr / 1000.
            corrected_pmde[i] = pmde[i] - pmdecCorr / 1000.
    return corrected_pmra, corrected_pmde


def full_preprocess(mag_lim, gaia_epoch, hipp_epoch, batch_size=None, read_local=False, data_path=None):
    """Completely preprocess the Gaia data from the Archive, into something usable for 
    crossmatching with Hipparcos
    1. Extract all proper Gaia data with M<mag_lim
    2. Apply PM correction (Cantat-Gaudin & Brandt 2021)
    3. Batch data for processing speed
    4. Apply backpropagation in the archive
    5. Create complete table"""

    gaia_login()
    t0 = time.time()
    print(f"Starting Gaia Data Preprocessing with M_lim = {mag_lim}..")
    
    # 1.
    if not read_local:
        all_gaia_maglim = query_gaia_preprocess(mag_lim)
    else:
        all_gaia_maglim = Table.read(data_path) # Reading full data takes ~1h, way too long just download it.
    t1 = time.time()
    print(f"Gaia Archival Data imported. Time elapsed {t1 - t0:.2f}s \n "
          f"Total Gaia objects: {int(len(all_gaia_maglim))}")

    # 2.
    ra, dec = extract_sky_positions(all_gaia_maglim)
    pmra, pmde = extract_proper_motions(all_gaia_maglim)

    #pmra_corr, pmde_corr = edr3ToICRF(pmra, pmde, ra, dec, all_gaia_maglim['phot_g_mean_mag'])

    #all_gaia_maglim['pmra'] = pmra_corr
    #all_gaia_maglim['pmde'] = pmde_corr
    t2 = time.time()
    print(f"Proper Motion Correction Applied. Time elapsed {t2 - t1:.2f}s")

    # 3.
    #print(all_gaia_maglim)
    gaia_batches = batch_table(all_gaia_maglim, batch_size=batch_size)
    #print(gaia_batches)
    print(f"Gaia data batched into {len(gaia_batches)} separate batches\n"
          f"Starting data propagation..")
    del all_gaia_maglim

    # 4.
    gaia_login() # Do this as late as possible to make sure it does not expire
    propagate_batches_error(gaia_batches, gaia_epoch, hipp_epoch)
    t3 = time.time()
    print(f"All Gaia data propagated and saved. Propagation time taken: {t3 - t2:.2f}s")
    del gaia_batches

    # 5.
    gaia_batches_prop = read_tables(table_path='../results/gaia_astrometric_batch_*', multiple=True)
    gaia_tab_prop = vstack(gaia_batches_prop)  # astropy.table function. Should arrange everything automatically
    gaia_tab_prop.write('../../results/GaiaBaseCat_SIMPLE.fits', format='fits')
    print(f"Completed. Total runtime: {time.time() - t0:.2f}s")


def main():
    mag_lim = 14  # determined with H-G relations
    gaia_epoch = 2016.0
    hipp_epoch = 1991.25
    batch_size = int(1e5) # Note: Code currently does not work if only one batch is made
    read_local = False
    data_path = '../../data/gaia_process_maglim14.vot'

    full_preprocess(mag_lim, gaia_epoch, hipp_epoch, batch_size, read_local, data_path)


if __name__ == '__main__':
    main()

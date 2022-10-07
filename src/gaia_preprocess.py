import numpy as np
from astrometry_equations import sind, cosd

def edr3ToICRF (pmra ,pmdec ,ra ,dec ,G):
    """
    Input: source position , coordinates ,
    and G magnitude from Gaia EDR3.
    Output: corrected proper motion.
    """
    if G >=13:
        return pmra , pmdec
    

    table1="""  0.0 9.0 18.4 33.8 -11.3
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

    Gmin = table1 [0]
    Gmax = table1 [1]

    #pick the appropriate omegaXYZ for the source’s magnitude:
    omegaX = table1 [2][( Gmin <=G)&(Gmax >G)][0]
    omegaY = table1 [3][( Gmin <=G)&(Gmax >G)][0]
    omegaZ = table1 [4][( Gmin <=G)&(Gmax >G)][0]

    pmraCorr = -1* sind(dec)*cosd(ra)*omegaX - sind(dec)*sind(ra)*omegaY + cosd(dec)*omegaZ
    pmdecCorr = sind(ra)*omegaX -cosd(ra)*omegaY

    return pmra -pmraCorr /1000. , pmdec - pmdecCorr /1000.


def full_preprocess(mag_lim, gaia_epoch, hipp_epoch, batch_size=None):
    """Completely preprocess the Gaia data from the Archive, into something usable for 
    crossmatching with Hipparcos
    1. Extract all proper Gaia data with M<mag_lim
    2. Apply PM correction (Cantat-Gaudin & Brandt 2021)
    3. Batch data for processing speed
    4. Apply backpropagation in the archive
    5. Create complete table"""

    #1.
    all_gaia_maglim = query_gaia_preprocess(mag_lim)
    
    #2.
    ra, dec  = extract_sky_positions(all_gaia_maglim)
    pmra, pmde = extract_proper_motions(all_gaia_maglim)

    pmra_corr, pmde_corr = edr3ToICRF(pmra, pmde, ra, dec, all_gaia_maglim['phot_g_mean_mag'])

    all_gaia_maglim['pmra'] = pmra_corr
    all_gaia_maglim['pmde'] = pmde_corr

    #3.
    gaia_batches = batch_table(all_gaia_maglim, batch_size=batch_size)

    #4. 
    propagate_batches_errros(gaia_batches, gaia_epoch, hipp_epoch)

    #5.
    gaia_batches_prop = read_tables(....)
    gaia_tab_prop = combine_tables(gaia_batches_prop)

    save_table(gaia_tab_prop, save_path='../results', save_name='GaiaBaseCat')
    
    #FUNCTIONS TO MAKE
    # gaia_queries: query_gaia_preprocess, propagate_batches_error
    # table_functions: read_tables, combine_tables
    #



def main():
    mag_lim = 14 # determined with H-G relations
    gaia_epoch = 2016.0
    hipp_epoch = 1991.25
    batch_size = 1e5

    full_preprocess(mag_lim, gaia_epoch, hipp_epoch, batch_size)
    #print(edr3ToICRF(1, 2, 1, 2, 10))

if __name__ == '__main__':
    main()

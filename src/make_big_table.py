import numpy as np
from astropy.io import fits

def get_hipp_gaia_data(crossmatch_idxs, HippCat='Hipparcos_mix', GaiaCat='GaiaBaseCat+PMC+EIB'):
    # Hipparcos - Gaia indexes of the cross match 
    hipp_cat = fits.open(f'{dpath}{HippCat}.fits')[1].data
    gaia_cat = fits.open(f'{dpath}{GaiaCat}.fits')[1].data

    # Make a table matching the XM order containing all propagated Gaia astrometry
    pos_gaia_table_full = copy.deepcopy(gaia_cat)
    pos_gaia_table_full = pos_gaia_table_full[crossmatch_idxs[:,2]]

    # Indexing this table is slightly more work because we do not have the table idxs stored
    pos_hipp_table_full = copy.deepcopy(hipp_cat)

    hip_indexes = np.zeros(crossmatch_idxs.shape[0], dtype=np.int64)
    for i in range(crossmatch_idxs.shape[0]):
        hip_indexes[i] = int(np.where(hipp_cat['hip'] == crossmatch_idxs[i][0])[0][0])
    pos_hipp_table_full = pos_hipp_table_full[hip_indexes]

    return pos_hipp_table_full, pos_gaia_table_full



def make_table():
    # GET ALL DATA #


    ###

    # Gpred - Gmatch
    
def main():
    make_table()

if __name__ == '__main__':
    main()

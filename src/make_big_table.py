import numpy as np
from astropy.io import fits
import copy
from astrometry_equations import predict_G_bv, compute_error_normalized_distance

def get_hipp_gaia_data(crossmatch_idxs, HippCat='Hipparcos_mix', GaiaCat='GaiaBaseCat+PMC+EIB', dpath='../../data/'):
    # Hipparcos - Gaia indexes of the cross match 
    hipp_cat = fits.open(f'{dpath}{HippCat}.fits')[1].data
    gaia_cat = fits.open(f'{dpath}{GaiaCat}.fits')[1].data
    print('start gaia')
    # Make a table matching the XM order containing all propagated Gaia astrometry
    pos_gaia_table_full = copy.deepcopy(gaia_cat)
    pos_gaia_table_full = pos_gaia_table_full[crossmatch_idxs[:,2]]
    print('start hipp')
    # Indexing this table is slightly more work because we do not have the table idxs stored
    pos_hipp_table_full = copy.deepcopy(hipp_cat)

    hip_indexes = np.zeros(crossmatch_idxs.shape[0], dtype=np.int64)
    for i in range(crossmatch_idxs.shape[0]):
        hip_indexes[i] = int(np.where(hipp_cat['hip'] == crossmatch_idxs[i][0])[0][0])
    pos_hipp_table_full = pos_hipp_table_full[hip_indexes]

    return pos_hipp_table_full, pos_gaia_table_full


def make_table():
    dpath = '../../data/'
    rpath = '../results/'
    GaiaCat = 'GaiaBaseCat_all.fits'
    HippCat = 'Hipparcos_mix.fits'
    
    # Fields that we want to add
    fieldnames = ['G_Gpred', 'Hip_BV', 'Gaia_BpRp', 'D', 'distance', 'delta_pm_alpha', 'delta_pm_dec', 'delta_pm_tot']
    
    # GET ALL DATA #
    pos_xm = np.load(rpath+'final_crossmatch+PMC+EIB_10as__best_neighbour_likeliest_position.npy')
    mag_xm = np.load(rpath+'final_crossmatch+PMC+EIB_10as__best_neighbour_likeliest_magnitude.npy')
    nea_xm = np.load(rpath+'final_crossmatch+PMC+EIB_10as_best_neighbour_nearest.npy')
    print(np.array_equal(pos_xm[:,0], mag_xm[:,0]))
    print(np.array_equal(pos_xm[:,0], nea_xm[:,0]))
    print(np.array_equal(mag_xm[:,0], nea_xm[:,0]))
    input()
    hipp_pos, gaia_pos = get_hipp_gaia_data(pos_xm)
    hipp_mag, gaia_mag = get_hipp_gaia_data(mag_xm)
    hipp_nea, gaia_nea = get_hipp_gaia_data(nea_xm)
    ###

    xm_tabs = [pos_xm, mag_xm, nea_xm]
    hipp_dat = [hipp_pos, hipp_mag, hipp_nea]
    gaia_dat = [gaia_pos, gaia_mag, gaia_nea]

    Nsel = len(xm_tabs)
    # Use these to place the right data in the right place of the table
    tab_idxs = [0]
    cum_len = 0
    for i in range(Nsel):
        tab_idxs.append(cum_len+len(xm_tabs[i]))
        cum_len += len(xm_tabs[i])
    print(tab_idxs)

    print('Starting ID Extraction..')
    table_columns = []
    # Get sorting indices based on hip such that we have all hip indices in order
    # Get the indices separately because we have to save gaia indices separately
    hip = np.array([])
    gaiaid = np.array([])
    for i in range(Nsel):
        hip = np.append(hip, xm_tabs[i][:,0])
        gaiaid = np.append(gaiaid, xm_tabs[i][:,1])
    sort_idxs = np.argsort(hip)

    table_columns.append(fits.Column(name='hip_id', format='K', array=hip[sort_idxs]))
    table_columns.append(fits.Column(name='gaia source_id', format='K', array=gaiaid[sort_idxs]))

    table_data = np.zeros((hip.shape[0], len(fieldnames)), dtype=np.float64)
    
    print('Starting Parameter Extraction..')
    # Extract all variables from the tables
    for i in range(Nsel):
        sub_tab = table_data[tab_idxs[i]:tab_idxs[i+1], :] # pointer towards correct part of table
        hipp = hipp_dat[i]
        gaia = gaia_dat[i]
        
        for j in range(len(fieldnames)):
            # The below kind of works like a long list of elif statements
            match fieldnames[j]:
                case 'G_Gpred':
                    G_pred, _ = predict_G_bv(hipp['hp_mag'], hipp['b_v'])
                    sub_tab[:, j] = gaia['phot_g_mean_mag'] - G_pred
                case 'D':
                    pos_hipp = np.array([hipp['ra'], hipp['dec']]).T * 3.6e6
                    pos_gaia = np.array([gaia['ra_prop'], gaia['dec_prop']]).T * 3.6e6
                    unc_hipp = np.array([hipp['e_ra_rad']**2, hipp['e_de_rad']**2, hipp['ra_dec_corr']]).T
                    unc_gaia = np.array([gaia['e_ra_prop']**2, gaia['e_de_prop']**2, gaia['ra_dec_prop']]).T
                    sub_tab[:, j] = compute_error_normalized_distance(pos_hipp, pos_gaia, unc_hipp, unc_gaia, method='full')                   
                case 'distance':
                    pos_hipp = np.array([hipp['ra'], hipp['dec']]).T * 3.6e6
                    pos_gaia = np.array([gaia['ra_prop'], gaia['dec_prop']]).T * 3.6e6
                    unc_hipp = np.array([hipp['e_ra_rad']**2, hipp['e_de_rad']**2, hipp['ra_dec_corr']]).T
                    unc_gaia = np.array([gaia['e_ra_prop']**2, gaia['e_de_prop']**2, gaia['ra_dec_prop']]).T
                    sub_tab[:, j] = compute_error_normalized_distance(pos_hipp, pos_gaia, unc_hipp, unc_gaia, method='none')   
                case 'delta_pm_alpha':
                    sub_tab[:, j] = np.abs(hipp['pm_ra'] - gaia['pmra_prop'])
                case 'delta_pm_dec':
                    sub_tab[:, j] = np.abs(hipp['pm_de'] - gaia['pmdec_prop'])
                case 'delta_pm_tot':
                    sub_tab[:, j] = np.sqrt((hipp['pm_ra'] - gaia['pmra_prop'])**2 + (hipp['pm_de'] - gaia['pmdec_prop'])**2)
                # maybe just use the catalogue names for everything below and assign them all automatically
                case 'Hip_BV':
                    sub_tab[:, j] = hipp['b_v']
                case 'Gaia_BpRp':
                    sub_tab[:, j] = gaia['bp_rp']
            print(table_data)     

    # Now turn all these tables into fits columns (can we do this in one big go from a ndarray?
    for j in range(table_data.shape[1]):
        table_columns.append(fits.Column(name=fieldnames[j], format='D', array=table_data[:,j]))

    print('Starting HDUList Creation and Saving..')
    table_hdu = fits.BinTableHDU.from_columns(table_columns)

    table_hdu.writeto(rpath+'all_results_10as.fits', overwrite=True)

    
def main():
    make_table()


if __name__ == '__main__':
    main()

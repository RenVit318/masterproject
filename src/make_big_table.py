import numpy as np
from astropy.io import fits
import copy
from astrometry_equations import predict_G_bv, compute_error_normalized_distance
from data_queries import query_extra_data

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


def get_delta_logG(xm, 
                   lumclass_logG_dict = { # Based on slides. Find real citation!
                        'I': -0.5,
                        'II': 0.5,
                        'III': 1.5,
                        'IV': 3.,
                        'V': 4.5,
                        'VI': 5.5,
                        'VII': 8. }): # PROBABLY NOT GOOD!
    """Compute the difference in log Surface Gravity [cgs] between the Hipparcos and Gaia matches. This is a proxy
    for the luminosity class of an object. The Gaia catalogue directly reports this number based on ..., for the 
    Hipparcos catalogue we need to do some work and make some assumptions"""
    from MeanStars import MeanStars
    ms = MeanStars()
    # Gaia              
    _, logG_gaia = query_extra_data(xm, ['logg_gspphot'])
    # assume unclassified objects are MS stars/dwarfs
    print(len( logG_gaia[np.isnan(logG_gaia)]))
    logG_gaia[np.isnan(logG_gaia)] = lumclass_logG_dict['V'] 
    # Hipp 
    _, tab_hipp = query_extra_data(xm, ['sptype'], cat='hipp1', return_strings=True)
    spectype_hipp = tab_hipp['sptype']
    logG_hipp = np.zeros(len(spectype_hipp), dtype=np.float64)
    # annoyingly match_spec only works on one string at a time.
    for i in range(len(spectype_hipp)):
        
        SpecType = ms.matchSpecType(spectype_hipp[i])
        if SpecType is None: # No Spectral Type. Assume it is a dwarf
            logG_hipp[i] = lumclass_logG_dict['V']
        elif len(SpecType) > 3: # There were multiple spectral types. None in our dataset
            logG_sum = 0
            for i in range(2, len(SpecType), 3): # select each third element and take the mean
                logG_sum += lumclass_logG_dict[SpecType[i]]
            logG_hipp[i] = logG_sum / (len(SpecType)/3)

        else: # Normally behaving object
            _, _, lclass = SpecType
            logG_hipp[i] = lumclass_logG_dict[lclass]

    print(logG_gaia, logG_hipp)
    print(logG_gaia.shape, logG_hipp.shape)
    #return logG_gaia.T - logG_hipp # Transpose because logG_gaia is vertical
    return logG_gaia, logG_hipp


def make_table():
    dpath = '../../data/'
    rpath = '../results/'
    GaiaCat = 'GaiaBaseCat_Final'
    HippCat = 'Hipparcos_mix'
    
    # Fields that we want to add
    # small table
    #fieldnames = ['G_Gpred', 'Hip_BV', 'Gaia_BpRp', 'D', 'distance', 'delta_pm_alpha', 'delta_pm_dec', 'delta_pm_tot']
    #savename = 'all_results_10as_small'
    # big table
    fieldnames = ['method', 'G_Gpred', 'Hp_mag', 'G_mag', 'Hip_BV', 'Gaia_BpRp', 'D', 'distance', 'delta_pm_alpha', 'delta_pm_dec', 'delta_pm_tot', 'delta_pm_angle', 'delta_plx']
    savename = 'all_neighbours_ml_data_final'
    select_funcs = False

    # GET ALL DATA #
    if select_funcs:
        pos_xm = np.load(rpath+'final_crossmatch+PMC+EIB_10as__best_neighbour_likeliest_position.npy')
        mag_xm = np.load(rpath+'final_crossmatch+PMC+EIB_10as__best_neighbour_likeliest_magnitude.npy')
        nea_xm = np.load(rpath+'final_crossmatch+PMC+EIB_10as_best_neighbour_nearest.npy')

        hipp_pos, gaia_pos = get_hipp_gaia_data(pos_xm, HippCat=HippCat, GaiaCat=GaiaCat)
        hipp_mag, gaia_mag = get_hipp_gaia_data(mag_xm, HippCat=HippCat, GaiaCat=GaiaCat)
        hipp_nea, gaia_nea = get_hipp_gaia_data(nea_xm, HippCat=HippCat, GaiaCat=GaiaCat)
        ###

        xm_tabs = [pos_xm, mag_xm, nea_xm]
        hipp_dat = [hipp_pos, hipp_mag, hipp_nea]
        gaia_dat = [gaia_pos, gaia_mag, gaia_nea]
    else:
        xm = np.load(rpath+'final_crossmatch_10as_all_neighbours.npy')
        hipp, gaia = get_hipp_gaia_data(xm, HippCat=HippCat, GaiaCat=GaiaCat)

        xm_tabs = [xm]
        hipp_dat = [hipp]
        gaia_dat = [gaia]

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
    hip = np.array([], dtype=np.int64)
    gaiaid = np.array([], dtype=np.int64)
    for i in range(Nsel):
        hip = np.append(hip, xm_tabs[i][:,0])
        gaiaid = np.append(gaiaid, xm_tabs[i][:,1])
    sort_idxs = np.argsort(hip)

    print(gaiaid[sort_idxs][0])
    table_columns.append(fits.Column(name='hip_id', format='i8', array=hip[sort_idxs]))
    table_columns.append(fits.Column(name='gaia source_id', format='i8', array=gaiaid[sort_idxs]))
#    print(table_columns[1].data)
 #   input()
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
                    sub_tab[:, j] = gaia['pmra_prop'] - hipp['pm_ra']
                case 'delta_pm_dec':
                    sub_tab[:, j] = gaia['pmdec_prop'] - hipp['pm_de']
                case 'delta_pm_tot':
                    sub_tab[:, j] = np.sqrt(gaia['pmra_prop']**2 + gaia['pmdec_prop']**2) - np.sqrt(hipp['pm_ra']**2 + hipp['pm_de']**2)
                case 'delta_pm_angle':
                    hipp_pm_len = np.sqrt(hipp['pm_ra']**2 + hipp['pm_de']**2)  
                    gaia_pm_len = np.sqrt(gaia['pmra_prop']**2 + gaia['pmdec_prop']**2)
                    dot_prod = hipp['pm_ra'] * gaia['pmra_prop'] + hipp['pm_de'] * gaia['pmdec_prop']
                    sub_tab[:, j] = np.arccos(dot_prod / (hipp_pm_len * gaia_pm_len)) # takes on values between 0 and pi
                case 'delta_plx': 
                    sub_tab[:, j] = gaia['parallax'] - hipp['plx']
                case 'method':
                    sub_tab[:,j] = i
                case 'Delta_LogG':
                    sub_tab[:,j] = get_delta_logG(xm_tabs[i])
                # maybe just use the catalogue names for everything below and assign them all automatically
                case 'Hp_mag':
                    sub_tab[:, j] = hipp['hp_mag']
                case 'Hip_BV':
                    sub_tab[:, j] = hipp['b_v']
                case 'G_mag':
                    sub_tab[:, j] = gaia['phot_g_mean_mag']
                case 'Gaia_BpRp':
                    sub_tab[:, j] = gaia['bp_rp']  

    # Now turn all these tables into fits columns (can we do this in one big go from a ndarray?
    for j in range(table_data.shape[1]):
        table_columns.append(fits.Column(name=fieldnames[j], format='D', array=table_data[:,j][sort_idxs]))

    print('Starting HDUList Creation and Saving..')
    table_hdu = fits.BinTableHDU.from_columns(table_columns)

    table_hdu.writeto(rpath+savename+'.fits', overwrite=True)

    
def main():
    make_table()
    #pos_xm = np.load('../results/final_crossmatch+PMC+EIB_10as__best_neighbour_likeliest_position.npy')
    #get_delta_logG(pos_xm)


if __name__ == '__main__':
    main()

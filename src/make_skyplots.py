################
#
# Script to make all skyplots of Hipparcos - Gaia matches including
# their corresponding error ellipses, and error normalized distance
#
################

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from data_queries import get_hipp_gaia_data
from astrometry_equations import error_ellipses, compute_error_normalized_distance, predict_G_bv


def main():
    # Plotting arguments
    ellipse_kwargs = {
        'fc': None,
        'fill': False
    }
    deg_to_mas = 3600*1e3

    # Crossmatch indexes and all Hipparcos / Gaia data
    xm = np.load('../results/final_crossmatch_10as_all_neighbours.npy')#Q[:100]
    hipp_table, gaia_table = get_hipp_gaia_data(xm, HippCat='Hipparcos_mix', GaiaCat='GaiaBaseCat_Final')
    #ml_table = fits.open('../results/all_neighbours_ml_data_final.fits')

    # Get all <117,955 unique Hipparcos stars in the xm 
    uniques, idx_start, counts = np.unique(xm[:,0], return_index=True, return_counts=True)
    
    print(np.where(uniques>2))

    for i in range(len(uniques)):
        fig, ax = plt.subplots(1,1,figsize=(5,5))
        
        # Hipparcos
        # We can just use idx_start because all hip are the same anyway
        hipp = hipp_table[idx_start[i]]
        h_ra, h_dec = hipp['ra']*deg_to_mas, hipp['dec']*deg_to_mas

        # Make Hipp error circle
        h_covmat = [hipp['e_ra_rad']**2, hipp['e_de_rad']**2,
                    hipp['ra_dec_corr']*hipp['e_ra_rad']*hipp['e_de_rad']]
        h_ellipse = error_ellipses([0,0], h_covmat, [1], **ellipse_kwargs, ec='black')

        ax.scatter(0, 0, c='black', s=15)
        ax.add_patch(h_ellipse[0])
       
        # Gaia Star(s)
        # Select all indexes corresponding to neigbhours of hipp[i]
        gaia_idxs = [idx_start[i] + j for j in range(counts[i])]
        for j, idx in enumerate(gaia_idxs):
            gaia = gaia_table[idx]
            g_ra, g_dec = gaia['ra_prop']*deg_to_mas, gaia['dec_prop']*deg_to_mas
            g_covmat = [gaia['e_ra_prop']**2, gaia['e_de_prop']**2, 
                        gaia['ra_dec_prop']*gaia['e_ra_prop']*gaia['e_de_prop']] 

            g_ellipse = error_ellipses([g_ra-h_ra, g_dec-h_dec], g_covmat, [1], **ellipse_kwargs, ec=f'C{j}')
        
            # Compute some statistics between the Hipp and Gaia star.
            # angular separation in mas, error norm. distance and G-Gpred
            # Need to redefine covariance because of a differing definition with Anthony :(
            g_cov = np.array([gaia['e_ra_prop']**2, gaia['e_de_prop']**2, gaia['ra_dec_prop']])
            h_cov = np.array([hipp['e_ra_rad']**2, hipp['e_de_rad']**2, hipp['ra_dec_corr']])
            sep = compute_error_normalized_distance(np.array([h_ra, h_dec]), np.array([g_ra, g_dec]), h_cov, g_cov, method='none')
            D = compute_error_normalized_distance(np.array([h_ra, h_dec]), np.array([g_ra, g_dec]), h_cov, g_cov, method='full')
            Gpred, _ = predict_G_bv(hipp['hp_mag'], hipp['b_v']) 
            G_Gpred = gaia['phot_g_mean_mag'] - Gpred

            label = rf'Gaia {gaia["source_id"]}'+'\n'+rf'$\rho$={sep:.1f}; $D$={D:.1f}; $\Delta$G={G_Gpred:.1f}'

            ax.scatter(g_ra-h_ra, g_dec-h_dec, c=f'C{j}', label=label, s=15)   
            ax.add_patch(g_ellipse[0])

        # Labels and title
        ax.set_xlabel(r'$\Delta \alpha$ [mas]')
        ax.set_ylabel(r'$\Delta \delta$ [mas]')
        ax.set_title(f'Hip {uniques[i]}')
        

        plt.legend()

        if counts[i] == 3:
            plt.show()
        plt.savefig(f'../../dataproducts/skyplots/hip{uniques[i]}', bbox_inches='tight')
        

if __name__ == '__main__':
    main()

############################
#
# Code to check how similar our crossmatches on the 99,525 Marrese objects
# are to the crossmatch results from Marrese et al. (2018)
#
############################

import numpy as np
from data_queries import launch_job, get_data

def main():
    cs_radii = [1, 2, 3, 5, 10]
    selectors = ['nearest', 'likeliest_position', 'likeliest_magnitude']
    basename = '../results/marrese_crossmatch_'    

    # Get Marrese data
    query = "select source_id, original_ext_source_id from gaiadr3.hipparcos2_best_neighbour order by original_ext_source_id asc"
    job = launch_job(query)
    tab = get_data(job)
    marrese_hipp = tab['original_ext_source_id']
    marrese_gaia = tab['source_id']

    for r in cs_radii:
        print(f'\nStarting with all {r}as crossmatch results')
        all_matches = np.load(basename+f'{r}as_all_neighbours.npy').shape[0]
        print(f'In total we found {all_matches} neighbours')
        for func in selectors:
            my_xm = np.load(basename+f'{r}as_best_neighbour_{func}.npy')
            marrese_match = 0
            not_in_our_xm = 0
            print(f'Selector: {func}')
            print(f'\tNum. Unique: {my_xm.shape[0]}')

            for i in range(len(marrese_hipp)):
                xm_idx = np.where(my_xm[:,0] == marrese_hipp[i])[0]
                if len(xm_idx) > 0: # hip id in both xm cats
                    if my_xm[xm_idx[0]][1] == marrese_gaia[i]:
                        marrese_match += 1
                else: # we didn't find a match for this objects
                    not_in_our_xm += 1

            print(f'\tMatch with Marrese: {marrese_match} ({100*marrese_match/len(marrese_hipp):.2f}%)')
            print(f'\tNot Matched: {not_in_our_xm}')   




if __name__ == '__main__':
    main()

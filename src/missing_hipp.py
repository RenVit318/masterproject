import numpy as np
from astropy.io import fits
from astroquery.gaia import Gaia
from data_queries import query_extra_data


def find_indexes()
    xm = np.load('../results/final_crossmatch_10as_best_neighbour_nearest.npy')
    hip = fits.open('../../data/Hipparcos_mix.fits')[1].data

    # Make an array with in the first column the Hip IDs missing from our crossmatch
    # and in the second column the corresponding hip cat ids
    diff = int(len(hip['hip']) - xm.shape[0])
    no_match = np.zeros((diff, 2))
    j = 0
    for i in range(len(hip['hip'])):
        hid = hip['hip'][i]
        if hid not in xm[:,0]:
            no_match[j][0] = hid
            no_match[j][1] = i
            j += 1  
    np.savetxt('../results/missing_hipp.txt', no_match)


def main():
    #find_indexes()
    analyze()
    

if __name__ == '__main__':
    main()

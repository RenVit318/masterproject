import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from plotting_functions import set_styles


def main():
    set_styles()
    tab_name = 'all_results_10as.fits'
    hdul = fits.open('../results/'+tab_name)
    header = hdul[1].header
    data = hdul[1].data
    
    names = []
    for i in range(len(data[0])):
        names.append(header[f'TTYPE{i+1}'])
    do_log = [False, False, False, False, False, True, True, True, True, True]

    # Data adjustments
    data['Gaia_BpRP'][data['Gaia_BpRP'] > 1e15] = np.nan

    # make a ton of scatterplots
    for i in range(2,len(names)): # start at 2 because the first two field are just IDs
        xname = names[i]        
        x = data[xname]
        
        if do_log[i]:
            x = np.log10(x)
            xname = r'$^{10}\log$ ' + xname
        for j in range(i+1, len(names)):
            yname = names[j]
            y = data[yname]
            if do_log[j]:
                y = np.log10(y)
                yname = r'$^{10}\log$ ' + yname

            # Plotting
            fig, ax = plt.subplots(1,1)
            dist = ax.hexbin(x, y, cmap='viridis', mincnt=1, bins='log')
            cb = plt.colorbar(dist)
            
            cb.set_label('Log Count')
            ax.set_xlabel(xname)
            ax.set_ylabel(yname)
            ax.set_title(xname + ' vs. ' + yname)
            
            plt.savefig('../../figures/scatterplots/hexbin_'+xname+'_'+yname, bbox_inches='tight')
            plt.close()


if __name__ == '__main__':
    main()

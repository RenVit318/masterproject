import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


def plot_star_square(table, center_coords, box_size,
                     ra_identifier='ra', dec_identifier='dec', mag_identifier='phot_g_mean_mag'):
    """Extensive catch-all plotting functions to show stars in the sky in a certain region
    Inputs: table: some kind of table with named columns containing all information on the stars to plot
            center_coords: center of coordinates to plot in degrees
            box_size: size of the box plot in arc-seconds
    TODO: - Add error plotting, PM plotting, Brightness dependent plotting, Colour Plotting"""

    # Selection of which stars to plot within the chosen interval
    ra_mid, de_mid = center_coords[0], center_coords[1]
    ra = table[ra_identifier]
    de = table[dec_identifier]
    size = box_size / 3600.  # box size in degrees
    #star_idxs = #np.where((ra < ra_mid + size) and (ra > ra_mid - size) and (dec < de_mid + size) and (dec > de_mid - size))
    ra_lims = (ra<(ra_mid+size))*(ra>(ra_mid-size))
    de_lims = (de<(de_mid+size))*(de>(de_mid-size))
    star_mask = ra_lims*de_lims

    ra = ra[star_mask]
    de = de[star_mask]
    mag = table[mag_identifier][star_mask]
    flux = 10.**(-mag/2.5)
    print(f"Stars selected - plotting {len(ra)} stars..")

    print(flux)
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    ax.scatter(ra, de, s=3e3*flux, c='gold')

    ax.set_xlim(ra_mid-size, ra_mid+size)
    ax.set_ylim(de_mid-size, de_mid+size)
    ax.invert_xaxis()
    ax.set_xlabel("RA (J1991.25) [Deg]")
    ax.set_ylabel("Dec (J1991.25) [Deg]")

    plt.show()


def main():
    table = fits.open('../../data/gaia_stars_sel12_noerr.fits')[1].data

    center_coords = [82.5, -2]  # degrees
    box_size = 3*3600  # arcseconds
    plot_star_square(table, center_coords, box_size,
                     ra_identifier='ra_prop', dec_identifier='dec_prop')


if __name__ == '__main__':
    main()

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.patches import Ellipse

def plot_star_square(table, center_coords, box_size,
                     plot_color=False, plot_errors=None,
                     ra_identifier='ra', dec_identifier='dec',
                     era_identifier='e_ra_prop', ede_identifier='e_de_prop',
                     corr_rade_identifier='ra_dec_prop',
                     mag_identifier='phot_g_mean_mag', color_identifier='bp_rp'):
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

    ra_lims = (ra < (ra_mid + size)) * (ra > (ra_mid - size))
    de_lims = (de < (de_mid + size)) * (de > (de_mid - size))
    star_mask = ra_lims * de_lims

    ra = ra[star_mask]
    de = de[star_mask]

    # Magnitude symbol-size dependency
    mag = table[mag_identifier][star_mask]
    flux = 10 ** (-mag / 2.5)
    symbol_size = 4e6 * flux
    print(symbol_size)
    symbol_size = np.clip(symbol_size, 0.1, 50)

    print(f"Stars selected - plotting {len(ra)} stars..")
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    if plot_color:
        bp_rp = table[color_identifier][star_mask]
        bp_rp = np.clip(bp_rp, -5, 5)
        stars = ax.scatter(ra, de, s=symbol_size, c=bp_rp, cmap='coolwarm', vmin=-5, vmax=5)
        plt.colorbar(stars)
    else:
        # Standard Plotting
        ax.scatter(ra, de, s=symbol_size, c='gold')

    if plot_errors is not None:
        e_ra = table[era_identifier][star_mask] / (1e3 * 3600.) # RA_error is given in mas
        e_de = table[ede_identifier][star_mask] / (1e3 * 3600.)
        corr_rade = table[corr_rade_identifier]
        for idx in range(len(ra)):
            e_circle = Ellipse([ra[idx], de[idx]], e_ra[idx], e_de[idx], corr_rade[idx],
                                   ec='black', fc=None, fill=False)
            ax.add_patch(e_circle)

    ax.set_xlim(ra_mid - size, ra_mid + size)
    ax.set_ylim(de_mid - size, de_mid + size)
    ax.invert_xaxis()
    ax.set_xlabel("RA (J1991.25) [Deg]")
    ax.set_ylabel("Dec (J1991.25) [Deg]")

    plt.show()


def main():
    table = fits.open('../../data/GaiaBaseCat+PMC+EIB.fits')[1].data
    tab_idx = np.random.choice((table['source_id'][table['phot_g_mean_mag'] < 10]).shape[0])
    center_coords = [table['ra'][tab_idx], table['dec'][tab_idx]]
    #center_coords = [82.5, -2]  # degrees
    box_size = 5*3600  # arcseconds
    plot_star_square(table, center_coords, box_size, plot_color=True, plot_errors=True,
                     ra_identifier='ra_prop', dec_identifier='dec_prop')


if __name__ == '__main__':
    main()

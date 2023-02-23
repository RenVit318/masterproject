import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.patches import Ellipse
import matplotlib as mpl

def set_styles():
    # plt.style.available
    plt.style.use('default')
    mpl.rcParams['axes.grid'] = True
    mpl.rcParams['lines.linewidth'] = 1.5
    mpl.rcParams['font.size'] = 11
    mpl.rcParams['font.family'] = 'serif'
    plt.style.use('seaborn-darkgrid')

def plot_star_square(table, center_coords, box_size,
                     plot_color=False, plot_errors=None,
                     errors_in_rad=False,  # Specifically for Hipparcos data
                     ra_identifier='ra', dec_identifier='dec',
                     era_identifier='e_ra_prop', ede_identifier='e_de_prop',
                     corr_rade_identifier='ra_dec_prop',
                     mag_identifier='phot_g_mean_mag', color_identifier='bp_rp',
                     ax=None):
    """Extensive catch-all plotting functions to show stars in the sky in a certain region
    Inputs: table: some kind of table with named columns containing all information on the stars to plot
            center_coords: center of coordinates to plot in degrees
            box_size: size of the box plot in arc-seconds
    TODO: - Add PM plotting"""

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

    symbol_size = np.clip(symbol_size, 0.1, 50)

    print(f"Stars selected - plotting {len(ra)} stars..")
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    if plot_color:
        bp_rp = table[color_identifier][star_mask]
        bp_rp = np.clip(bp_rp, -5, 5)
        stars = ax.scatter(ra, de, s=symbol_size, c=bp_rp, cmap='coolwarm', vmin=-5, vmax=5)
        plt.colorbar(stars)
    else:
        # Standard Plotting
        ax.scatter(ra, de, s=symbol_size, c='gold')

    # ERROR PLOTTING #
    if plot_errors is not None:
        e_ra = table[era_identifier][star_mask] / (1e3 * 3600.)
        e_de = table[ede_identifier][star_mask] / (1e3 * 3600.)

        if plot_errors == 'full':
            corr_rade = table[corr_rade_identifier]
            for idx in range(len(ra)):
                e_circle = Ellipse([ra[idx], de[idx]], e_ra[idx], e_de[idx], corr_rade[idx],
                                   ec='black', fc=None, fill=False)
                ax.add_patch(e_circle)
        elif plot_errors == 'ellipse':
            for idx in range(len(ra)):
                e_circle = Ellipse([ra[idx], de[idx]], e_ra[idx], e_de[idx],
                                   ec='black', fc=None, fill=False)
                ax.add_patch(e_circle)

    ax.set_xlim(ra_mid - size, ra_mid + size)
    ax.set_ylim(de_mid - size, de_mid + size)
    ax.invert_xaxis()
    ax.set_xlabel("RA (J1991.25) [Deg]")
    ax.set_ylabel("Dec (J1991.25) [Deg]")

    return ax

# CODE FROM A.G.A. Brown
def error_ellipses(mu, covmat, sigma_levels, **kwargs):
    """
    Given a covariance matrix for a 2D Normal distribution calculate the uncertainty-ellipses and return
    matplotlib patches for plotting them.
    Parameters
    ----------
    mu : float array
        Mean of Normal distribution (2-vector)
    covmat : float array
        Covariance matrix stored as [sigma_x^2, sigma_y^2, sigma_xy]
    sigma_levels : float or 1-D array
        Equivalent n-sigma levels to draw
    Returns
    -------
    patches : list of matplotlib.patches.Ellipse
        List of matplotlib.patches.Ellipse objects
    Other parameters
    ----------------
    **kwargs :
        Extra arguments for matplotlib.patches.Ellipse
    """
    import matplotlib as mpl
    from scipy.special import erf

    sigmaLevels2D = -2.0 * np.log(
        1.0 - erf(np.array([sigma_levels]).flatten() / np.sqrt(2.0))
    )

    eigvalmax = 0.5 * (
        covmat[0]
        + covmat[1]
        + np.sqrt((covmat[0] - covmat[1]) ** 2 + 4 * covmat[2] ** 2)
    )
    eigvalmin = 0.5 * (
        covmat[0]
        + covmat[1]
        - np.sqrt((covmat[0] - covmat[1]) ** 2 + 4 * covmat[2] ** 2)
    )
    angle = np.arctan2((covmat[0] - eigvalmax), -covmat[2]) / np.pi * 180
    errEllipses = []
    for csqr in sigmaLevels2D:
        errEllipses.append(
            mpl.patches.Ellipse(
                mu,
                2 * np.sqrt(csqr * eigvalmax),
                2 * np.sqrt(csqr * eigvalmin),
                angle=angle,
                **kwargs
            )
        )

    return errEllipses


def main():
    table_gaia = fits.open('../../data/GaiaBaseCat+PMC+EIB.fits')[1].data
    table_hipp = fits.open('../../data/Hipparcos2.fits')[1].data

    # Select a star from the cross-match table
    xm_tab = np.load('../results/crossmatch+PMC+EIB_2as_all_neighbours.npy')
    tab_idx = np.where(table_gaia['source_id'] == np.random.choice(xm_tab[:, 1]))

    # tab_idx = np.random.choice((table['source_id'][table['phot_g_mean_mag'] < 10]).shape[0])
    center_coords = [table_gaia['ra'][tab_idx], table_gaia['dec'][tab_idx]]
    box_size = 5  # arcseconds

    # Plotting. Maybe convert the identifier arguments into two dictionaries
    ax = plot_star_square(table_gaia, center_coords, box_size, plot_color=True, plot_errors='full',
                          ra_identifier='ra_prop', dec_identifier='dec_prop')

    plot_star_square(table_hipp, center_coords, box_size, plot_color=False, plot_errors='ellipse',
                     era_identifier='e_ra_rad', ede_identifier='e_de_rad',
                     mag_identifier='hp_mag', color_identifier='b_v', ax=ax)

    plt.show()


if __name__ == '__main__':
    main()

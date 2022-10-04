from astroquery.gaia import Gaia

def gaia_login(user='rkievit', password=""):
    Gaia.login(user=user)

def read_gaia_hipp_data(gaia_path, hipp_path,
                        num_hipp='all',
                        **kwargs):
    """Read Gaia and Hipparcos data for cross-match purposes. """
    tab_g = fits.open(gaia_path)
    tab_h = fits.open(hipp_path)

    if not num_hipp == 'all':
        # NOTE: When cropping FITS Table like this, some header vals might be invalid
        tab_h[1].data = tab_h[1].data[:num_hipp]

    return tab_h[1].data, tab_g[1].data


def query_gaia_simple_conesearch(mag_lim=None,
                                 back_prop_gaia=False,
                                 epoch_1=1991.25, epoch_2=2016.0,
                                 apply_frame_rot=False,
                                 **kwargs):
    """Queries Gaia IDs, RA and Dec (back-propagated) for a simple cone search cross-match"""
    query = ""
    query = query + f"SELECT source_id, ra, dec, parallax, pmra, pmra_error, pmdec, pmdec_error, phot_g_mean_mag "
    if back_prop_gaia:
        query = query + ",array_element(a0, 1) as ra_prop, \
                          array_element(a0, 2) as dec_prop, \
                          array_element(a0, 3) as parallax_prop, \
                          array_element(a0, 4) as pmra_prop, \
                          array_element(a0, 5) as pmdec_prop, \
                          array_element(a0, 6) as rv_prop "

        query = query + f"FROM \
                         ( \
                         SELECT *, EPOCH_PROP(ra, dec, parallax, pmra, pmdec, radial_velocity, {epoch_2}, {epoch_1}) as a0  "

    query = query + "FROM gaiadr3.gaia_source AS gaia "

    if mag_lim is not None:
        query = query + f"WHERE phot_g_mean_mag < {mag_lim} "

    if back_prop_gaia:
        query = query + ') as p'

    print(query)

    job = Gaia.launch_job_async(query=query)

    return job.get_results()


def main():
    gaia_login()
    table = query_gaia_simple_conesearch(back_prop_gaia=True, mag_lim=5)
    return table


if __name__ == '__main__':
    main()
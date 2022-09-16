from astroquery.gaia import Gaia

username = 'rkievit'
Gaia.login(user=username)

gaia_table_name = 'user_' + username + '.gaia_sel12_prop'
xmatch_table_name = 'xm_hipp_gaia12'
search_radius = 1.0 #arcsecond

# astroquery function that automatically does circular (?) positional cross match between two catalogues
Gaia.cross_match(full_qualified_table_name_a=gaia_table_name,
                 full_qualified_table_name_b='external.hipparcos_newreduction',
                 results_table_name=xmatch_table_name,
                 verbose=True)


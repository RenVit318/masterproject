from astroquery.gaia import Gaia

username = 'rkievit'
Gaia.login(user=username)

gaia_table_name = 'user_' + username + '.gaia_sel12_prop' # Adjusted Gaia3 table containing all objects with mag<12 and backpropagated to J1991.25
#hipp_table_name = 'user_' + username + '.hipp_small' # Adjusted Hipparcos2 table selecting 150 random objects
hipp_table_name = 'public.hipparcos_newreduction'
xmatch_table_name = 'xm_hipp_gaia12_full'
search_radius = 1.0 #arcsecond

query = "SELECT hipp.hip, gaia.source_id, DISTANCE(\
         POINT(hipp.ra, hipp.dec),\
         POINT(gaia.ra_prop, gaia.dec_prop)) * 3600. AS dist_arcsec\
         FROM user_rkievit.gaia_sel12_prop AS gaia\
         JOIN public.hipparcos_newreduction AS hipp\
         ON 1 = CONTAINS(\
           POINT(hipp.ra, hipp.dec),\
           CIRCLE(gaia.ra_prop, gaia.dec_prop, 1. / 3600.) )"

print(query)


# Simple 1" Cross-Match
conesearch_xm_job = Gaia.launch_job_async(query=query,
                                          name='hipp2_gaia12_1as_conesearch')

print("Job Completed, Now Uploading...")

job = Gaia.upload_table_from_job(conesearch_xm_job)


# Mark the ra and dec columns as containing the sky coordinates
#Gaia.update_user_table(table_name=gaia_table_name,
#                       list_of_changes=[["ra","flags","Ra"], ["dec","flags","Dec"]])
#Gaia.update_user_table(table_name=hipp_table_name,
#                       list_of_changes=[["ra","flags","Ra"], ["dec","flags","Dec"]])





# astroquery function that automatically does circular (?) positional cross match between two catalogues
# automatically also writes the cross-match table to the archive
#Gaia.cross_match(full_qualified_table_name_a=hipp_table_name,
#                 full_qualified_table_name_b=gaia_table_name,
#                 results_table_name=xmatch_table_name,
#                verbose=True)


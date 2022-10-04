import copy

def extract_sky_positions(tab):
    """Extract RA and Dec arrays from labeled tables. Note, names must be ra and dec"""
    ra = tab['ra']
    de = tab['dec']
    return ra, de


def convert_to_ids(xm_table, tab_g, tab_h):
    """Takes an array with [hip_array_number, gaia_array_number, distance]
    and converts it to     [hip.hip, gaia.source_id, distance]"""
    tab = copy.deepcopy(xm_table)
    for i in range(xm_table.shape[0]):
        tab[i][0] = tab_h[1].data['hip'][int(xm_table[i][0])]
        tab[i][1] = tab_g[1].data['source_id'][int(xm_table[i][1])]

    return tab


def main():
    pass


if __name__ == '__main__':
    main()
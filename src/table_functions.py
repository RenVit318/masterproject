import copy
import numpy as np 

def extract_sky_positions(tab):
    """Extract RA and Dec arrays from labeled tables. Note, names must currently be ra and dec"""
    ra = tab['ra']
    de = tab['dec']
    return ra, de


def extract_proper_motions(tab):
    """Extract pmra and pmde arrays from labeled tables. Note, names currently must be pmra and pmde"""
    pmra = tab['pmra']
    pmde = tab['pmde']
    return pmra, pmde


def convert_to_ids(xm_table, tab_g, tab_h):
    """Takes an array with [hip_array_number, gaia_array_number, distance]
    and converts it to     [hip.hip, gaia.source_id, distance]"""
    tab = copy.deepcopy(xm_table)
    for i in range(xm_table.shape[0]):
        tab[i][0] = tab_h[1].data['hip'][int(xm_table[i][0])]
        tab[i][1] = tab_g[1].data['source_id'][int(xm_table[i][1])]

    return tab


def batch_table(table, num_batches=None, batch_size=None):
    """Split a table up in to either a number of chunks, or a max. chunk size"""
    idxs = np.arange(len(tab))
    if num_batches is not None and batch_size is not None:
        raise NotImplementedError("Please provide either only num_batches or batch_size")
    elif num_batches is not None:
        batches_idxs = np.split(idxs, num_batches)
    elif batch_size is not None
        num_batches = np.ceil(table.shape[0]/batch_size)
        batches_idxs = np.split(idxs, num_batches)
    else:
        raise ValueError("Please provide either num_batches or batch_size")

    batches = []
    for idxs in batches_idxs:
        batches.append(table[idxs])    

    return batches
        

def combine_tables(batches):
    """Combine a list of tables back into one full table"""
    pass


def main():
    pass


if __name__ == '__main__':
    main()

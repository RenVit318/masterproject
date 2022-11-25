############################
#
# All functions related to manipulation and extraction of tables
# Most code made to work with type: astropy.table.Table 
#
############################

import copy
import numpy as np


def extract_sky_positions(tab, ra_id='ra', de_id='dec'):
    """Extract RA and Dec arrays from labeled tables. Note, names must currently be ra and dec"""
    ra = np.array(tab[ra_id], dtype=np.float64)  # casting into float64 in a numpy array shouldn't break anything
    de = np.array(tab[de_id], dtype=np.float64)  # but does fix issues with numba
    return ra, de


def extract_proper_motions(tab, pmra_id='pmra', pmde_id='pmdec'):
    """Extract pmra and pmde arrays from labeled tables. Note, names currently must be pmra and pmdec"""
    pmra = np.array(tab[pmra_id], dtype=np.float64)  # same comment as in extract_sky_positions
    pmde = np.array(tab[pmde_id], dtype=np.float64)
    return pmra, pmde


def convert_to_ids(xm_table, tab_g, tab_h):
    """Takes an array with [hip_array_number, gaia_array_number, distance]
    and converts it to     [hip.hip, gaia.source_id] and [distance]"""
    tab = np.zeros((xm_table.shape[0], 3), dtype=np.int64)
    print(tab.shape)

    for i in range(xm_table.shape[0]):
        tab[i][0] = int(tab_h['hip'][int(xm_table[i][0])])  # tab_h[1].data['hip'] <- old. Don't use because we only extract tab[1].data
        tab[i][1] = np.int64(tab_g['source_id'][int(xm_table[i][1])])  # in read gaia_hipp_data. This should also work better with queried gaia data
        tab[i][2] = int(xm_table[i][1]) # Save this for proper comparison to its own GaiaBaseCat

    return tab, xm_table[:, 2]


def batch_table(table, num_batches=None, batch_size=None):
    """Split a table up in to either a number of chunks, or a max. chunk size"""
    idxs = np.arange(len(table))
    if num_batches is not None and batch_size is not None:
        raise NotImplementedError("Please provide either only num_batches or batch_size")
    elif num_batches is not None:
        batches_idxs = np.array_split(idxs, int(num_batches))
    elif batch_size is not None:
        num_batches = int(np.ceil(len(table)/batch_size))
        batches_idxs = np.array_split(idxs, num_batches)
    else:
        raise ValueError("Please provide either num_batches or batch_size")

    batches = []
    for idxs in batches_idxs:
        batches.append(table[idxs])   # Not the most efficient, but works relatively quickly
    return batches
        

def read_tables(table_path, multiple=False):
    """Either read in a single table, or read all tables following a 'ls'-like search, with * and ?"""
    from astropy.table import Table
    import glob

    if not multiple:
        return Table.read(table_path)
    elif multiple:
        all_table_paths = glob.glob(table_path)
        all_tables = []
        for t_path in all_table_paths:
            t = Table.read(t_path)
            all_tables.append(t)
        return all_tables


def main():
    pass


if __name__ == '__main__':
    main()

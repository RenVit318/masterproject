############################
#
# Functions for the second step in the crossmatch pipeline
# Multiple methods of describing the best neighbour for a
# Hipparcos star with multiple possible matches
#
############################


def nearest():
    """"""
    pass


def likeliest_position():
    pass


def brightest():
    pass


def select_best_neighbour(tab_xm, best_match_selection='nearest'):
    """From a table containing all possible cross-matches, selects the best neighbour
       for each hipparcos object with >1 reported neighbour"""

    # Check all objects with >1 candidate
    # https://stackoverflow.com/questions/30003068/how-to-get-a-list-of-all-indices-of-repeated-elements-in-a-numpy-array
    hip_sort_idxs = np.argsort(tab_xm[:, 0])  # I think, check this shape

    # Check which selection type
    if best_match_selection == 'nearest':
        select_func = nearest
    elif best_match_selection == 'brightest':
        select_func = brightest
    elif best_match_selection == 'likeliest_position':
        select_func = likeliest_position

    # Call selection function


def main():
    r_path = '../results/'
    xm_tab_name = ''
    tab_xm = np.load(r_path + xm_tab_name)
    best_match_selection = 'nearest'
    select_best_neigbhour()


if __name__ == '__main__':
    main()

import numpy as np

def sind(x):
    return np.sin(np.radians(x))
def cosd(x):
    return np.cos(np.radians(x))

def angular_distance(ra_a, dec_a, ra_b, dec_b, nu=None, unit='radians', nu_thresh=1e-3):
    """Calculates the angular distance between point A and point B without any assumptions
    For the equation, and its derivation see https://en.wikipedia.org/wiki/Angular_distance#General_case
    IMPORTANT: ra and dec are expected to be in RADIANS, if they're not use unit to specify"""

    if unit == 'radians':
        pass
    elif unit == 'degrees':
        # Convert degrees to radians
        [ra_a, dec_a, ra_b, dec_b] = np.radians([ra_a, dec_a, ra_b, dec_b])
    else:
        raise ValueError(f"Unit {unit} is not recognized")

    theta = np.zeros(ra_a.shape[0])

    for i in range(theta.shape[0]):
        if nu is not None:
            if nu[i] > nu_thresh: # 1e-3 is chosen based on where the errors appear to start 'kicking' in
                # Full equation of displacement
                sin_factor = np.sin(dec_a[i])*np.sin(dec_b[i])
                cos_factor = np.cos(dec_a[i])*np.cos(dec_b[i])*np.cos(ra_a[i]-ra_b[i])
                theta[i] = np.arccos(sin_factor + cos_factor)
                continue
        else:
            # Small angular distance approximation to account for floating point errors
            factor1 = (ra_a[i] - ra_b[i]) * np.cos(dec_a[i])
            factor2 = dec_a[i] - dec_b[i]
            theta[i] = np.sqrt( (factor1 ** 2) + (factor2 ** 2))

    return theta

def predict_G_bv(hp_mag, b_v):
    """Gaia magnitude prediction based on B-V colour from the Gaia handbook"""
    return -0.02392 - 0.4069 * b_v + 0.04569 * (b_v**2) - 0.0452 * (b_v**3) + hp_mag, 0.02417

def predict_G_vi(hp_mag, v_i):
    """Gaia magnitude prediction based on V-I colour from the Gaia handbook"""
    return +0.01546 - 0.4308 * v_i - 0.01872 * (v_i**2) + hp_mag, 0.08181

def rayleigh(x, sigma):
    return (x/(sigma**2)) * np.exp(-0.5*((x/sigma)**2.))


def compute_error_normalized_distance(pos1, pos2, unc1, unc2, method='full', coordinates='spherical', unit='mas'):
    """Calculates the 2-dimensional error normalized distance, D, between two objects with uncertainties,
    following a variety of possible methods as listed below.

        None:        No error normalization performed, this function just returns the distance between x,y
        Simple:      Uncertainties are approximated as a circle with radius sigma^2 = 0.5(sigma_x^2 + sigma_y^2)
        Directional: Computes D separately in the x and y directions first, and then combines using Pythagoras
        Full:        Utilizes the full covariance matrix and distance vector

    Inputs:
        pos1: Nx2 array containing the [x,y] coordinates of the first object
        pos2: See above, but for the second object
        unc1: Nx3 array containing all (co)variance terms of object 1 necessary for the chosen methods.
              List as [sigma_x^2, sigma_y^2, sigma_xy], not given values are treated as 0. # Maybe this should be np.nan?
        unc2: See above, but for the second object
        method: Choose from none, simple, directional, full. For descriptions see above
        coordinates: If set to spherical, transforms x |-> x cos(y) for distance calculation to account for spherical effects near poles
        unit: Check to ensure trigonometry is performed properly in the case of spherical coordinates

    Outputs:
        error_norm_distances: array of size N describing the error normalized distance for all provided object pairs
    """

    # Bad inputs check
    if not method in ['none', 'simple', 'directional', 'full']:
        raise ValueError(f'method {method} not implemented in this function')
    if not ((pos1.shape == pos2.shape) and (unc1.shape == unc2.shape)):
        raise ValueError(f'position or uncertainty arrays are of different shapes')

    if len(pos1.shape) == 1: # Then N=1
        x1, y1 = pos1
        sigma_x1, sigma_y1, sigma_xy1 = unc1
        x2, y2 = pos2
        sigma_x2, sigma_y2, sigma_xy2 = unc2
    elif len(pos1.shape) == 2:
        # Extract all columns, and store in separate variables
        x1, y1 = pos1[:,0], pos1[:,1]
        sigma_x1, sigma_y1, sigma_xy1 = unc1[:,0], unc1[:,1], unc1[:,2]
        x2, y2 = pos2[:,0], pos2[:,1]
        sigma_x2, sigma_y2, sigma_xy2 = unc2[:,0], unc2[:,1], unc2[:,2]
    else:
        raise ValueError(f'Array with dimension {len(pos1.shape)} not supported')

    # Transform x to x cos(y) to account for spherical effects
    if coordinates == 'spherical':
        if unit == 'rad':
            y1_rad, y2_rad = y1, y2
        elif unit == 'deg':
            y1_rad, y2_rad = np.radians(y1), np.radians(y2)
        elif unit == 'mas':
            y1_rad, y2_rad = np.radians(y1/3.6e6), np.radians(y2/3.6e6) # 3600 * 1e3 = 3.6e6
        else:
            raise NotImplementedError(f'unit {unit} not recognized')
        x1_corr = x1 * np.cos(y1_rad)
        x2_corr = x2 * np.cos(y2_rad)
    else:
        x1_corr = x1
        x2_corr = x2

    # Only method that computes D not from the combined distance
    if method == 'directional':
        total_x_unc = np.sqrt(sigma_x1 + sigma_x2)
        total_y_unc = np.sqrt(sigma_y1 + sigma_y2)

        x_dist_norm = (x1_corr - x2_corr)/total_x_unc
        y_dist_norm = (y1 - y2)/total_y_unc

        return np.sqrt(x_dist_norm**2. + y_dist_norm**2.)

    elif method == 'full':
        # chi^2 = v^T S^-1 v
        # Setup distance vector and covariance matrix
        if len(pos1.shape) > 1:
            delta_distance = np.stack(((x1_corr-x2_corr), (y1-y2)), axis=1)

            cov_mat_top = np.stack((sigma_x1 + sigma_x2, sigma_xy1 * np.sqrt(sigma_x1*sigma_y1) + sigma_xy2 * np.sqrt(sigma_x2*sigma_y2)), axis=1)
            cov_mat_bot = np.stack((sigma_xy1 * np.sqrt(sigma_x1*sigma_y1) + sigma_xy2 * np.sqrt(sigma_x2*sigma_y2), sigma_y1 + sigma_y2), axis=1)

            cov_matrix = np.stack((cov_mat_top, cov_mat_bot), axis=1)
        else:
            delta_distance = np.array([x1_corr - x2_corr, y1 - y2])
            #print(f'v:{delta_distance}')
            cov_matrix = np.array([[sigma_x1+sigma_x2, sigma_xy1*np.sqrt(sigma_x1*sigma_y1)+sigma_xy2*np.sqrt(sigma_x2*sigma_y2)], [sigma_xy1*np.sqrt(sigma_x1*sigma_y1)+sigma_xy2*np.sqrt(sigma_x2*sigma_y2), sigma_y1+sigma_y2]])
            #print(f'S: {cov_matrix}')
            #cov_matrix = np.array([[sigma_x1+sigma_x2, sigma_xy1+sigma_xy2], [sigma_xy1+sigma_xy2, sigma_y1+sigma_y2]])
        cov_matrix_inverse = np.linalg.inv(cov_matrix) # Inverts each of the N matrices individually
        #print(f'S^-1: {cov_matrix_inverse}')

        # Calculate chi squared
        if len(pos1.shape) > 1:
            # https://stackoverflow.com/questions/35894631/multiply-array-of-vectors-with-array-of-matrices-return-array-of-vectors
            # Figure out how einsum exactly works (and better name?)
            vt_Sinverse = np.einsum('ij,ijk->ik', delta_distance, cov_matrix_inverse)
            chi_squared = np.einsum('ij,ij->i', vt_Sinverse, delta_distance)
        else:
            vt_Sinverse = np.matmul(delta_distance, cov_matrix_inverse)
            chi_squared = np.matmul(vt_Sinverse, delta_distance)
            #print(f'chi_squared {chi_squared}')

        return np.sqrt(np.abs(chi_squared))

    # Start by computing distances
    distances = np.sqrt((x1_corr - x2_corr)**2. + (y1 - y2)**2.)

    if method == 'none':
        return distances
    elif method == 'simple':
        avg_unc1_sq = 0.5*(sigma_x1 + sigma_y1)
        avg_unc2_sq = 0.5*(sigma_x2 + sigma_y2)
        combined_unc = np.sqrt(avg_unc1_sq + avg_unc2_sq)

        return distances/combined_unc







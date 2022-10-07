import numpy as np

def sind(x):
    return np.sin(np.radians(x))
def cosd(x):
    return np.cos(np.radians(x))

def angular_distance(ra_a, dec_a, ra_b, dec_b, nu, unit='radians', nu_thresh=1e-3):
    """Calculates the angular distance between point A and point B without any assumptions
    For the equation, and its derivation see https://en.wikipedia.org/wiki/Angular_distance#General_case
    IMPORTANT: ra and dec are expected to be in RADIANS, if they're not use unit to specify"""

    if unit =='radians':
        pass
    elif unit=='degrees':
        # Convert degrees to radians
        [ra_a, dec_a, ra_b, dec_b] = np.radians([ra_a, dec_a, ra_b, dec_b])
    else:
        raise ValueError(f"Unit {unit} is not recognized")


    theta = np.zeros(ra_a.shape[0])

    for i in range(theta.shape[0]):
        if nu[i] > nu_thresh: # 1e-3 is chosen based on where the errors appear to start 'kicking' in
            # Full equation of displacement
            sin_factor = np.sin(dec_a[i])*np.sin(dec_b[i])
            cos_factor = np.cos(dec_a[i])*np.cos(dec_b[i])*np.cos(ra_a[i]-ra_b[i])
            theta[i] = np.arccos(sin_factor + cos_factor)
        else:
            # Small angular distance approximation to account for floating point errors
            factor1 = (ra_a[i] - ra_b[i]) * np.cos(dec_a[i])
            factor2 = dec_a[i] - dec_b[i]
            theta[i] = np.sqrt( (factor1 ** 2) + (factor2 ** 2))

    return theta






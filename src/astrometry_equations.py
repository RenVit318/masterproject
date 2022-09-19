import numpy as np


def angular_distance(ra_a, dec_a, ra_b, dec_b, unit='radians'):
    """Calculates the angular distance between point A and point B without any assumptions
    For the equation, and its derivation see https://en.wikipedia.org/wiki/Angular_distance#General_case
    IMPORTANT: ra and dec are expected to be in RADIANS, if they're not use unit to specify"""

    print(np.min(ra_a))

    if unit =='radians':
        pass
    elif unit=='degrees':
        # Convert degrees to radians
        [ra_a, dec_a, ra_b, dec_b] = np.radians([ra_a, dec_a, ra_b, dec_b])
    else:
        raise ValueError(f"Unit {unit} is not recognized")


    sin_factor = np.sin(dec_a)*np.sin(dec_b)
    print(np.min(sin_factor))
    cos_factor = np.cos(dec_a)*np.cos(dec_b)*np.cos(ra_a-ra_b)
    print(np.min(cos_factor))
    theta = np.arccos(sin_factor + cos_factor)
    print(np.min(theta))

    return theta


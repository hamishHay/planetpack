
def test():
    print("HELLO WORLD")

    return 1

def ecc_ocean_diss(body, e, h, cd):
    """
    Function implements the scaling laws from Chen et al. (2014) for the
    eccentricity tidal potential, given an eccentricity e, ocean thickness h,
    and drag coefficient cd. body input must be a planetPy type body.
    """

    import planetPy.body as p_body
    import numpy as np

    # Check that body is an instance of class planetPy.body
    if type(body) != p_body.Moon and type(body) != p_body.Planet:
        raise ValueError("Argument 1 must have a type of", p_body.Moon, "or", p_body.Planet)

    # Get constants
    den = 1000.0
    rot = body.mean_motion
    R = body.radius
    g = body.gravity

    # Calculate v
    v = 0.13 * cd * rot**3.0 * R**5.0 * e / (g * h**2.0)

    # Calculate E_o
    E0 = den * rot**6.0 * R**8.0 / (g**2.0 * h)

    # Calculate Edot_o
    Edot0 = E0 * v / R**2.0

    # Calculate dissipated energy [W]
    Edot = 519.0 / 35.0 * np.pi * Edot0 * e**2.0

    return Edot

def obliq_ocean_diss(body, obl, h, cd):
    """
    Function implements the scaling laws from Chen et al. (2014) for the
    obliquity tidal potential, given an eccentricity e, ocean thickness h,
    and drag coefficient cd. body input must be a planetPy type body.
    """

    import planetPy.body as p_body
    import numpy as np

    # Check that body is an instance of class planetPy.body
    if type(body) != p_body.Moon and type(body) != p_body.Planet:
        raise ValueError("Argument 1 must have a type of", p_body.Moon, "or", p_body.Planet)

    # Get constants
    den = 1000.0
    rot = body.mean_motion
    R = body.radius
    g = body.gravity

    # Calculate v
    v = np.sqrt(1. + ((200./3. * 0.4 * cd * g * obl)/(rot**2.0 * R))**2.0)
    v = np.sqrt(-1. + v)
    v *= rot**3.0 * R**4.0 / (20 * np.sqrt(2.) * g * h)

    # Calculate u0
    u0 = rot**3.0 * R**3.0 / (g * h)

    # Calculate dissipated energy [W]
    Edot = 12. * np.pi * den * h * rot**2.0 * R**2.0 * v* obl**2.0
    Edot /= (1. + (400 * v**2.0)/(u0**2.0 * R**2.0))

    return Edot

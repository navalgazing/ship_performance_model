import numpy as np


def calc_ice_resistance(ice_conc, ice_thick, ship_speed_ms, lpp, beam, draft, lpar, lbow, stem_angle, ice_density=None, g=None):
    # Set default values for gravity and ice density
    if ice_density is None:
        ice_density = 900

    if g is None:
        g = 9.805

    # Set coefficients for Finnish-Swedish rules
    f_1 = 0.23
    f_2 = 4.58
    f_3 = 1.47
    f_4 = 0.29
    g_1 = 18.9
    g_2 = 0.67
    g_3 = 1.55

    fn = ship_speed_ms / np.sqrt(g * lpp)

    local_ice_thickness = ice_conc * ice_thick

    if local_ice_thickness < 0.3:
        # Luofeng Huang's method
        # TODO this is a coefficient specific to COSCO, change if this is used on other ships
        cr = 0.20

        ice_resistance = cr * ice_density * beam * ice_thick * ship_speed_ms ** 2 * ice_conc ** 1.5 * fn ** (-1)

    else:
        # Finnish-Swedish Ice Class method
        c_1 = f_1 * ((beam * lpar * ice_thick) / (((2 * draft) / beam) + 1)) + (1 + 0.021 * stem_angle) * (f_2 * beam * ice_thick ** 2 + f_3 * lbow * ice_thick ** 2 + f_4 * beam * lbow * ice_thick)
        c_2 = (1 + 0.063 * stem_angle) * (g_1 * ice_thick ** 1.5 + g_2 * beam * ice_thick) + g_3 * ice_thick * (1 + ((1.2 * draft) / beam)) * (beam ** 2 / np.sqrt(lpp))

        ice_resistance = c_1 + (c_2 * ship_speed_ms) * 1000

    return ice_resistance

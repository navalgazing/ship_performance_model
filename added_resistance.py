import numpy as np
import math
import scipy as sp


# ITTC added resistance

# Due to the effects of wind (ITTC method)
# TODO add units to all documentation / check those already there
def wind_resitance(a_xv, c_x, v_wr, rho_a):
    """
    Caculates the wind resistance using the method from ITTC Recommended procedure and guidelines - Analysis of Speed/Power Trial Data (7.5-04-01-01.2 section 4.3)

    :param a_xv: area of maximum transverse section area (m^2)
    :param c_x: wind resistance coefficient at relative wind direction
    :param v_wr: relative wind speed (m/s) TODO Check units
    :param rho_a: the mass density of air (kg/m^3)
    :return r_aa: resistance increase due to relative wind ()
    """

    r_aa = 0.5 * rho_a * v_wr ** 2 * c_x * a_xv

    # TODO is there an approximation for c_x?
    return r_aa


###########################################################################################
#### Fujiwara method for added resistance in waves ####
###########################################################################################

def wind_resistance_fujiwara(psi_wr, a_yv, a_xv, a_od, loa, b, c_mc, h_br, h_c):
    """
    A general regression formula based on model tests in wind tunnels for various ships developed by Fujiwara et al. Fujiwara, T., Ueno, M. and Ikeda, Y.: "A
    New Estimation Method of Wind Forces and Moments acting on Ships on the basis of Physical Component Models", J. JASNAOE, Vol.2, 2005.

    :param psi_wr: Relative wind direction; 0 means heading winds.
    :param a_yv: Projected lateral area above the waterline.
    :param a_xv: Area of maximum transverse section exposed to the winds.
    :param a_od: Lateral projected area of superstructures etc. on deck.
    :param loa: Length overall.
    :param b: Ship breadth.
    :param c_mc: Horizontal distance from midship section to centre of lateral projected area a_yv.
    :param h_br: Height of top of superstructure (bridge etc.).
    :param h_c: Height from waterline to centre of lateral projected area a_yv.
    :return c_aa: Wind resistance coefficient.
    """
    # TODO check coefficient indexing is correct
    # Converting degrees to radians
    psi_wr = abs(psi_wr)
    psi_wr = math.radians(psi_wr)
    if math.radians(0) <= psi_wr < math.radians(90):
        c_aa = wind_less_than_90(psi_wr, a_yv, a_xv, a_od, loa, b, c_mc, h_br)

    elif psi_wr == math.radians(90):
        c_aa = wind_equal_to_90(a_yv, a_xv, a_od, loa, b, c_mc, h_br, h_c)

    elif math.radians(90) < psi_wr <= math.radians(180):
        c_aa = wind_greater_than_90(psi_wr, a_yv, a_xv, a_od, loa, b, h_br, h_c)

    return c_aa


def wind_equal_to_90(a_yv, a_xv, a_od, loa, b, c_mc, h_br, h_c):
    """
    Calculates the wind resistance coefficient when relative wind direction is 90 degrees speed.e. beam winds.

    :param a_yv: Projected lateral area above the waterline.
    :param a_xv: Area of maximum transverse section exposed to the winds.
    :param a_od: Lateral projected area of superstructures etc. on deck.
    :param loa: Length overall.
    :param b: Ship breadth.
    :param c_mc: Horizontal distance from midship section to centre of lateral projected area a_yv.
    :param h_br: Height of top of superstructure (bridge etc.).
    :param h_c: Height from waterline to centre of lateral projected area a_yv.
    :return: Wind resistance coefficient.
    """
    c_aa = 0.5 * (wind_greater_than_90(100, a_yv, a_xv, a_od, loa, b, h_br, h_c) + wind_less_than_90(80, a_yv, a_xv, a_od, loa, b, c_mc, h_br))
    return c_aa


def wind_greater_than_90(psi_wr, a_yv, a_xv, a_od, loa, b, h_br, h_c):
    """
    Calculates the wind resistance coefficient when relative wind direction is greater than 90 degrees speed.e. between beam and following winds.

    :param psi_wr: Relative wind direction; 0 means heading winds.
    :param a_yv: Projected lateral area above the waterline.
    :param a_xv: Area of maximum transverse section exposed to the winds.
    :param a_od: Lateral projected area of superstructures etc. on deck.
    :param loa: Length overall.
    :param b: Ship breadth.
    :param h_br: Height of top of superstructure (bridge etc.).
    :param h_c: Height from waterline to centre of lateral projected area a_yv.
    :return: Wind resistance coefficient.
    """
    c_xli = 1.901 + -12.727 * (a_yv / (loa * h_br)) + -24.407 * (a_xv / a_yv) + 40.310 * (b / loa) + 5.481 * (a_xv / (b * h_br))
    c_lf = -0.018 + 5.091 * (b / loa) + -10.367 * (h_c / loa) + 3.011 * (a_od / loa ** 2) + 0.341 * (a_xv / b ** 2)
    c_alf = 0.314 + 1.117 * (a_od / a_yv)
    c_aa = c_aa_calc(c_lf, c_xli, c_alf, psi_wr)
    return c_aa


def wind_less_than_90(psi_wr, a_yv, a_xv, a_od, loa, b, c_mc, h_br):
    """
    Calculates the wind resistance coefficient when relative wind direction is less than 90 degrees speed.e. between beam and heading winds.

    :param psi_wr: Relative wind direction; 0 means heading winds.
    :param a_yv: Projected lateral area above the waterline.
    :param a_xv: Area of maximum transverse section exposed to the winds.
    :param a_od: Lateral projected area of superstructures etc. on deck.
    :param c_mc: Horizontal distance from midship section to centre of lateral projected area a_yv.
    :param loa: Length overall.
    :param b: Ship breadth.
    :param h_br: Height of top of superstructure (bridge etc.).
    :return c_aa: Wind resistance coefficient.
    """
    c_lf = 0.922 - 0.507 * (a_yv / (loa * b)) - 1.162 * (c_mc / loa)
    c_xli = -0.458 - 3.245 * (a_yv / (loa * h_br)) + 2.313 * (a_xv / (b * h_br))
    c_alf = 0.585 + 0.906 * (a_od / a_yv) - 3.239 * (b / loa)
    c_aa = c_aa_calc(c_lf, c_xli, c_alf, psi_wr)
    return c_aa


def c_aa_calc(c_lf, c_xli, c_alf, psi_wr):
    """

    :param c_lf: Coefficient of longitudinal flow drag. TODO is this correct
    :param c_xli: Coefficient of induced drag. TODO is this correct?
    :param c_alf: TODO what is this?
    :param psi_wr: Relative wind direction; 0 means heading winds.
    :return c_aa: Wind resistance coefficient.
    """
    c_aa = c_lf * np.cos(psi_wr) + c_xli * (np.sin(psi_wr) - 0.5 * np.sin(psi_wr) * np.cos(psi_wr) ** 2) * np.sin(psi_wr) * np.cos(psi_wr) + c_alf * np.sin(psi_wr) * np.cos(psi_wr) ** 3
    return c_aa


##################################################################################
#### STAWave-1 added wave resistance calculation ####
##################################################################################

def added_wave_resistance_stawave1(rho, h, b, l_bwl, g=9.805):
    """

    :param rho: Density of water (kg/m**3)
    :param h: Significant wave height (m)
    :param b: Beam (m)
    :param l_bwl: Length of the bow on the waterline to 95% of maximum beam (m)
    :param g: Acceleration due to gravity (m/s**2)
    :return: Added resistance due to waves (kN) TODO units? Probably kN?
    """
    r_awl = 1 / 16 * rho * g * h ** 2 * b * np.sqrt(b / l_bwl)
    return r_awl


#################################################################################


def added_wave_resistance_liu_et_al_2016(rho, beam, c_b, l_pp, l_e, k_yy, f_n, zeta_a, omega, _lambda, g=9.805):
    """
    l_e: distance from the forward perpendicular to the position of 99% maximum beam (m)
    :return: Added resistance due to waves (kN)
    """

    if f_n == 0:
        return 0

    if math.isnan(f_n):
        return

    E = math.atan(beam / (2.0 * l_e))

    r_awr = (2.25 / 2) * rho * g * beam * zeta_a ** 2 * np.sin(E) ** 2 * (1 + 5 * np.sqrt(l_pp / _lambda) * f_n) * (0.87 / c_b) ** (1.0 + 4.0 * np.sqrt(f_n))

    a1 = 60.3 * (c_b ** 1.34) * (0.87 / c_b) ** (1 + f_n)

    if f_n < 0.12:
        a2 = 0.0072 + 0.1676 * f_n
    elif f_n >= 0.12:
        a2 = (f_n ** 1.5) * np.exp(-3.5 * f_n)

    if f_n < 0.05:
        omega_bar = ((np.sqrt(l_pp / g) * np.cbrt(k_yy) * (0.05 ** 0.143)) / 1.17) * omega
    elif f_n >= 0.05:
        omega_bar = ((np.sqrt(l_pp / g) * np.cbrt(k_yy / l_pp) * (f_n ** 0.143)) / 1.17) * omega

    if c_b < 0.75:
        if omega_bar < 1:
            b1 = 11.0
            d1 = 14.0
        else:
            b1 = -8.5
            d1 = -566.0 * ((l_pp / beam) ** -2.66) * 6.0
    elif c_b >= 0.75:
        if omega_bar < 1:
            b1 = 11.0
            d1 = 566.0 * ((l_pp / beam) ** -2.66)
        else:
            b1 = -8.5
            d1 = -566.0 * ((l_pp / beam) ** -2.66) * 6.0

    r_awm = 4.0 * rho * g * zeta_a ** 2 * beam ** 2 / l_pp * omega_bar ** b1 * np.exp(b1 / d1 * (1 - omega_bar ** d1)) * a1 * a2

    r_aw = r_awm + r_awr

    return r_aw


################################################################################

def added_wave_resistance_sta2(rho, g, zeta_a, beam, draft, l_pp, f_n, k_yy, omega, c_b, v_s_metrespersecond, _lambda):
    v_s_kts = v_s_metrespersecond * 1.94384
    f_n = v_s_metrespersecond / np.sqrt(l_pp * g)

    omega_bar = (np.sqrt(l_pp / g) * np.cbrt(k_yy/l_pp)) * f_n**0.143 * omega / 1.17  # / (1.17 * f_n ** (-0.143))) * omega

    a1 = 60.3 * (c_b ** 1.34)

    a2 = f_n**1.5 * np.exp(-3.5 * f_n)

    if omega_bar < 1:
        b1 = 11.0
    else:
        b1 = -8.5

    if omega_bar < 1:
        d1 = 14.0
    else:
        d1 = -566 * (l_pp / beam) ** -2.66

    raw_bar_omega = (omega_bar ** b1) * np.exp((b1/d1) * (1 - (omega_bar ** d1))) * a1 * (f_n ** 1.5) * np.exp(-3.5 * f_n)

    r_awm = 4 * rho * g * (zeta_a ** 2) * (beam ** 2) / l_pp * a1 *a2 * omega_bar**b1 * np.exp((b1/d1)*(1-omega_bar**d1))  # * raw_bar_omega

    k = (2 * np.pi) / _lambda

    f1 = 0.692 * (v_s_kts/np.sqrt(draft * g)) ** 0.769 + 1.81 * c_b ** 6.95

    alpha1omega = ((np.pi**2 * sp.special.i1(1.5 * k * draft)**2) / (np.pi ** 2 * sp.special.i1(1.5 * k * draft)**2 + sp.special.k1(1.5 * k * draft)**2)) * f1

    r_awr = 0.5 * rho * g * zeta_a**2 * beam * alpha1omega

    r_aw = r_awr + r_awm

    return r_aw

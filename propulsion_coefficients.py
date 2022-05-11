import numpy as np
import sys


def calc_relative_rotative_efficiency(num_of_props, bar, c_p, lcb, pitch_diameter_ratio=None):
    """
    This estimates the relative rotative efficiency using Holtrop's method. See "Ship Resistance and Propulsion - Molland, Turnock & Hudson" Section 16.3.2.1 and 16.3.2.2

    :param num_of_props: Number of propellers. This can be either 1 or 2.
    :param bar: Blade area ratio = expanded blade area / propeller disc area.
    :param c_p: Prismatic coefficient.
    :param lcb: Longitudinal centre of buoyancy. This is the position of the centre of buoyancy forward of 0.5 LWL as a percentage of LWL.
    :param pitch_diameter_ratio: Propeller pitch / propeller diameter.
    :return: Relative rotative efficiency.
    """

    if num_of_props == 1:
        eta_r = 0.9922 - 0.05908 * bar + 0.07424 * (c_p - 0.0225 * lcb)  # From Holtrop, for single screw ONLY TODO need to check where LCB is measured from
    elif num_of_props == 2:
        if pitch_diameter_ratio is None:
            sys.exit("Pitch diameter ratio is required when the number of propellers = 2")
        eta_r = 0.9737 + 0.111 * (c_p - 0.0225 * lcb) - 0.06325 * pitch_diameter_ratio  # from Holtrop, for twin screw ONLY TODO add if statement for these
    else:
        sys.exit("The number of propellers must be either 1 or 2")
    return eta_r


def calc_wake_fraction(c_b):
    """
    This estimates the wake fraction using Taylor's method. See "Ship Resistance and Propulsion - Molland, Turnock & Hudson" Section 8.8.2.1

    :param c_b: Block coefficient
    :return: Wake fraction
    """
    # simplest version used speed.e. Taylor TODO implement a better version
    w_t = 0.5 * c_b - 0.05
    return w_t


def calc_thrust_deduction_factor(b, l, t, d, c_p, lcb, c_stern):
    """
    This estimates the thrust deduction factor using Holtrop's method. See "Ship Resistance and Propulsion - Molland, Turnock & Hudson" Section 8.8.2.2
    :param b: Moulded breadth (m)
    :param l: Waterline length (m)
    :param t: Moulded draft (m)
    :param d: Propeller diameter (m)
    :param c_p: Prismatic coefficient
    :param lcb: Longitudinal centre of buoyancy. This is the position of the centre of buoyancy forward of 0.5 LWL as a percentage of LWL.
    :param c_stern: Coefficient depending on stern type. C_stern values depending on after body form are as follows....Pram with gondola: c_stern = -25. V-shaped sections: c_stern = -10. Normal section shape: c_stern = 0. U-shaped sections with Hogner stern: c_stern = 10.
    :return: Thrust deduction factor.
    """
    # TODO this is valid for single screw only. Add another method for twin screw.
    t = 0.25014 * (b / l) ** 0.28956 * np.sqrt(b * t / d) ** 0.2624 / (1 - c_p + 0.0225 * lcb) ** 0.01762 + 0.0015 * c_stern
    return t


def calc_hull_efficiency(t, w_t):
    """
    Calculate the hull efficiency from the wake fraction and thrust deduction factor.

    :param t: Thrust deduction factor.
    :param w_t: Wake fraction.
    :return: Hull efficiency.
    """
    eta_h = (1 - t) / (1 - w_t)
    return eta_h


def calc_kt_kq_j(pitch_diameter_ratio, bar, blade_num, prop_diameter, revs_per_second, v_a):
    """
    Calculates the open water efficiency based on using a Wageningen series propeller.See "Ship Resistance and Propulsion - Molland, Turnock & Hudson" Section 16.2.1.1.

    :param pitch_diameter_ratio: Propeller pitch / propeller diameter.
    :param bar: Blade area ratio (expanded blade area / propeller disc area).
    :param blade_num: Number of propeller blades.
    :param prop_diameter: Propeller diameter (m).
    :param revs_per_second: Propeller rate of revolution (rps).
    :param v_a: Wake speed(ship speed * [1 - w_t]) (m/s).
    :return: Open water efficiency.
    """

    # Import polynomial coefficients for Kt and Kq calculations
    wageningen_polynomial_coefficients = np.genfromtxt('wageningen_polynomial_coefficients.csv', delimiter=',')

    # Calculating J based on propeller speed (RPS) and diameter
    j = v_a / (revs_per_second * prop_diameter)

    # Calculating Kq and Kt from the polynomial coefficients
    k_t = np.sum(
        wageningen_polynomial_coefficients[0:39, 0] * (j ** wageningen_polynomial_coefficients[0:39, 1]) * (pitch_diameter_ratio ** wageningen_polynomial_coefficients[0:39, 2]) * (
                bar ** wageningen_polynomial_coefficients[0:39, 3]) * (
                blade_num ** wageningen_polynomial_coefficients[0:39, 4]))
    k_q = np.sum(
        wageningen_polynomial_coefficients[0:47, 5] * (j ** wageningen_polynomial_coefficients[0:47, 6]) * (pitch_diameter_ratio ** wageningen_polynomial_coefficients[0:47, 7]) * (
                bar ** wageningen_polynomial_coefficients[0:47, 8]) * (
                blade_num ** wageningen_polynomial_coefficients[0:47, 9]))

    return k_t, k_q, j  # TODO this is a bodge, when more time sort out how these functions are split up


def calc_open_water_efficiency(pitch_diameter_ratio, bar, blade_num, prop_diameter, revs_per_second, v_a):
    """
    Calculates the open water efficiency based on using a Wageningen series propeller.See "Ship Resistance and Propulsion - Molland, Turnock & Hudson" Section 16.2.1.1.

    :param pitch_diameter_ratio: Propeller pitch / propeller diameter.
    :param bar: Blade area ratio (expanded blade area / propeller disc area).
    :param blade_num: Number of propeller blades.
    :param prop_diameter: Propeller diameter (m).
    :param revs_per_second: Propeller rate of revolution (rps).
    :param v_a: Wake speed(ship speed * [1 - w_t]) (m/s).
    :return: Open water efficiency.

    #TODO the limits of applicability for wageningen bwl series is num_blades 2-7, BAR 0.3-1.05, P/D 0.6-1.4
    """

    # Import polynomial coefficients for Kt and Kq calculations
    wageningen_polynomial_coefficients = np.genfromtxt('wageningen_polynomial_coefficients.csv', delimiter=',')

    # Calculating J based on propeller speed (RPS) and diameter
    # if (revs_per_second * prop_diameter) == 0:
    # print('Zero erro RPS is  ', revs_per_second, 'prop_diameter is  ', prop_diameter)
    j = v_a / (revs_per_second * prop_diameter)
    if j <= 0:
        print('J is <= 0, rps is ', revs_per_second, 'v_a is  ', v_a)
        return 0
    if j > 1.5:
        print('J > 1.5, j is ', j, ' rps is ', revs_per_second, ' v_a is  ', v_a)
    # Calculating Kq and Kt from the polynomial coefficients
    k_t = np.sum(
        wageningen_polynomial_coefficients[0:39, 0] * (j ** wageningen_polynomial_coefficients[0:39, 1]) * (pitch_diameter_ratio ** wageningen_polynomial_coefficients[0:39, 2]) * (
                bar ** wageningen_polynomial_coefficients[0:39, 3]) * (
                blade_num ** wageningen_polynomial_coefficients[0:39, 4]))
    k_q = np.sum(
        wageningen_polynomial_coefficients[0:47, 5] * (j ** wageningen_polynomial_coefficients[0:47, 6]) * (pitch_diameter_ratio ** wageningen_polynomial_coefficients[0:47, 7]) * (
                bar ** wageningen_polynomial_coefficients[0:47, 8]) * (
                blade_num ** wageningen_polynomial_coefficients[0:47, 9]))

    eta_o = (j * k_t) / (2 * np.pi * k_q)
    return eta_o


def calc_qpc(l, b, t, c_p, c_b, lcb, c_stern, v_s, num_of_props, pitch_diameter_ratio, bar, blade_num, prop_diameter, revs_per_second):
    """

    :param l: Waterline length (m)
    :param b: Moulded breadth (m)
    :param t: Moulded draft (m)
    :param c_p: Prismatic coefficient
    :param c_b: Block coefficient
    :param lcb: Longitudinal centre of buoyancy. This is the position of the centre of buoyancy forward of 0.5 LWL as a percentage of LWL.
    :param c_stern: Coefficient depending on stern type. C_stern values depending on afterbody form are as follows....Pram with gondola: c_stern = -25. V-shaped sections: c_stern = -10. Normal section shape: c_stern = 0. U-shaped sections with Hogner stern: c_stern = 10.
    :param v_s: Ship speed (m/s) # TODO is this correct?
    :param num_of_props: Number of propellers. This can be either 1 or 2.
    :param pitch_diameter_ratio:
    :param bar: Blade area ratio (expanded blade area / propeller disc area).
    :param blade_num: Number of propeller blades.
    :param prop_diameter: Propeller diameter (m).
    :param revs_per_second: Propeller rate of revolution (rps).
    :return: Quasi Propulsive Coefficient

    Returns zero if QPC is above 1 or below zero. # TODO is there a better way of doing this?
    #  TODO the limits of applicability for wageningen bwl series is num_blades 2-7, BAR 0.3-1.05, P/D 0.6-1.4
    """
    # TODO check all these parameters are correct, speed.e. is d the same as prop_diameter?

    t = calc_thrust_deduction_factor(b, l, t, prop_diameter, c_p, lcb, c_stern)
    w_t = calc_wake_fraction(c_b)
    v_a = v_s * (1 - w_t)  # units = m/s

    eta_h = calc_hull_efficiency(t, w_t)
    eta_o = calc_open_water_efficiency(pitch_diameter_ratio, bar, blade_num, prop_diameter, revs_per_second, v_a)
    eta_r = calc_relative_rotative_efficiency(num_of_props, bar, c_p, lcb, pitch_diameter_ratio)

    qpc = eta_h * eta_o * eta_r

    if qpc >= 1 or qpc < 0:
        qpc1 = qpc
        qpc = 0
        print('QPC set to zero, original qpc was   ', qpc1)
        # TODO is there a better way of handling this?

    if (revs_per_second * 60) / v_s < 5:
        qpc = 0  # TODO this is a bodge must remove before final version, or think of a better way of handling this ---> Why is this here???? I think the better way of doing it is setting bounds on acceptable advance coefficient, p/d etc.

    return qpc


def calc_effective_power(qpc, eta_t, installed_power):
    """

    :param qpc: Quasi Propulsive Coefficient also called eta subscript d TODO format this as greek character
    :param eta_t: Transmission efficiency (typically 0.98 for engines aft, 0.95 for geared main engines)
    :param installed_power: Installed power (TODO units?
    :return: Effective power (TODO units?
    """

    return installed_power * eta_t * qpc

# TODO rearrange code to pass n_prop parameter to a higher level function then call functions for single screw thrust deduction and rotative efficiency depending on n_prop
# TODO this is returning really strange values if rpm or speed are out of a certain range. How can these exceptions be handled properly?

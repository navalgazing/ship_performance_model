import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

# TODO why don't these methods need the operating RPM? This is surely needed to get the correct fuel consumption in the case of a heavy or lightly loaded propeller?

def delta_sfoc(pc_mep, gradient, intercept):
    """
    This function will return the change in Specific Fuel Oil Consumption,
    relative to L1 on the engine layout diagram, for a given percentage
    of MEP.
    The equation for the relationship takes the form of a straight line:
    y = mx+c
    Where m and c are obtained from graphs published in MAN Project Guides.
    The values of m and c will vary depending on the engine.
    gradient = m
    intercept = c
    """
    return gradient*pc_mep+intercept


def calc_mep_pc_mp(mcr_p, mcr_rpm, l1_kw_c, l3_kw_c, n_piston, rpm_min, rpm_max):
    """
    Calculate percentage MEP, given MCR and power at L1
    mcr_p = Maximum Continuous Rated Power [kW]
    mcr_rpm = Speed at MCR [rpm]
    l1_kw_c = Power at L1 [kW] per cylinder
    n_piston = number of cylinders
    N.B. One needs to exercise caution when lifting data from loglog plots!
    """
    # calculate slope of constant MEP line from L3 to L1
    mep_slope = np.log10(l1_kw_c/l3_kw_c)/np.log10(rpm_max/rpm_min)

    # calculate where constant MEP line intersects the ordinate axis
    mep_o = mcr_p/pow(mcr_rpm, mep_slope)

    # calculate where the constant MEP line through MP intersects the
    # 100% rpm line
    mep_y = mep_o*pow(rpm_max, mep_slope)

    # calculate the % MEP of point MP (wrt L1)
    mep_pc_mp = 100.0-((l1_kw_c*n_piston-mep_y)/(l1_kw_c*n_piston)*100.0)

    return mep_pc_mp


def sfoc_at_point(sfoc_gradient, sfoc_intercept, pc_mcr, mcr_p, run_p,
                  l1_sfoc, mep_pc_mp):
    """
    Calculate SFOC at the given running point
    MP = power at MCR
    RP = power at running point
    pc_mcr = percent of SCMR corresponding to the MEP/SFOC line
    """

    # calculate curve of SFOC [g/kWh] vs MCR [kW]
    # calculate SFOC at 100% SMCR (i.e. point M) on the engine layout diagram
    x_one = pc_mcr[0]/100.0*mcr_p
    y_one = l1_sfoc + delta_sfoc(mep_pc_mp, sfoc_gradient[0], sfoc_intercept[0])

    # calculate SFOC at 70% SMCR
    x_two = pc_mcr[1]/100.0*mcr_p
    y_two = l1_sfoc+delta_sfoc(mep_pc_mp, sfoc_gradient[1], sfoc_intercept[1])

    # calculate SFOC at 50% SMCR
    x_three = pc_mcr[2]/100.0*mcr_p
    y_three = l1_sfoc+delta_sfoc(mep_pc_mp, sfoc_gradient[2], sfoc_intercept[2])

    # create a spline through the data points
    x_i = np.array([x_three, x_two, x_one])
    y_i = np.array([y_three, y_two, y_one])
    order = 2  # spline order: 1 linear, 2 quadratic, 3 cubic ...
    spline = InterpolatedUnivariateSpline(x_i, y_i, k=order)

    # correction for very low MCR
    # this is a crude percentage difference calculation based upon a spline
    # through the 100, 70 and 50 % points, and the curve provided in MAN
    # project guide
    x_four = 0.4*mcr_p
    y_four = spline(0.4*mcr_p)-0.602/100.0*spline(0.4*mcr_p)

    # calculate new spline, corrected for low MCR points
    xi_corrected = np.array([x_four, x_three, x_two, x_one])
    yi_corrected = np.array([y_four, y_three, y_two, y_one])
    order = 2  # spline order: 1 linear, 2 quadratic, 3 cubic ...
    spline_corrected = InterpolatedUnivariateSpline(xi_corrected, yi_corrected, k=order)

    # calculate running power as a percentage of MCR
    run_p_pc = 100.0-(mcr_p-run_p)/mcr_p*100

    # interpolate spline for value of sfoc
    sfoc = spline_corrected(run_p_pc/100.0*mcr_p)

    return sfoc

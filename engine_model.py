# Whole Ship Model (WSM)
# John Calleya

# Version 0.9.1 (beta)
# - Updated Interface with layout and emissions,
# - Added layout considerations (enough for figure)
# - Added emissions file but not included in model yet
# - Updated performance and degradation file
# Next steps (0.9.2)
# - Add routing and final items to interface (not including new ship definitions)
# - Link-up Added Resistance calculation
# SFOC correction (0.9.3)
# - Corrected SFOC when gearbox is used.
# - Small issue about variables being used in parts of the code.
# - Now passing name of the engine and number of pistons. This is for a better understanding
# with the outputs. Less time will be spent in debugging.
# - Fixed issue with large two-stroke engines not being selected by adding extra cylinders
# in the database.
# Auxiliary engine correction (0.9.4)
# Total installed power was being passed as the total power at 85%MCR and not at 100%MCR
# Auxiliary engine correction (0.9.5)
# Total installed power is entering the two engines section which is adding up the sfoc, not correct for
# Glotram assumptions. General assumptions for auxiliary engines need to be revised.

# Author of GEM.py: David Trodden (Newcastle University)
# GEM.py is a SCC file
# About GEM.py:
# File for engine model that uses MAN database.


# -*- coding: utf-8 -*-
"""
Created on Thu May 14 11:05:48 2015
@author: David Trodden <David.Trodden@ncl.ac.uk>
"""

import sys
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

# TODO have a look at https://link.springer.com/article/10.1007/s00773-017-0523-1 and how they model the engine performance

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
    return gradient * pc_mep + intercept


def mep_pc_l1(mcr_p, mcr_rpm, l1_kw_c, l3_kw_c, n_piston, rpm_min, rpm_max):
    """
    Calculate percentage MEP, given MCR and power at L1
    mcr_p = Maximum Continuous Rated Power [kW]
    mcr_rpm = Speed at MCR [rpm]
    l1_kw_c = Power at L1 [kW] per cylinder
    n_piston = number of cylinders
    N.B. One needs to exercise caution when lifting data from loglog plots!
    """
    # calculate slope of constant MEP line from L3 to L1
    mep_slope = np.log10(l1_kw_c / l3_kw_c) / np.log10(rpm_max / rpm_min)

    # calculate where constant MEP line intersects the ordinate axis
    mep_o = mcr_p / pow(mcr_rpm, mep_slope)

    # calculate where the constant MEP line through MP intersects the
    # 100% rpm line
    mep_y = mep_o * pow(rpm_max, mep_slope)

    # calculate the % MEP of point MP (wrt L1)
    mep_pc_mp = 100.0 - ((l1_kw_c * n_piston - mep_y) / (l1_kw_c * n_piston) * 100.0)

    return mep_pc_mp


def sfoc_at_point(sfoc_gradient, sfoc_intercept, pc_mcr, mcr_p, run_p, l1_sfoc, mep_pc_mp):
    """
    Calculate FOC at the given running point
    MP = power at MCR
    RP = power at running point
    pc_mcr = percent of SCMR corresponding to the MEP/SFOC line
    """

    # calculate curve of SFOC [g/kWh] vs MCR [kW]
    #
    # calculate SFOC at 100% SMCR (speed.e. point M) on the engine layout diagram
    x_one = pc_mcr[0] / 100.0 * mcr_p
    y_one = l1_sfoc + delta_sfoc(mep_pc_mp, sfoc_gradient[0], sfoc_intercept[0])

    # calculate SFOC at 70% SMCR
    x_two = pc_mcr[1] / 100.0 * mcr_p
    y_two = l1_sfoc + delta_sfoc(mep_pc_mp, sfoc_gradient[1], sfoc_intercept[1])

    # calculate SFOC at 50% SMCR
    x_three = pc_mcr[2] / 100.0 * mcr_p
    y_three = l1_sfoc + delta_sfoc(mep_pc_mp, sfoc_gradient[2], sfoc_intercept[2])

    # create a spline through the data points
    x_i = np.array([x_three, x_two, x_one])
    y_i = np.array([y_three, y_two, y_one])
    order = 2  # spline order: 1 linear, 2 quadratic, 3 cubic ...
    spline = InterpolatedUnivariateSpline(x_i, y_i, k=order)

    # correction for very low MCR
    # this is a crude percentage difference calculation based upon a spline
    # through the 100, 70 and 50 % points, and the curve provided in MAN
    # project guide
    x_four = 0.4 * mcr_p
    y_four = spline(0.4 * mcr_p) - 0.602 / 100.0 * spline(0.4 * mcr_p)
    # y_four = s(0.4*MP) # No correction!

    # calculate new spline, corrected for low MCR points
    xi_corrected = np.array([x_four, x_three, x_two, x_one])
    yi_corrected = np.array([y_four, y_three, y_two, y_one])
    order = 2  # spline order: 1 linear, 2 quadratic, 3 cubic ...
    spline_corrected = InterpolatedUnivariateSpline(xi_corrected, yi_corrected, k=order)

    # calculate running power as a percentage of MCR
    run_p_pc = 100.0 - (mcr_p - run_p) / mcr_p * 100

    # interpolate spline for value of sfoc
    sfoc = spline_corrected(run_p_pc / 100.0 * mcr_p)
    return sfoc


class PoweringSpecs:
    """
    Propulsion and Auxiliary Engine Powering Calculations
    """

    def __init__(self, q_trial, rpm_trial, q_run, rpm_run, hotel_load_design,
                 hotel_load_service, pto, eta_pto, cpp, sea_margin,
                 main_engine_margin, light_running_factor, aux_engine_margin,
                 main_engine_type, aux_engine_type, fuel_type_main,
                 fuel_type_aux, green_technologies):

        # attributes that are passed through initialisation
        self.q_trial = q_trial
        self.rpm_trial = rpm_trial
        self.q_run = q_run
        self.rpm_run = rpm_run
        self.hotel_load_design = hotel_load_design
        self.hotel_load_service = hotel_load_service
        self.pto = pto
        self.eta_pto = eta_pto
        self.cpp = cpp
        self.sea_margin = sea_margin
        self.main_engine_margin = main_engine_margin
        self.f_lr = light_running_factor
        self.aux_engine_margin = aux_engine_margin
        self.main_engine_type = main_engine_type
        self.aux_engine_type = aux_engine_type
        self.fuel_type_main = fuel_type_main
        self.fuel_type_aux = fuel_type_aux
        self.green_technologies = green_technologies

        # TODO include shaft and gearbox efficiency as initialisation variables
        self.eta_shaft = 0.97  # shafting efficiency factor
        self.eta_gearbox = 0.96  # gearbox efficiency factor

        self.gearbox_ratio = 1.0

        # attributes that are calculated
        self.rpm_service = 0.0
        self.rpm_mcr = 0.0
        self.p_run = 0.0
        self.p_trial = 0.0
        self.p_mcr = 0.0
        self.p_service = 0.0

        self.hotel_load_max = 0.0
        self.total_aux_engines_required = 0
        self.total_aux_mcr = 0.0

        self.main_engine_designation = "not assigned"
        self.aux_engine_designation = "not assigned"

        self.sfoc_main_at_run = 0.0
        self.sfoc_aux_at_run = 0.0

        self.specific_co2_main_engine = 0.0
        self.specific_co2_aux_engine = 0.0
        self.specific_co2_total = 0.0
        self.co2_main_engine = 0.0
        self.co2_aux_engine = 0.0
        self.co2_total = 0.0

        self.main_engine_mass = 0.0
        self.main_engine_volume = 0.0
        self.aux_engine_mass = 0.0
        self.aux_engine_volume = 0.0

        # changes in emissions from use of different green technologies
        self.delta_CO2 = 0.0
        self.delta_SOX = 0.0
        self.delta_NOX = 0.0

        # some input checks
        if self.main_engine_type < 1 or self.main_engine_type > 2:
            print("ERROR: Main Engine Type format incorrect. Aborting...")
            sys.exit(1)

    def main_engine_requirements(self):
        """
        This method calculates main propulsion engine requirements from
        given input
        """

        # calculate propeller design point, in kW.
        # This is the trial-conditions, or calm water scenario
        # This is equivalent to light-running, or propeller design point, PD
        # This includes the power required by any shaft generator (PTO)
        p_delivered_trial = self.q_trial * 2.0 * np.pi * self.rpm_trial / 60.0 + self.pto

        # this is the delivered power to the propeller, so need to account for
        # shafting, gear-box and pto if fitted
        self.p_trial = p_delivered_trial + (1.0 - self.eta_shaft) * p_delivered_trial

        if self.pto != 0.0:  # pto fitted
            self.p_trial = self.p_trial + (1.0 - self.eta_pto) * self.p_trial
        # p.trial is now brake power required for trial conditions

        # calculate light-running propeller curve
        # this is based on assuming power is proportional to rpm^3
        # this could be revised for different ship types
        # c.f. "Basic Principles of Ship Propulsion" by MAN B&W
        # This is not actually used here!
        # CL = PT/(pow(RPMT, 3.0))

        # calculate heavy-running propeller curve
        # this is based on assuming power is proportional to rpm^3
        # this could be revised for different ship types
        # c.f. "Basic Principles of Ship Propulsion" by MAN B&W
        heavy_prop_c = self.p_trial / (pow((self.rpm_trial - self.f_lr * self.rpm_trial), 3.0))

        # calculate service propulsion point, in kW, as stipulated by the
        # sea-margin
        # N.B. this is different to the actual current running service point
        self.p_service = heavy_prop_c * pow(self.rpm_trial, 3.0) + self.sea_margin * heavy_prop_c * pow(self.rpm_trial, 3.0)

        # calculate the rpm at the service propulsion point
        # this is based on assuming power is proportional to rpm^3
        # this could be revised for different ship types
        # c.f. "Basic Principles of Ship Propulsion" by MAN B&W
        self.rpm_service = pow((self.p_service / heavy_prop_c), (1.0 / 3.0))

        # calculate specified maximum continuous rated power (SMCR), as
        # stipulated by the engine-margin
        self.p_mcr = self.p_service * (self.main_engine_margin + 1.0)

        # calculate rpm at SMCR
        # this is based on assuming power is proportional to rpm^3
        # this could be revised for different ship types
        # c.f. "Basic Principles of Ship Propulsion" by MAN B&W
        self.rpm_mcr = pow((self.p_mcr / heavy_prop_c), (1.0 / 3.0))

        # calculate power in current running conditions
        self.p_run = self.q_run * 2.0 * np.pi * self.rpm_run / 60.0 + self.pto

    def estimate_co2(self):
        """
        :rtype: object
        :return: TODO complete this
        """
        # CO_2 emission factor based upon stoichiometric combustion
        # CO2_EF = 3.1141
        # TODO remove bit that uses fuel list - or change in some way
        co2_main_emission_factor = fuel_list[self.fuel_type_main][2]
        co2_aux_emission_factor = fuel_list[self.fuel_type_aux][2]

        # estimate specific CO_2 emissions [g/kWh]
        self.specific_co2_main_engine = self.sfoc_main_at_run * co2_main_emission_factor
        self.specific_co2_aux_engine = self.sfoc_aux_at_run * co2_aux_emission_factor
        self.specific_co2_total = self.specific_co2_main_engine + self.specific_co2_aux_engine

        # calculate fuel oil consumption [draft/h]
        foc_main = self.sfoc_main_at_run * self.p_run / 1.0E+06
        foc_aux = self.sfoc_aux_at_run * self.hotel_load_service / 1.0E+06

        # estimate absolute CO_2 emissions [draft/h]
        self.co2_main_engine = foc_main * co2_main_emission_factor
        self.co2_aux_engine = foc_aux * co2_aux_emission_factor
        self.co2_total = self.co2_main_engine + self.co2_aux_engine

    def delta_emissions(self):
        """
        Estimate change in emissions through use of various technologies.
        Changes are in %, so for example, total reduction in NOX would be:
        total_reduction_nox = nox_emissions + delta_nox*nox_emissions
        source for estimations: http://cleantech.cnss.no/
        values are averaged unless otherwise stated
        :rtype: object
        """

        # deltas are cumulative so extreme caution is needed to avoid mutually
        # exclusive technologies!
        # TODO make checks on the above
        if self.green_technologies[0]:  # Low Sulphur Fuel
            self.delta_SOX = self.delta_SOX + (-0.800)

        if self.green_technologies[1]:  # Scrubber fitted?
            self.delta_SOX = self.delta_SOX + (-0.925)

        if self.green_technologies[2]:  # Direct Water Injection
            self.delta_NOX = self.delta_NOX + (-0.600)
            self.delta_CO2 = self.delta_CO2 + (+0.02)

        if self.green_technologies[3]:  # Exhaust Gas Recirculation
            self.delta_NOX = self.delta_NOX + (-0.525)

        if self.green_technologies[4]:  # Selective Catalytic Reduction
            self.delta_NOX = self.delta_NOX + (-0.945)

        if self.green_technologies[5]:  # Humid Air Motor
            self.delta_NOX = self.delta_NOX + (-0.500)

        if self.green_technologies[6]:  # Combustion Air Saturation System
            self.delta_NOX = self.delta_NOX + (-0.450)

        if self.green_technologies[7]:  # Water in fuel emulsion
            self.delta_NOX = self.delta_NOX + (-0.350)

        if self.green_technologies[8]:  # Internal Engine Modification
            self.delta_NOX = self.delta_NOX + (-0.350)

    def dimensions(self):
        """
        Estimate mass and volume of main and aux engines.
        Source: "Design of Propulsion and Electric Power Generation Systems"
                H.K. Woud and D. Stapersma
                IMarEST, 2002.
        Values are averaged from source.
        """

        # main engine
        if self.main_engine_type == 1:  # slow speed two-stroke
            self.main_engine_mass = 38.5 * self.p_mcr  # [kg]
            self.main_engine_volume = 0.035 * self.p_mcr  # [m^3]
        elif self.main_engine_type == 2:  # four-stroke medium speed
            self.main_engine_mass = 12.5 * self.p_mcr  # [kg]
            self.main_engine_volume = 0.016 * self.p_mcr  # [m^3]
        else:  # not implemented yet
            print("ERROR: Main engine type not yet implemented in dimension calculations. Aborting...")
            sys.exit(1)

        # aux engine(s)
        if self.aux_engine_type == 1:  # slow speed two-stroke (!)
            self.aux_engine_mass = 38.5 * self.total_aux_mcr  # [kg]
            self.aux_engine_volume = 0.035 * self.total_aux_mcr  # [m^3]
        elif self.aux_engine_type == 2:  # four-stroke high speed (different to medium!)
            self.aux_engine_mass = 4.15 * self.total_aux_mcr  # [kg]
            self.aux_engine_volume = 0.0054 * self.total_aux_mcr  # [m^3]
        else:  # not implemented yet
            print("ERROR: Auxiliary engine type not yet implemented in dimension calculations. Aborting...")
            sys.exit(1)


def run_gem(case_study, q_trial, rpm_trial, q_run, rpm_run, hotel_load_design,
            hotel_load_service, pto, eta_pto, cpp, sea_margin,
            main_engine_margin, light_running_factor, aux_engine_margin,
            main_engine_type, aux_engine_type, fuel_type_main, fuel_type_aux,
            green_technologies):
    # TODO what is the variable "case study"?
    """
    This function is the main routine which pulls together all the functions
    which make up the Generic Engine Model
    """

    # ship1 instance
    ship1 = PoweringSpecs(q_trial, rpm_trial, q_run, rpm_run,
                          hotel_load_design, hotel_load_service, pto, eta_pto,
                          cpp, sea_margin, main_engine_margin,
                          light_running_factor, aux_engine_margin,
                          main_engine_type, aux_engine_type, fuel_type_main,
                          fuel_type_aux, green_technologies)

    # calculate main engine powering requirements
    ship1.main_engine_requirements()

    # estimate C0_2 emissions
    ship1.estimate_co2()

    # estimate deltas in emissions from using different green technologies
    ship1.delta_emissions()

    # estimate mass and volume of main and auxiliary engine(s)
    ship1.dimensions()

    # return necessary output
    gem_output = [ship1.p_mcr,  # installed main engine power MCR [kW]
                  ship1.sfoc_main_at_run,  # main engine SFOC at running point [draft/kWh]
                  ship1.main_engine_mass,  # mass of main engine [kg]
                  ship1.main_engine_volume,  # volume of main engine [m^3]
                  ship1.specific_co2_main_engine,  # specific CO2 emissions from main engine [g/kWh]
                  ]

    return gem_output

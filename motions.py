import numpy as np


def calculate_jensen_acceleration(lwl, bwl, draft, vol_disp, velocity, heading=0, wave_amplitude=1, long_position=0): # TODO change this to take the heading from somewhere
    # from Jensen et al. 2004 https://www.sciencedirect.com/science/article/pii/S0029801803001082 and Henrique's code
    #  on view-source:http://shiplab.hials.org/app/shipmotion/

    g = 9.807
    f_n = velocity / np.sqrt(g * lwl)
    wave_frequency = []
    vertical_displacement = []
    vertical_acceleration = []
    beta = heading * np.pi / 180
    c_b = vol_disp / (lwl * bwl * draft)
    b = bwl * c_b # here, b is a variable specific to the method given by Jensen. Explained in his paper in equation 3.1



    for i in range(0,1000):
        wave_frequency.append(0.05 + 1.95*i/1000)  # from Henrique's code, I think this is to calculate results for
        # frequencies between 0 and 2
        wave_number = (wave_frequency[i] ** 2) / g  # this is k in the Jensen paper
        eff_wave_number = np.abs(wave_number * np.cos(beta))
        smith_factor = np.exp(-1 * eff_wave_number * draft)
        alpha = 1 - f_n * np.sqrt(wave_number * lwl) * np.cos(beta)
        encounter_frequency = alpha * wave_frequency[i]
        sectional_hydro_damping = 2 * np.sin(0.5 * wave_number * b * alpha ** 2) * np.exp(-1 * wave_number)  # this is
        #  A in the Jensen paper
        f = np.sqrt(((1 - wave_number * draft) ** 2) + ((sectional_hydro_damping ** 2) / (wave_number * b * alpha ** 3)) ** 2)

        # Below are F and G the forcing functions
        F = smith_factor * f * (2 / (eff_wave_number * lwl)) * np.sin(eff_wave_number * lwl / 2)
        G = smith_factor * f * (24 / ((eff_wave_number * lwl) ** 2) * lwl) * (np.sin(eff_wave_number * lwl / 2) - (eff_wave_number * lwl / 2) * np.cos(eff_wave_number * lwl / 2))

        eta = np.sqrt((1 - 2 * wave_number * draft * alpha ** 2) ** 2 + (sectional_hydro_damping ** 2 / wave_number * b * alpha ** 2) ** 2) ** -1

        frf_heave = eta * F * wave_amplitude  # Henrique's code multiplies this by wave amplitude but Jensen doesn't in his paper, why?
        frf_pitch = eta * G * wave_amplitude  # Henrique's code multiplies this by wave amplitude but Jensen doesn't in his paper, why?

        vertical_displacement.append(np.sqrt(frf_heave ** 2 + (long_position ** 2 * frf_pitch ** 2)))
        vertical_acceleration.append(encounter_frequency ** 2 * vertical_displacement[i])

    return max(vertical_acceleration)


#  TODO
#  Add contraints on dimensions/dimensional ratios from Jensen paper with a suitable error handling procedure and user popup
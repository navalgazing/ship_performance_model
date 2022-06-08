from desdeo-problem import variable_builder
from desdeo-problem import ScalarObjective
from desdeo-problem import MOProblem
from main import calc_p_i_required #import the objective functions
from motions import calculate_jensen_acceleration





var_names = []
initial_values = []
lower_bounds = []
upper_bounds = []

variables = variable_builder(var_names, initial_values, lower_bounds, upper_bounds)

f1 = ScalarObjective("power consumption", calc_p_i_required, maximize=False)
f2 = ScalarObjective("vertical acceleration", calculate_jensen_acceleration, maximize=False)
#f3 = ScalarObjective("speed through water", )

problem = MOProblem(objectives=[f1, f2], variables=variables)

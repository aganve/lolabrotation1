import matplotlib
import lib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from necro_uncal_new import model
import scipy.interpolate
from pysb.integrate import *
from simplepso.pso import PSO
import collections​
model.enable_synth_deg()

# Declaring our observable: phosphorylated MLKL
obs_names = ['MLKLp_obs']
mlklp_obs = 'MLKLp_obs'

# Defining a few helper functions to use
def normalize(trajectories):
    """Rescale a matrix of model trajectories to 0-1"""
    ymin = trajectories.min(0)
    ymax = trajectories.max(0)
    return (trajectories - ymin) / (ymax - ymin)

# Linspace refers to generating linearly spaced values
# Creating a vector of complex numbers of 11 evenly spaced points between 0 and 960
# E.g., np.linspace(1.0, 5.0, num=10)
t = np.linspace(0, 960, num=11)
solver1 = ScipyOdeSimulator(model, tspan=t)

# Creating a grid of values, all of the same type, and is indexed by a list of non-negative integers (i.e., tuple)
x10 = np.array([0,2,4,6,8,10,12,14,16,18, 20])
y10 = np.array([0., 0.001, 0.02, 0.03, 0.04, 0.06, .09, .21, .40, .65, .81])

ydata_norm = y10

rate_params = model.parameters_rules()
param_values = np.array([p.value for p in model.parameters])
rate_mask = np.array([p in rate_params for p in model.parameters])

#
original_values = np.array([p.value for p in model.parameters])

# We search in log10 space for the parameters
# We will use a best guess starting position for the model, up or down 1 order of magnitude
log10_original_values = np.log10(original_values[rate_mask])

# Here we define the cost function, which finds a line of best fit for our inputs (X values) and outputs (Y values)
# The event we are finding the cost of is the difference between estimated values, or the difference between the hypothesis and the real values
# Perform analysis using the squared error cost function
# Returns a tuple: a sequence that cannot be changed unlike lists; uses parantheses
def obj_function(params):
    params_tmp = np.copy(params)
    param_values[rate_mask] = 10 ** params_tmp  # don't need to change
    result = solver1.run(param_values=param_values)
    ysim_array1 = result.observables['MLKLp_obs'][:]
    ysim_norm1 = normalize(ysim_array1)
​
    e1 = np.sum((ydata_norm - ysim_norm1) ** 2)
​
    return e1,

  counter = 0
    # Here, we initial the class
    # We must provide the cost function and a starting value
    # Provide the bounds, speed, particle, and iterations count
    for i in range(1000):
        optimizer = PSO(cost_function=obj_function,start = log10_original_values, verbose=True)
        # We also must set bounds. This can be a single scalar or an array of len(start_position)
        optimizer.set_bounds(parameter_range=2)
        optimizer.set_speed(speed_min=-.25, speed_max=.25)
        optimizer.run(num_particles=75, num_iterations=25)
        best_pars[i] = optimizer.best
        print(optimizer.best)
        # print(i, counter)
    np.save('optimizer_best_75_50_100TNF',best_pars)
​
if '__main__' == __name__:
    run_example()





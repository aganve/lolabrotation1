from pylab import *
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
from pysb.macros import catalyze
import matplotlib.pyplot as plt
import numpy as np
from pysb.simulator import ScipyOdeSimulator
from necroptosis import model
import pandas as pd

# print('model species')
# print(model.species)
# print('model parameters')
# print(model.parameters)
# print('model odes')
# print(list(model.odes))

# load optimized parameters at 25 particles, and TNF conc of 100 ng; calibration ran on 9/27/19
opt_params = np.load('necro_optimizer_best_25_100_927_TNF100.npy')
# define length of optimized parameters
n_pars = len(opt_params)
all_pars = np.zeros((n_pars, len(model.parameters)))

# rate parameters are model parameter rules (i.e., kinetic rate constants)
rate_params = model.parameters_rules()
# parameter values are an array of model parameters (i.e., initial conditions and kinetic rate constants)
param_values = np.array([p.value for p in model.parameters])
# mask rate parameters (i.e., kinetic rate constants) from our model parameters
rate_mask = np.array([p in rate_params for p in model.parameters])

# for every loop...
for i in range(len(opt_params)):
    par = opt_params[i]
    param_values[rate_mask] = 10 ** par #unlogging param values
    all_pars[i] = param_values

# print(all_pars[:100])
# all_pars = all_pars[:5]
x100 = np.array([30, 90, 270, 480, 600, 720, 840, 960]) # time in minutes
y100 = np.array([0.00885691708746097,0.0161886154261265,0.0373005242261882,0.2798939020159581,0.510, .7797294067, 0.95,1]) # normalized values

tspan = np.linspace(0, 1440, 101) # time span of simulation (start, stop, step)
# Scipy is a Python library for scientific computing. Scipy is the solver for our model
result = ScipyOdeSimulator(model, tspan=tspan).run(param_values=all_pars) # run solver of model
# Simulation results
df = result.dataframe

plt.figure(figsize=(15,10))
plt.subplot(241)
for i in range(len(all_pars)):
    plt.plot(tspan, df.loc[i]['TNF_obs'].iloc[:], lw = 1.5, label = 'TNF_obs') #plot observable
# plt.scatter(x100, y100)
plt.xlabel('Time min', fontsize=14)
plt.ylabel('Amount produced [molecules]', fontsize=14)
plt.title('TNF_obs')
# plt.legend(loc = 'best', fontsize = 12) #add legend to plot
# plt.savefig('E+S_to_P.pdf') #to save figure

# plt.figure()
plt.subplot(242)
for i in range(len(all_pars)):
    plt.plot(tspan, df.loc[i]['CI_inactive_obs'].iloc[:], lw = 1.5, label = 'CI_inactive_obs') #plot observable
# plt.scatter(x100, y100)
plt.xlabel('Time min', fontsize=14)
plt.ylabel('Amount of produced [molecules]', fontsize=14)
plt.title('CI_inactive_obs')
# plt.legend(loc = 'best', fontsize = 12) #add legend to plot
# plt.savefig('E+S_to_P.pdf') #to save figure

# plt.figure()
plt.subplot(243)
for i in range(len(all_pars)):
    plt.plot(tspan, df.loc[i]['TNF_CI_CIi_obs'].iloc[:], lw = 1.5, label = 'TNF_CI_CIi_obs') #plot observable
# plt.scatter(x100, y100)
plt.xlabel('Time min', fontsize=14)
plt.ylabel('Amount of produced [molecules]', fontsize=14)
plt.title('TNF_CI_CIi_obs')
# plt.legend(loc = 'best', fontsize = 12) #add legend to plot
# plt.savefig('E+S_to_P.pdf') #to save figure

# plt.figure()
plt.subplot(244)
for i in range(len(all_pars)):
    plt.plot(tspan, df.loc[i]['CIa_obs'].iloc[:], lw = 1.5, label = 'CIa_obs') #plot observable
# plt.scatter(x100, y100)
plt.xlabel('Time min', fontsize=14)
plt.ylabel('Amount of produced [molecules]', fontsize=14)
plt.title('CIa_obs')
# plt.legend(loc = 'best', fontsize = 12) #add legend to plot
# plt.savefig('E+S_to_P.pdf') #to save figure

# plt.figure()
plt.subplot(245)
for i in range(len(all_pars)):
    plt.plot(tspan, df.loc[i]['CIIa_obs'].iloc[:], lw = 1.5, label = 'CIIa_obs') #plot observable
# plt.scatter(x100, y100)
plt.xlabel('Time min', fontsize=14)
plt.ylabel('Amount of produced [molecules]', fontsize=14)
plt.title('CIIa_obs')
# plt.legend(loc = 'best', fontsize = 12) #add legend to plot
# plt.savefig('E+S_to_P.pdf') #to save figure

# plt.figure()
plt.subplot(246)
for i in range(len(all_pars)):
    plt.plot(tspan, df.loc[i]['CII_MLKL_bind_obs'].iloc[:], lw = 1.5, label = 'CII_MLKL_bind_obs') #plot observable
# plt.scatter(x100, y100)
plt.xlabel('Time min', fontsize=14)
plt.ylabel('Amount of produced [molecules]', fontsize=14)
plt.title('CII_MLKL_bind_obs')
# plt.legend(loc = 'best', fontsize = 12) #add legend to plot
# plt.savefig('E+S_to_P.pdf') #to save figure

# plt.figure()
plt.subplot(247)
for i in range(len(all_pars)):
    plt.plot(tspan, df.loc[i]['MLKLu_obs'].iloc[:], lw = 1.5, label = 'MLKLu_obs') #plot observable
# plt.scatter(x100, y100)
plt.xlabel('Time min', fontsize=14)
plt.ylabel('Amount of produced [molecules]', fontsize=14)
plt.title('MLKLu_obs')
# plt.legend(loc = 'best', fontsize = 12) #add legend to plot
# plt.savefig('E+S_to_P.pdf') #to save figure

# plt.figure()
plt.subplot(248)
for i in range(len(all_pars)):
    plt.plot(tspan, df.loc[i]['MLKLp_obs'].iloc[:], lw = 1.5, label = 'MLKLp_obs') #plot observable
# plt.scatter(x100, y100)
plt.xlabel('Time min', fontsize=14)
plt.ylabel('Amount of produced [molecules]', fontsize=14)
plt.title('MLKLp_obs')
# plt.legend(loc = 'best', fontsize = 12) #add legend to plot
# plt.savefig('E+S_to_P.pdf') #to save figure
plt.show()

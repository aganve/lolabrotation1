from pylab import *
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
from pysb.macros import catalyze
import matplotlib.pyplot as plt
import numpy as np
from pysb.simulator import ScipyOdeSimulator

# instantiate a model
Model()

# declare monomers
Monomer('TNF', ['bCI'])
Monomer('CI', ['bTNF','state'], {'state':['i','a']})
Monomer('CII', ['bMLKL','state'], {'state': ['a']})
Monomer('MLKL', ['bCII', 'state'], {'state': ['u','p']})

# input the parameter values
Parameter('kf1', 1.0e-6)
Parameter('kr1', 1.0e-3)
Parameter('kf2', 1.0e-6)
Parameter('kf3', 1.0e-6)
Parameter('kf4', 1.0e-6)
Parameter('kr4', 1.0e-3)
Parameter('kf5', 1.0e-6)

# now input the rules
Rule('TNF_CI_bind_CI_inactive', TNF(bCI=None) + CI(bTNF=None, state='i') | TNF(bCI=1) % CI(bTNF=1, state='i'), kf1, kr1)
Rule('TNF_CI_dissociate_CI_active', TNF(bCI=1) % CI(bTNF=1, state='i') >> TNF(bCI=None) + CI(bTNF=None, state='a'), kf2)
Rule('CI_CII_conversion', CI(bTNF=None, state='a') >> CII(bMLKL=None, state='a'), kf3)
Rule('CII_MLKL_bind', CII(bMLKL=None, state='a') + MLKL(bCII=None, state='u') | CII(bMLKL=1, state='a') % MLKL(bCII=1, state='u'), kf4, kr4)
Rule('CII_MLKL_dissociate_MLKL_phosphorylated', CII(bMLKL=1, state='a') % MLKL(bCII=1, state='u') >> MLKL(bCII=None, state='p'), kf5)

# initial conditions
Parameter('TNF_0', 2326)
Parameter('CI_0', 9000)
Parameter('MLKL_0', 5544)
Initial(TNF(bCI=None), TNF_0)
Initial(CI(bTNF=None, state='i'), CI_0)
Initial(MLKL(bCII=None, state='u'), MLKL_0)

# observable
Observable('TNF_obs', TNF(bCI=None))
Observable('CI_inactive_obs', CI(bTNF=None, state='i'))
Observable('TNF_CI_obs', TNF(bCI=1) % CI(bTNF=1, state='i'))
Observable('CIa_obs', CI(bTNF=None, state='a'))
Observable('CIIa_obs', CII(bMLKL=None, state='a'))
Observable('CII_MLKL_bind_obs', CII(bMLKL=1, state='a') % MLKL(bCII=1, state='u'))
Observable('MLKLu_obs', MLKL(bCII=None, state='u'))
Observable('MLKLp_obs', MLKL(bCII=None, state='p'))

generate_equations(model)

tspan = np.linspace(0, 1440, 1441) # time span of simulation (start, stop, step)
result = ScipyOdeSimulator(model, tspan=tspan).run() # run solver of model

# plt.figure()
# plt.plot(tspan, result.observables['MLKLp_obs'][:], lw = 1.5, label = 'MLKLp_obs') #plot observable
# plt.xlabel('Time [hours]', fontsize=14)
# plt.ylabel('Phosphorylated MLKL amount [molecules]', fontsize=14)
# plt.title('Amount of MLKLp over time')
# plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# plt.savefig('pMLKL to paramset14.pdf') # to save figure
# plt.show()


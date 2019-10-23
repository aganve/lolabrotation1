from pylab import *
from pysb.core import * # brings in all of the Python classes needed to define a PySB model
from pysb.bng import *
from pysb.integrate import *
from pysb.macros import equilibrate, catalyze
import matplotlib.pyplot as plt
import numpy as np
from pysb.simulator import ScipyOdeSimulator # Simulating an ODE with SciPy

# instantiate a model
Model()

# Cys2 = oxidized dimer form of the amino acid cysteine
# XC = System XC: cystine/glutamate transporter located in the PM
# CR = cystine reductase; reduces Cys2 to Cys in the PM
# Cys = cysteine
# Glu_Cys = glutamyl-cysteine; intermediate produced in the sequence of reactions converting cysteine to glutathione (GHS)
# GSS = glutathione synthetase; converts Glut_Cys to GSH
# GSH = reduced glutathione; cofactor for GPX4 enzyme that reduces lipid peroxides to lipid alcohols
# GSSG = oxidized glutathione; GPX4 oxidizes GSH to GSSG when reducing lipid peroxides/H2O2 to lipid alcohols/H2O

# declare monomers
Monomer('Cys2', ['bCR']) # binds cystine reductase (CR) which reduces 1 molecule of cystine to 2 molecules of cysteine
Monomer('CR', ['bCys2'])
Monomer('Cys', ['bGCL'])
Monomer('GCL', ['bCys'])
Monomer('GlutCys', ['bGSS'])
Monomer('GSS', ['bGlutCys'])
Monomer('GSH', ['b1GPX4', 'state'], {'state': ['red', 'ox']}) # reduced glutathione (GSH); oxidized form of GSH is GSSG
Monomer('GPX4', ['b1GSH', 'b2Liperox', 'state'], {'state': ['i', 'a']})
Monomer('Liperox', ['b2GPX4', 'state'], {'state': ['rad', 'red']}) # radical and reduced form

# One main 3d compartment
# Main 3d compartment is environment (Env) and three membranes inside: plasma membrane (PM), cytosol (Cyto) and lysosome (Lysos)
Parameter('VEnv', 1)
Parameter('VPM', 1)
Parameter('VCyto', 1)
Compartment('Env', None, 3, VEnv) # Env doesn't have a parent; it's 3D in volume
Compartment('PM', Env, 2, VPM) # PM's parent is the Env; it's 2D in volume
Compartment('Cyto', PM, 3, VCyto) # Cyto's parent is the PM; it's 3D in volume; proteins can move in the x, y, and z direction, therefore 3D

# input the parameter values
Parameter('kf1', 1.0e-6) # kf is a slower reaction than kr
Parameter('kr1', 1.0e-3)
Parameter('kf2', 1.0e-6)
Parameter('kr2', 1.0e-3) # kr is a faster reaction than kf
Parameter('kc2', 1.0e-1) # catalysis occurs more quickly than forward and reverse reactions
Parameter('kf3', 1.0e-6)
Parameter('kr3', 1.0e-3)
Parameter('kc3', 1.0e-1)
Parameter('kf4', 1.0e-6)
Parameter('kr4', 1.0e-3)
Parameter('kc4', 1.0e-1)
Parameter('kf5', 1.0e-6)
Parameter('kr5', 1.0e-3)
Parameter('kf6', 1.0e-6)
Parameter('kr6', 1.0e-3)
Parameter('kf7', 1.0e-6)

# now input the rules
# transport Cys2 (cystine) from the Env into the Cyto
equilibrate(Cys2(bCR=None) ** Env, Cys2(bCR=None) ** Cyto, [kf1, kr1])
# Rule('Cys2_transport_via_XC', Cys2() ** Env + XC() ** PM >> Cys2() ** Cyto + XC() ** PM, kf1)

# once in the Cyto, Cys2 (cystine) is reduced to Cys via CR catalysis
catalyze(CR() ** Cyto, 'bCys2', Cys2() ** Cyto, 'bCR', Cys(bGCL=None) ** Cyto, [kf2, kr2, kc2])
# Rule('Cys2_CR_bind', Cys2(bCR=None) ** Cyto + CR(bCys2=None) ** Cyto | Cys2(bCR=1) ** Cyto % CR(bCys2=1) ** Cyto, kf2, kr2)

# convert Cys to GlutCys via GCL catalysis in the Cyto
catalyze(GCL() ** Cyto, 'bCys', Cys() ** Cyto, 'bGCL', GlutCys(bGSS=None) ** Cyto, [kf3, kr3, kc3])
# Rule('Cys_GSH_conversion', Cys(bCys2=None) ** Cyto | GSH(bGPX4=None) ** Cyto, kf3, kr2)

# convert Glut-Cys to GSH via GSS catalysis in the Cyto
catalyze(GSS() ** Cyto, 'bGlutCys', GlutCys() ** Cyto, 'bGSS', GSH(b1GPX4=None, state='red') ** Cyto, [kf4, kr4, kc4])

# GPX4 binds GSH and is activated; GSH stays in its reduced formed
Rule('GSH_binds_GPX4_activated', GSH(b1GPX4=None, state='red') ** Cyto + GPX4(b1GSH=None, b2Liperox=None, state='i') ** Cyto | GSH(b1GPX4=1, state='red') ** Cyto % GPX4(b1GSH=1, b2Liperox=None, state='a') ** Cyto, kf5, kr5)

# GPX4-GSH complex binds lipid peroxide; GSH is oxidized to GSSG in the reaction and liperox is reduced to a lipid alcohol
Rule('Liperox_binds', GSH(b1GPX4=1, state='red') ** Cyto % GPX4(b1GSH=1, b2Liperox=None, state='a') ** Cyto + Liperox(b2GPX4=None, state='rad') ** Cyto | GSH(b1GPX4=1, state='ox') ** Cyto % GPX4(b1GSH=1, b2Liperox=2, state='a') ** Cyto % Liperox(b2GPX4=2, state='red') ** Cyto, kf6, kr6)

# GPX4-GSH-peroxide complex dissociates
Rule('GPX4_GSH_liperox_dissociate', GSH(b1GPX4=1, state='ox') ** Cyto % GPX4(b1GSH=1, b2Liperox=2, state='a') ** Cyto % Liperox(b2GPX4=2, state='red') ** Cyto >> GSH(b1GPX4=None, state='ox') ** Cyto + GPX4(b1GSH=None, b2Liperox=None, state='i') ** Cyto + Liperox(b2GPX4=None, state='red') ** Cyto, kf7)

# Initial Conditions
Parameter('Cys2_0', 2300)
Parameter('CR_0', 2300)
Parameter('GCL_0', 2300)
Parameter('GSS_0', 2300)
Parameter('GPX4_0', 1000)
Parameter('Liperox_0', 2300)
Initial(Cys2(bCR=None) ** Env, Cys2_0)
Initial(CR(bCys2=None) ** Cyto, CR_0)
Initial(GCL(bCys=None) ** Cyto, GCL_0)
Initial(GSS(bGlutCys=None) ** Cyto, GSS_0)
Initial(GPX4(b1GSH=None, b2Liperox=None, state='i') ** Cyto, GPX4_0)
Initial(Liperox(b2GPX4=None, state='rad') ** Cyto, Liperox_0)

# Observables
Observable('Cys2_Env_obs', Cys2(bCR=None) ** Env)
Observable('Cys2_Cyto_obs', Cys2(bCR=None) ** Cyto)
Observable('CR_obs', CR(bCys2=None) ** Cyto)
Observable('Cys_obs', Cys(bGCL=None) ** Cyto)
Observable('GCL_obs', GCL(bCys=None) ** Cyto)
Observable('Glut_Cys_obs', GlutCys(bGSS=None) ** Cyto)
Observable('GSS_obs', GSS(bGlutCys=None) ** Cyto)
Observable('GSH_unbound_red_obs', GSH(b1GPX4=None, state='red') ** Cyto)
Observable('GSH_unbound_ox_obs', GSH(b1GPX4=None, state='ox') ** Cyto)
Observable('GPX4_unbound_inactive_obs', GPX4(b1GSH=None, b2Liperox=None, state='i') ** Cyto)
Observable('GPX4_bound_active_obs', GSH(b1GPX4=1, state='red') ** Cyto % GPX4(b1GSH=1, b2Liperox=None, state='a') ** Cyto)
Observable('Liperox_rad_unbound_obs', Liperox(b2GPX4=None, state='rad') ** Cyto)
Observable('Liperox_red_unbound_obs', Liperox(b2GPX4=None, state='red') ** Cyto)
Observable('GPX4_GSH_Liperox_bound_obs', GSH(b1GPX4=1, state='ox') ** Cyto % GPX4(b1GSH=1, b2Liperox=2, state='a') ** Cyto % Liperox(b2GPX4=2, state='red') ** Cyto)

generate_equations(model)
generate_network(model)
print(model.species)

tspan = np.linspace(0, 1440, 1441) # time span of simulation (start, stop, step)
result = ScipyOdeSimulator(model, tspan=tspan).run() # run solver of model

plt.figure()
plt.plot(tspan, result.observables['Liperox_red_unbound_obs'][:], lw = 1.5, label = 'Liperox_red_unbound_obs') #plot observable
plt.xlabel('Time [min]', fontsize=14)
plt.ylabel('Amount of Lipid alcohol produced [# of molecules]', fontsize=14)
plt.title('Lipid alcohol molecules')
plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# plt.savefig('GSH:GPX4 binding.pdf') # to save figure
plt.show()

plt.figure()
plt.plot(tspan, result.observables['Liperox_rad_unbound_obs'][:], lw = 1.5, label = 'Liperox_rad_unbound_obs') #plot observable
plt.xlabel('Time [min]', fontsize=14)
plt.ylabel('Amount of Liperoxides produced [# of molecules]', fontsize=14)
plt.title('Liperoxide molecules')
plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# plt.savefig('GSH:GPX4 binding.pdf') # to save figure
plt.show()




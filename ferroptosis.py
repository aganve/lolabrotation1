from pylab import *
from pysb.core import * # brings in all of the Python classes needed to define a model
from pysb.bng import *
from pysb.integrate import *
from pysb.macros import catalyze
import matplotlib.pyplot as plt
import numpy as np
from pysb.simulator import ScipyOdeSimulator # Simulating an ODE with SciPy

# instantiate a model
Model()

# Cys2 = oxidized dimer form of the amino acid cysteine
# XC = System XC: cystine/glutamate transporter
# CR = cystine reductase; reduces Cys2 to Cys
# Cys = cysteine
# Glu_Cys = glutamyl-cysteine; intermediate produced in the sequence of reactions converting cysteine to glutathione (GHS)
# GSS = glutathione synthetase; converts Glut_Cys to GSH
# GSH = glutathione; cofactor for GPX4

# declare monomers
Monomer('Cys2', ['b']) # binds cystine reductase (CR) which reduces 1 molecule of cystine to 2 molecules of cysteine
Monomer('XC', ['b']) # no direct binding between Cys2 and XC
Monomer('CR', ['b'])
Monomer('Cys', ['b'])
Monomer('GCL', ['b'])
Monomer('GlutCys', ['b'])
Monomer('GSS', ['b'])
Monomer('GSH', ['bGPX4'])
Monomer('GPX4', ['bGSH'], {'state': ['i', 'a']})
Monomer('Peroxide', ['b'])
Monomer('LipidAlcohol', ['b'])

# Monomer('GCL', ['bCys'])
# Monomer('Glu_Cys', ['bGSS'])
# Monomer('GSS', ['bGlu_Cys'])

# One main 3d compartment and two 2d membranes inside it
# Main 3d compartment is environment and two 2d membranes inside: plasma membrane (PM) and cytosol (Cyto)
Parameter('VEnv', 1)
Parameter('VPM', 1)
Parameter ('VCyto', 1)
Compartment('Env', None, 3, VEnv) # Env doesn't have a parent; it's 3D in volume
Compartment('PM', Env, 2, VPM) # PM's parent is the Env; it's 2D
Compartment('Cyto', PM, 3, VCyto) # Cyto's parent is the PM; it's 3D; proteins can move in the x, y, and z direction, therefore 3D

# input the paramater values
Parameter('kf1', 1.0e-6)
Parameter('kf2', 1.0e-6)
Parameter('kr2', 1.0e-3)
Parameter('kf3', 1.0e-6)
Parameter('kf4', 1.0e-6)
Parameter('kf5', 1.0e-6)
Parameter('kf6', 1.0e-6)
Parameter('kf7', 1.0e-6)
Parameter('kr1', 1.0e-3)

Parameter('kr3', 1.0e-3)

# now input the rules

# transport Cys2 across the PM via XC transporter
Rule('Cys2_transport_via_XC', Cys2() ** Env + XC() ** PM >> Cys2() ** Cyto + XC() ** PM, kf1)
# Rule('Cys2_CR_bind', Cys2(bCR=None) ** Cyto + CR(bCys2=None) ** Cyto | Cys2(bCR=1) ** Cyto % CR(bCys2=1) ** Cyto, kf2, kr2)

# convert Cys2 to Cys via CR catalysis
catalyze(CR, Cys2, Cys, [kf2, kfr2])

# convert Cys to GlutCys via GCL catalysis
catalyze(GCL, Cys, GlutCys,[kf3, kr3])
# Rule('Cys_GSH_conversion', Cys(bCys2=None) ** Cyto | GSH(bGPX4=None) ** Cyto, kf3, kr2)

# convert Glut-Cys to GSH via GSS catalyze
catalyze(GSS, GlutCys, GSH(b=None), [kf4, kr4])

# GPX4 binds GSH, a necessary cofactor to reduce lipid peroxides to lipid alcohols
Rule('GSH_GPX4_bind', GSH(bGPX4=None) ** Cyto + GPX4(bGSH=None, state='i') ** Cyto | GSH(bGPX4=1, state='i') ** Cyto % GPX4(bGSH=1) ** Cyto, kf5, kr5)

# GPX4 is activated when GSH binds
Rule('GPX4_active', GSH(bGPX4=1, state='i') ** Cyto % GPX4(bGSH=1) >> GSH(bGPX4=None) ** Cyto + GPX4(bGSH=None, state='a') ** Cyto, kf6)

# convert lipid peroxide to lipid alcohol via GPX4 catalysis
catalyze(GPX4(bGSH=None, state='a') ** Cyto, Peroxide, LipidAlcohol, [kf7])
# Rule('GPX4_bind_peroxide', GPX4(bGSH=None, state='a') ** Cyto + Peroxide(bGPX4=None) ** Cyto >> GPX4(bGSH=1, state='a') ** Cyto % Peroxide(bGPX4=1) ** Cyto, kf6)


# Initial Conditions
Parameter('Cys2_0', 2300)
Parameter('XC_0', 2000)
Parameter('GPX4_0', 1000)
Initial(Cys2(bXC=None) ** Env, Cys2_0)
Initial(XC(bCys2=None) ** PM, XC_0)
Initial(GPX4(bGSH=None) ** Cyto, GPX4_0)

# Observables
Observable('Cys2_Env_obs', Cys2(bXC=None) ** Env)
Observable('XC_obs', XC(bCys2=None) ** PM)
Observable('Cys2_XC_bind_obs', Cys2(bXC=1) ** PM % XC(bCys2=1) ** PM)
Observable('Cys2_Cyto_obs', Cys2(bXC=None) ** Cyto)
Observable('Cys_obs', Cys(bCys2=None) ** Cyto)
Observable('GSH_obs', GSH(bGPX4=None) ** Cyto)
Observable('GPX4_obs', GPX4(bGSH=None) ** Cyto)
Observable('GSH_GPX4_bind_obs', GSH(bGPX4=1) ** Cyto % GPX4(bGSH=1) ** Cyto)

tspan = np.linspace(0, 1440, 1441) # time span of simulation (start, stop, step)
result = ScipyOdeSimulator(model, tspan=tspan).run() # run solver of model

plt.figure()
plt.plot(tspan, result.observables['GSH_GPX4_bind_obs'][:], lw = 1.5, label = 'GSH_GPX4_bind') #plot observable
plt.xlabel('Time [min]', fontsize=14)
plt.ylabel('Amount of GSH-GPX4 produced [# of molecules]', fontsize=14)
plt.title('GSH-GPX4 complexes')
plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# plt.savefig('GSH:GPX4 binding.pdf') # to save figure
plt.show()

plt.figure()
plt.plot(tspan, result.observables['GPX4_obs'][:], lw = 1.5, label = 'GPX4') #plot observable
plt.xlabel('Time [min]', fontsize=14)
plt.ylabel('Amount of GPX4 utilized [# of molecules]', fontsize=14)
plt.title('Amount of GPX4 utilized over time')
plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# plt.savefig('GPX4.pdf') # to save figure
plt.show()

plt.figure()
plt.plot(tspan, result.observables['GSH_obs'][:], lw = 1.5, label = 'GSH') #plot observable
plt.xlabel('Time [min]', fontsize=14)
plt.ylabel('Amount of GSH utilized [# of molecules]', fontsize=14)
plt.title('Amount of GSH utilized over time')
plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# plt.savefig('GSH.pdf') # to save figure
plt.show()

plt.figure()
plt.plot(tspan, result.observables['Cys_obs'][:], lw = 1.5, label = 'Cys') #plot observable
plt.xlabel('Time [min]', fontsize=14)
plt.ylabel('Amount of Cys produced [# of molecules]', fontsize=14)
plt.title('Amount of Cys produced over time')
plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# plt.savefig('Cys.pdf') # to save figure
plt.show()

plt.figure()
plt.plot(tspan, result.observables['Cys2_Cyto_obs'][:], lw = 1.5, label = 'Cys2') #plot observable
plt.xlabel('Time [min]', fontsize=14)
plt.ylabel('Amount of Cys2 converted [# of molecules]', fontsize=14)
plt.title('Amount of Cys2 consumed over time')
plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# plt.savefig('Cys2.pdf') # to save figure
plt.show()

plt.figure()
plt.plot(tspan, result.observables['Cys2_XC_bind_obs'][:], lw = 1.5, label = 'Cys2_XC_bind') #plot observable
plt.xlabel('Time [min]', fontsize=14)
plt.ylabel('Amount of Cys2-XC complexes [# of molecules]', fontsize=14)
plt.title('Amount of Cys2-XC complexes over time')
plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# plt.savefig('Cys2_XC_bind.pdf') # to save figure
plt.show()

plt.figure()
plt.plot(tspan, result.observables['XC_obs'][:], lw = 1.5, label = 'XC') #plot observable
plt.xlabel('Time [min]', fontsize=14)
plt.ylabel('Amount of XC proteins being bound [# of molecules]', fontsize=14)
plt.title('Amount of XC proteins being bound over time')
plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# plt.savefig('XC.pdf') # to save figure
plt.show()

plt.figure()
plt.plot(tspan, result.observables['Cys2_Env_obs'][:], lw = 1.5, label = 'Cys2_Env') #plot observable
plt.xlabel('Time [min]', fontsize=14)
plt.ylabel('Amount of Cys2 proteins binding to XC [# of molecules]', fontsize=14)
plt.title('Amount of Cys2 protein binding to XC over time')
plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# plt.savefig('Cys2_Env.pdf') # to save figure
plt.show()


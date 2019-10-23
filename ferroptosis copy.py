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
Monomer('GPX4', ['b1GSH', 'b2PUFA_PL', 'state'], {'state': ['i', 'a']})
Monomer('Fe3', ['b1Transferrin','b2STEAP3'])
Monomer('Transferrin', ['b1Fe3', 'b2TFRC'])
Monomer('TFRC', ['b2Transferrin'])
Monomer('STEAP3', ['b2Fe3'])
Monomer('Fe2', ['bH2O2'])
Monomer('H2O2', ['bFe2'])
Monomer('HydroxyRad', ['b1']) # hydroxy radical
Monomer('PUFA_PL', ['b1', 'b2GPX4', 'state'], {'state': ['WT', 'rad', 'peroxy_rad', 'liperox', 'lipid_alcohol']})
Monomer('OXY', ['b1'])
Monomer('FFA', ['b1', 'state'], {'state': ['WT', 'rad']})

# One main 3d compartment
# Main 3d compartment is environment (Env) and three membranes inside: plasma membrane (PM), cytosol (Cyto) and lysosome (Lysos)
Parameter('VEnv', 1)
Parameter('VPM', 1)
Parameter('VCyto', 1)
Parameter('VLysos_PM', 1)
Parameter('VLysos_Cyto', 1)
Compartment('Env', None, 3, VEnv) # Env doesn't have a parent; it's 3D in volume
Compartment('PM', Env, 2, VPM) # PM's parent is the Env; it's 2D in volume
Compartment('Cyto', PM, 3, VCyto) # Cyto's parent is the PM; it's 3D in volume; proteins can move in the x, y, and z direction, therefore 3D
Compartment('Lysos_PM', Cyto, 2, VLysos_PM) # Lysos_PM's parent is the Cyto; it's 2D in volume
Compartment('Lysos_Cyto', Lysos_PM, 3, VLysos_Cyto)

# input the parameter values
Parameter('kf1', 1.0e-2) # kf is a slower reaction than kr
Parameter('kr1', 1.0e-3)
Parameter('kf2', 1.0e-2)
Parameter('kr2', 1.0e-3) # kr is a faster reaction than kf
Parameter('kc2', 1.0e-1) # catalysis occurs more quickly than forward and reverse reactions
Parameter('kf3', 1.0e-2)
Parameter('kr3', 1.0e-3)
Parameter('kc3', 1.0e-1)
Parameter('kf4', 1.0e-2)
Parameter('kr4', 1.0e-3)
Parameter('kc4', 1.0e-1)
Parameter('kf5', 1.0e-2)
Parameter('kr5', 1.0e-3)
Parameter('kf6', 1.0e-2)
Parameter('kr6', 1.0e-3)
Parameter('kf7', 1.0e-2)
Parameter('kf8', 1.0e-2)
Parameter('kf9', 1.0e-2)
Parameter('kf10', 1.0e-2)
Parameter('kf11', 1.0e-2)
Parameter('kf12', 1.0e-2)
Parameter('kf13', 1.0e-2)
Parameter('kf14', 1.0e-2)
Parameter('kf15', 1.0e-2)
Parameter('kf16', 1.0e-2)
Parameter('kf17', 1.0e-2)
Parameter('kf18', 1.0e-2)
Parameter('kf19', 1.0e-2)
Parameter('kf20', 1.0e-2)
Parameter('kf21', 1.0e-2)
Parameter('kf22', 1.0e-2)
Parameter('kf23', 1.0e-2)
Parameter('kr7', 1.0e-3)
Parameter('kr8', 1.0e-3)
Parameter('kr9', 1.0e-3)
Parameter('kr10', 1.0e-3)
Parameter('kr11', 1.0e-3)
Parameter('kr12', 1.0e-3)
Parameter('kr13', 1.0e-3)
Parameter('kr14', 1.0e-3)
Parameter('kr15', 1.0e-3)
Parameter('kr16', 1.0e-3)
Parameter('kr17', 1.0e-3)
Parameter('kr18', 1.0e-3)
Parameter('kr19', 1.0e-3)
Parameter('kr20', 1.0e-3)
Parameter('kr21', 1.0e-3)
Parameter('kr22', 1.0e-3)
Parameter('kr23', 1.0e-3)
Parameter('kc13', 1.0e-1)
Parameter('kc15', 1.0e-1)

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
Rule('GSH_binds_GPX4_activated', GSH(b1GPX4=None, state='red') ** Cyto + GPX4(b1GSH=None, b2PUFA_PL=None, state='i') ** Cyto | GSH(b1GPX4=1, state='red') ** Cyto % GPX4(b1GSH=1, b2PUFA_PL=None, state='a') ** Cyto, kf5, kr5)

# Now that GPX4 is activated with GSH binding, now we look at another competing pathway that regulates iron metabolism...

# Fe3+ binds transferrin in the Env
Rule('Fe3_binds_Transferrin', Fe3(b1Transferrin=None, b2STEAP3=None) ** Env + Transferrin(b1Fe3=None, b2TFRC=None) ** Env | Fe3(b1Transferrin=1, b2STEAP3=None) ** Env % Transferrin(b1Fe3=1, b2TFRC=None) ** Env, kf8, kr8)

# Transferrin-Fe complex binds transferrin receptor (TFRC) in the PM
Rule('Transferrin_Fe_complex_binds_TFRC', Fe3(b1Transferrin=1, b2STEAP3=None) ** Env % Transferrin(b1Fe3=1, b2TFRC=None) ** Env + TFRC(b2Transferrin=None) ** PM | Fe3(b1Transferrin=1, b2STEAP3=None) ** PM % Transferrin(b1Fe3=1, b2TFRC=2) ** PM % TFRC(b2Transferrin=2) ** PM, kf9, kr9) # can you have 3 molecules binding?

# TFRC/transferrin-iron complex is endocytosed from the PM into the Lysos_Cyto
equilibrate(Fe3(b1Transferrin=1, b2STEAP3=None) ** PM % Transferrin(b1Fe3=1, b2TFRC=2) ** PM % TFRC(b2Transferrin=2) ** PM, Fe3(b1Transferrin=1, b2STEAP3=None) ** Lysos_Cyto % Transferrin(b1Fe3=1, b2TFRC=2) ** Lysos_Cyto % TFRC(b2Transferrin=2) ** Lysos_Cyto, [kf10, kr10])

# TFRC/transferrin-iron complex is localized to the lysosome
# equilibrate(Fe3(b1Transferrin=1, b2STEAP3=None) ** Cyto % Transferrin(b1Fe3=1, b2TFRC=1) ** Cyto % TFRC(b2Transferrin=1) ** Cyto, Fe3(b1Transferrin=1, b2STEAP3=None) ** Lysos % Transferrin(b1Fe3=1, b2TFRC=1) ** Lysos % TFRC(b2Transferrin=1) ** Lysos, [kf11, kr11])

# TFRC/transferrin-iron complex dissociates due to acidic conditions of the lysosome
Rule('TFRC_Transferrin_Fe3_complex_dissociates', Fe3(b1Transferrin=1, b2STEAP3=None) ** Lysos_Cyto % Transferrin(b1Fe3=1, b2TFRC=2) ** Lysos_Cyto % TFRC(b2Transferrin=2) ** Lysos_Cyto | Fe3(b1Transferrin=None, b2STEAP3=None) ** Lysos_Cyto + Transferrin(b1Fe3=None, b2TFRC=None) ** Lysos_Cyto + TFRC(b2Transferrin=None) ** Lysos_Cyto, kf12, kr12)

# Fe3+ is oxidized to Fe2+ via STEAP3 catalysis
catalyze(STEAP3() ** Lysos_Cyto, 'b2Fe3', Fe3() ** Lysos_Cyto, 'b2STEAP3', Fe2(bH2O2=None) ** Lysos_Cyto, [kf13, kr13, kc13])

# Fe2+ is exported from the lysosome (Lyso) to the cytosol (Cyto)
equilibrate(Fe2(bH2O2=None) ** Lysos_Cyto, Fe2(bH2O2=None) ** Cyto, [kf14, kr14])

# Now the Fenton Reaction proceeds, where we produce hydroxy radicals by reactin Fe2+ with H2OH. Hydroxy radicals react with lipids to form lipid radicals!

# Fe2+ reacts with hydrogen peroxide (H202) via the Fenton reaction to form an hydroxyl radical (HydroxyRad) in the cytosol (Cyto)
catalyze(Fe2() ** Cyto, 'bH2O2', H2O2() ** Cyto, 'bFe2', HydroxyRad(b1=None) ** Cyto, [kf15, kr15, kc15])

# Hydroxyl radical (HydroxyRad) is translocated from the cytosol (Cyto) to the plasma membrane (PM)
equilibrate(HydroxyRad(b1=None) ** Cyto, HydroxyRad(b1=None) ** PM, [kf16, kr16])

# Now the process ensues where lipid peroxides form but are subsequently reduced to lipid alcohols via GPX4 action...

# Hydroxyl radical (HydroxyRad) binds PL_PUFA (phospholipid-polyunsaturated fatty acid) at the plasma membrane (PM)
Rule("Hydroxyrad_PL_PUFA_binds", HydroxyRad(b1=None) ** PM + PUFA_PL(b1=None, b2GPX4=None, state="WT") ** PM | HydroxyRad(b1=1) ** PM % PUFA_PL(b1=1, b2GPX4=None, state="WT") ** PM, kf17, kr17)

# HydroxyRad-PL_PUFA complex dissociates at the plasma membrane (PM) and PL_PUFA becomes a lipid radical -- PL_PUFA_radical
Rule("PL_PUFA_radical_formation", HydroxyRad(b1=1) ** PM % PUFA_PL(b1=1, b2GPX4=None, state="WT") ** PM | HydroxyRad(b1=None) ** PM + PUFA_PL(b1=None, b2GPX4=None, state="rad") ** PM, kf18, kr18)

# PL_PUFA_radical binds molecular OXY (oxygen) at the PM
Rule("PL_PUFA_OXY_bind", PUFA_PL(b1=None, b2GPX4=None, state="rad") ** PM + OXY(b1=None) ** PM | PUFA_PL(b1=1, b2GPX4=None, state="rad") ** PM % OXY(b1=1) ** PM, kf19, kr19)

# PL_PUFA_radical-OXY complex dissociate and PL_PUFA_radical becomes a PL_PUFA_peroxy_radical at the PM
Rule("PL_PUFA_peroxy_radical_formation", PUFA_PL(b1=1, b2GPX4=None, state="rad") ** PM % OXY(b1=1) ** PM | PUFA_PL(b1=None, b2GPX4=None, state="peroxy_rad") ** PM + OXY(b1=None) ** PM, kf20, kr20)

# PL_PUFA_peroxy_radical binds a FFA at the PM
Rule('PL_PUFA_peroxy_radical_FFA_WT_bond', PUFA_PL(b1=None, b2GPX4=None, state="peroxy_rad") ** PM + FFA(b1=None, state="WT") ** PM | PUFA_PL(b1=1, b2GPX4=None, state="peroxy_rad") ** PM % FFA(b1=1, state="WT") ** PM, kf21, kr21)

# PL_PUFA_peroxy_radical-FFA complex dissociates at the PM; PL_PUFA_peroxy_radical becomes a lipid peroxide (Liperox) and the FFA is converted to a lipid radical
Rule('PL_PUFA_liperox_FFA_radical_formation', PUFA_PL(b1=1, b2GPX4=None, state="peroxy_rad") ** PM % FFA(b1=1, state="WT") ** PM | PUFA_PL(b1=None, b2GPX4=None, state="liperox") ** PM + FFA(b1=None, state="rad") ** PM, kf22, kr22)

# GPX4-GSH complex moves from the Cyto to the PM
equilibrate(GSH(b1GPX4=1, state='red') ** Cyto % GPX4(b1GSH=1, b2PUFA_PL=None, state='a') ** Cyto, GSH(b1GPX4=1, state='red') ** PM % GPX4(b1GSH=1, b2PUFA_PL=None, state='a') ** PM, [kf23, kr23])

# GPX4-GSH complex binds lipid peroxide; GSH is oxidized to GSSG in the reaction and liperox is reduced to a lipid alcohol
Rule('Liperox_binds', GSH(b1GPX4=1, state='red') ** PM % GPX4(b1GSH=1, b2PUFA_PL=None, state='a') ** PM + PUFA_PL(b2GPX4=None, state="liperox") ** PM | GSH(b1GPX4=1, state='ox') ** PM % GPX4(b1GSH=1, b2PUFA_PL=2, state='a') ** PM % PUFA_PL(b2GPX4=2, state="lipid_alcohol") ** PM, kf6, kr6)

# GPX4-GSH-peroxide complex dissociates
Rule('GPX4_GSH_lipid_alcohol_dissociate', GSH(b1GPX4=1, state='ox') ** PM % GPX4(b1GSH=1, b2PUFA_PL=2, state='a') ** PM % PUFA_PL(b2GPX4=2, state="lipid_alcohol") ** PM >> GSH(b1GPX4=None, state='ox') ** PM + GPX4(b1GSH=None, b2PUFA_PL=None, state='i') ** PM + PUFA_PL(b2GPX4=None, state="lipid_alcohol") ** PM, kf7)

# FFA reacts with PL_PUFA_peroxy_radical to form a lipid peroxide (Liperox) at the PM
# catalyze(FFA() ** PM,'b3PUFA_PL', PUFA_PL(b1H2O2=None, b2OXY=None, state="peroxy_rad") ** PM,'b3FFA', Liperox(b2GPX4=None, state='ox'), [kf21, kr21, kc21]) # CHECK THIS AGAIN, CATALYSIS RXN?
# how do you show FFA bind to PUFA_PL when PUFA_PL Is the peroxy radical state?

# Initial Conditions
Parameter('Cys2_0', 2300)
Parameter('CR_0', 2300)
Parameter('GCL_0', 2300)
Parameter('GSS_0', 2300)
Parameter('GPX4_0', 1000)
Parameter('Fe3_0', 2300)
Parameter('Transferrin_0', 2300)
Parameter('TFRC_0', 2300)
Parameter('STEAP3_0', 2300)
Parameter('H2O2_0', 2300)
Parameter('PUFA_PL_0', 2300)
Parameter('OXY_0', 2300)
Parameter('FFA_0', 2300)
Initial(Cys2(bCR=None) ** Env, Cys2_0)
Initial(CR(bCys2=None) ** Cyto, CR_0)
Initial(GCL(bCys=None) ** Cyto, GCL_0)
Initial(GSS(bGlutCys=None) ** Cyto, GSS_0)
Initial(GPX4(b1GSH=None, b2PUFA_PL=None, state='i') ** Cyto, GPX4_0)
Initial(Fe3(b1Transferrin=None, b2STEAP3=None) ** Env, Fe3_0)
Initial(Transferrin(b1Fe3=None, b2TFRC=None) ** Env, Transferrin_0)
Initial(TFRC(b2Transferrin=None) ** PM, TFRC_0)
Initial(STEAP3(b2Fe3=None) ** Lysos_Cyto, STEAP3_0)
Initial(H2O2(bFe2=None) ** Cyto, H2O2_0)
Initial(PUFA_PL(b1=None, b2GPX4=None, state="WT") ** PM, PUFA_PL_0)
Initial(OXY(b1=None) ** PM, OXY_0)
Initial(FFA(b1=None, state="WT") ** PM, FFA_0)

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
Observable('GPX4_unbound_inactive_obs', GPX4(b1GSH=None, b2PUFA_PL=None, state='i') ** Cyto)
Observable('GPX4_bound_active_obs', GSH(b1GPX4=1, state='red') ** Cyto % GPX4(b1GSH=1, b2PUFA_PL=None, state='a') ** Cyto)
Observable('Fe3_Env_obs', Fe3(b1Transferrin=None, b2STEAP3=None) ** Env)
Observable('Transferrin_Env_obs', Transferrin(b1Fe3=None, b2TFRC=None) ** Env)
Observable('Fe3_Transferrin_complex_Env_obs', Fe3(b1Transferrin=1, b2STEAP3=None) ** Env % Transferrin(b1Fe3=1, b2TFRC=None) ** Env)
# Observable('Fe3_Transferrin_complex_Cyto_obs', Fe3(b1Transferrin=1, b2STEAP3=None) ** Cyto % Transferrin(b1Fe3=1, b2TFRC=1) ** Cyto % TFRC(b2Transferrin=1) ** Cyto)
Observable('TFRC_PM_obs', TFRC(b2Transferrin=None) ** PM)
Observable('Fe3_Transferrin_TFRC_complex_PM_obs', Fe3(b1Transferrin=1, b2STEAP3=None) ** PM % Transferrin(b1Fe3=1, b2TFRC=2) ** PM % TFRC(b2Transferrin=2) ** PM)
Observable('Fe3_Transferrin_TFRC_complex_Lysos_Cyto_obs', Fe3(b1Transferrin=1, b2STEAP3=None) ** Lysos_Cyto % Transferrin(b1Fe3=1, b2TFRC=2) ** Lysos_Cyto % TFRC(b2Transferrin=2) ** Lysos_Cyto)
Observable('Fe3_Lysos_Cyto_obs', Fe3(b1Transferrin=None, b2STEAP3=None) ** Lysos_Cyto)
Observable('Transferrin_Lysos_Cyto_obs', Transferrin(b1Fe3=None, b2TFRC=None) ** Lysos_Cyto)
Observable('TFRC_Lysos_Cyto_obs', TFRC(b2Transferrin=None) ** Lysos_Cyto)
Observable('STEAP3_obs', STEAP3(b2Fe3=None) ** Lysos_Cyto)
Observable('Fe2_Lysos_Cyto_obs', Fe2(bH2O2=None) ** Lysos_Cyto)
Observable('Fe2_Cyto_obs', Fe2(bH2O2=None) ** Cyto)
Observable('H2O2_obs', H2O2(bFe2=None) ** Cyto)
Observable('HydroxyRad_Cyto_obs', HydroxyRad(b1=None) ** Cyto)
Observable('HydroxyRad_PM_obs', HydroxyRad(b1=None) ** PM)
Observable('PUFA_PL_WT_obs', PUFA_PL(b1=None, b2GPX4=None, state="WT") ** PM)
Observable('HydroxyRad_PUFA_PL_bound_obs', HydroxyRad(b1=1) ** PM % PUFA_PL(b1=1, b2GPX4=None, state="WT") ** PM)
Observable('PUFA_PL_lipid_radical_obs', PUFA_PL(b1=None, b2GPX4=None, state="rad") ** PM)
Observable('OXY_obs', OXY(b1=None) ** PM)
Observable('PUFA_PL_lipid_radical_OXY_bound_obs', PUFA_PL(b1=1, b2GPX4=None, state="rad") ** PM % OXY(b1=1) ** PM)
Observable('PUFA_PL_peroxy_rad_obs', PUFA_PL(b1=None, b2GPX4=None, state="peroxy_rad") ** PM)
Observable('FFA_WT_obs', FFA(b1=None, state="WT") ** PM)
Observable('PUFA_PL_peroxy_rad_FFA_WT_bound_obs', PUFA_PL(b1=1, b2GPX4=None, state="peroxy_rad") ** PM % FFA(b1=1, state="WT") ** PM)
Observable('PUFA_PL_liperox_obs', PUFA_PL(b1=None, b2GPX4=None, state="liperox") ** PM)
Observable('FFA_rad_obs', FFA(b1=None, state="rad") ** PM)
Observable('GPX4_GSH_Lipid_alcohol_complex_obs', GSH(b1GPX4=1, state='ox') ** Cyto % GPX4(b1GSH=1, b2PUFA_PL=2, state='a') ** Cyto % PUFA_PL(b2GPX4=2, state='lipid_alcohol') ** Cyto)
Observable('PUFA_PL_lipid_alcohol_obs', PUFA_PL(b2GPX4=None, state="lipid_alcohol") ** Cyto)

# generate_equations(model)
# generate_network(model)
# print(model.species)

tspan = np.linspace(0, 1440, 1441) # time span of simulation (start, stop, step)
result = ScipyOdeSimulator(model, tspan=tspan).run() # run solver of model

plt.figure()
plt.plot(tspan, result.observables['PUFA_PL_liperox_obs'][:], lw = 1.5, label = 'PUFA_PL_liperox_obs') #plot observable
plt.xlabel('Time [min]', fontsize=14)
plt.ylabel('Amount of PUFA_PL_liperox_obs produced [# of molecules]', fontsize=14)
plt.title('PUFA_PL_liperox_obs molecules')
plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# plt.savefig('GSH:GPX4 binding.pdf') # to save figure
plt.show()

plt.figure()
plt.plot(tspan, result.observables['FFA_rad_obs'][:], lw = 1.5, label = 'FFA_rad_obs') #plot observable
plt.xlabel('Time [min]', fontsize=14)
plt.ylabel('Amount of PFFA_rad_obs produced [# of molecules]', fontsize=14)
plt.title('FFA_rad_obs molecules')
plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# plt.savefig('GSH:GPX4 binding.pdf') # to save figure
plt.show()

# plt.figure()
# plt.plot(tspan, result.observables['Cys2_Env_obs'][:], lw = 1.5, label = 'Cys2_Env_obs molecules') #plot observable
# plt.xlabel('Time [min]', fontsize=14)
# plt.ylabel('Amount of Cys2_Env_obs produced [# of molecules]', fontsize=14)
# plt.title('Cys2_Env_obs molecules')
# plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# # plt.savefig('GSH:GPX4 binding.pdf') # to save figure
# plt.show()
#
# plt.figure()
# plt.plot(tspan, result.observables['Cys2_Cyto_obs'][:], lw = 1.5, label = 'Cys2_Cyto_obs') #plot observable
# plt.xlabel('Time [min]', fontsize=14)
# plt.ylabel('Amount of Cys2_Cyto_obs produced [# of molecules]', fontsize=14)
# plt.title('Cys2_Cyto_obs molecules')
# plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# # plt.savefig('GSH:GPX4 binding.pdf') # to save figure
# plt.show()
#
# plt.figure()
# plt.plot(tspan, result.observables['CR_obs'][:], lw = 1.5, label = 'CR_obs') #plot observable
# plt.xlabel('Time [min]', fontsize=14)
# plt.ylabel('Amount of CR_obs produced [# of molecules]', fontsize=14)
# plt.title('CR_obs molecules')
# plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# # plt.savefig('GSH:GPX4 binding.pdf') # to save figure
# plt.show()
#
# plt.figure()
# plt.plot(tspan, result.observables['Cys_obs'][:], lw = 1.5, label = 'Cys_obs') #plot observable
# plt.xlabel('Time [min]', fontsize=14)
# plt.ylabel('Amount of Cys_obs produced [# of molecules]', fontsize=14)
# plt.title('Cys_obs_obs molecules')
# plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# # plt.savefig('GSH:GPX4 binding.pdf') # to save figure
# plt.show()
#
# plt.figure()
# plt.plot(tspan, result.observables['GCL_obs'][:], lw = 1.5, label = 'GCL_obs') #plot observable
# plt.xlabel('Time [min]', fontsize=14)
# plt.ylabel('Amount of GCL_obs produced [# of molecules]', fontsize=14)
# plt.title('GCL_obs molecules')
# plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# # plt.savefig('GSH:GPX4 binding.pdf') # to save figure
# plt.show()
#
# plt.figure()
# plt.plot(tspan, result.observables['Glut_Cys_obs'][:], lw = 1.5, label = 'Glut_Cys_obs') #plot observable
# plt.xlabel('Time [min]', fontsize=14)
# plt.ylabel('Amount of Glut_Cys_obs produced [# of molecules]', fontsize=14)
# plt.title('Glut_Cys_obs molecules')
# plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# # plt.savefig('GSH:GPX4 binding.pdf') # to save figure
# plt.show()
#
# plt.figure()
# plt.plot(tspan, result.observables['GSS_obs'][:], lw = 1.5, label = 'GSS_obs') #plot observable
# plt.xlabel('Time [min]', fontsize=14)
# plt.ylabel('Amount of GSS_obs produced [# of molecules]', fontsize=14)
# plt.title('GSS_obs molecules')
# plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# # plt.savefig('GSH:GPX4 binding.pdf') # to save figure
# plt.show()
#
# plt.figure()
# plt.plot(tspan, result.observables['GSH_unbound_red_obs'][:], lw = 1.5, label = 'GSH_unbound_red_obs') #plot observable
# plt.xlabel('Time [min]', fontsize=14)
# plt.ylabel('Amount of GSH_unbound_red_obs produced [# of molecules]', fontsize=14)
# plt.title('GSH_unbound_red_obs molecules')
# plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# # plt.savefig('GSH:GPX4 binding.pdf') # to save figure
# plt.show()
#
# plt.figure()
# plt.plot(tspan, result.observables['GSH_unbound_ox_obs'][:], lw = 1.5, label = 'GSH_unbound_ox_obs') #plot observable
# plt.xlabel('Time [min]', fontsize=14)
# plt.ylabel('Amount of GSH_unbound_ox_obs produced [# of molecules]', fontsize=14)
# plt.title('GSH_unbound_ox_obs molecules')
# plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# # plt.savefig('GSH:GPX4 binding.pdf') # to save figure
# plt.show()
#
# plt.figure()
# plt.plot(tspan, result.observables['GPX4_unbound_inactive_obs'][:], lw = 1.5, label = 'GPX4_unbound_inactive_obs') #plot observable
# plt.xlabel('Time [min]', fontsize=14)
# plt.ylabel('Amount of GPX4_unbound_inactive_obs produced [# of molecules]', fontsize=14)
# plt.title('GPX4_unbound_inactive_obs molecules')
# plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# # plt.savefig('GSH:GPX4 binding.pdf') # to save figure
# plt.show()
#
# plt.figure()
# plt.plot(tspan, result.observables['GPX4_bound_active_obs'][:], lw = 1.5, label = 'GPX4_bound_active_obs') #plot observable
# plt.xlabel('Time [min]', fontsize=14)
# plt.ylabel('Amount of GPX4_bound_active_obs produced [# of molecules]', fontsize=14)
# plt.title('GPX4_bound_active_obs molecules')
# plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# # plt.savefig('GSH:GPX4 binding.pdf') # to save figure
# plt.show()
#
# plt.figure()
# plt.plot(tspan, result.observables['Fe3_Env_obs'][:], lw = 1.5, label = 'Fe3_Env_obs') #plot observable
# plt.xlabel('Time [min]', fontsize=14)
# plt.ylabel('Amount of Fe3_Env_obs produced [# of molecules]', fontsize=14)
# plt.title('Fe3_Env_obs molecules')
# plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# # plt.savefig('GSH:GPX4 binding.pdf') # to save figure
# plt.show()
#
# plt.figure()
# plt.plot(tspan, result.observables['Fe3_Env_obs'][:], lw = 1.5, label = 'Fe3_Env_obs') #plot observable
# plt.xlabel('Time [min]', fontsize=14)
# plt.ylabel('Amount of Fe3_Env_obs produced [# of molecules]', fontsize=14)
# plt.title('Fe3_Env_obs molecules')
# plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# # plt.savefig('GSH:GPX4 binding.pdf') # to save figure
# plt.show()
#
# plt.figure()
# plt.plot(tspan, result.observables['Transferrin_Env_obs'][:], lw = 1.5, label = 'Fe3_Env_obs') #plot observable
# plt.xlabel('Time [min]', fontsize=14)
# plt.ylabel('Amount of Fe3_Env_obs produced [# of molecules]', fontsize=14)
# plt.title('Fe3_Env_obs molecules')
# plt.legend(loc = 'best', fontsize = 12) # add legend to plot
# # plt.savefig('GSH:GPX4 binding.pdf') # to save figure
# plt.show()
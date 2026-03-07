# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: cantera_env
#     language: python
#     name: cantera_env
# ---

# %%
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

# %%
full_path = '/projects/westgroup/nora/Code/'
gas_RMG = ct.Solution(full_path+'projects/rebasing_PFAS/models/PFAS/insights_from_Weber/PFOA+Air/fix_reg_spec_in_core/lower_tolerance/forbidden_group_birad_recomb/chemkin/copies/copy_chem0271.cti')
gas_Weber = ct.Solution(full_path+'projects/rebasing_PFAS/models/PFAS/insights_from_Weber/chemkin_for_Weber_2025_12_12/PFAS_O2_H2O_N2O_publish.cti')

# %%
# fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))


#species order = CF2, H, CF, HF

for spec_index_RMG, spec_index_Web in zip([16, 13, 18, 5], [5, 39, 41, 3]):
    
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))

    print(gas_RMG.species()[spec_index_RMG], gas_Weber.species()[spec_index_Web])
    H_RMG = [gas_RMG.species()[spec_index_RMG].thermo.h(temp)/4184/1000 for temp in range(200,1200)]
    Cp_RMG = [gas_RMG.species()[spec_index_RMG].thermo.cp(temp)/4184/1000 for temp in range(200,1200)]
    S_RMG = [gas_RMG.species()[spec_index_RMG].thermo.s(temp)/4184/1000 for temp in range(200,1200)]


    H_Web = [gas_Weber.species()[spec_index_Web].thermo.h(temp)/4184/1000 for temp in range(200,1200)]
    Cp_Web = [gas_Weber.species()[spec_index_Web].thermo.cp(temp)/4184/1000 for temp in range(200,1200)]
    S_Web = [gas_Weber.species()[spec_index_Web].thermo.s(temp)/4184/1000 for temp in range(200,1200)]

    
    # axes[0].plot(range(200,1200), H_RMG, c='red')
    # axes[0].plot(range(200,1200), H_Web, c='blue')
    H_Web_900_1000 = []
    for val, temp in zip(H_Web, range(200,1200)):
        if temp in range(900,1000):
             H_Web_900_1000.append(val)

    H_RMG_900_1000 = []
    for val, temp in zip(H_RMG, range(200,1200)):
        if temp in range(900,1000):
             H_RMG_900_1000.append(val)                
                
                
    axes[0].plot(range(900,1000), H_Web_900_1000, c='red')
    axes[0].plot(range(900,1000), H_RMG_900_1000, c='blue')
    
    axes[1].plot(range(200,1200), Cp_RMG, c='red')
    axes[1].plot(range(200,1200), Cp_Web, c='blue')
    
    if gas_Weber.species()[spec_index_Web].name=='H':
        print(Cp_Web[0], Cp_Web[-1])
        print(Cp_RMG[0], Cp_RMG[-1]) 
        
    #     axes[1].set_yticks(np.arange(20000, 55000e-7, 2500e-7))
    
    axes[2].plot(range(200,1200), S_RMG, c='red')
    axes[2].plot(range(200,1200), S_Web, c='blue')
    
    axes[0].set_title(gas_Weber.species()[spec_index_Web])

# %% jupyter={"outputs_hidden": true}
#46 CF2 + H <=> CF + HF
#4 CF2(15) + H(12) <=> CF(17) + HF(4)

for rxn_index_RMG, rxn_index_Web in zip([4], [46]):
    


# %%
gas_RMG.reactions()[4]

# %%
gas_RMG.TPX = 650, ct.one_atm, {"PFOA(1)": 4.02e-4, #402 ppm of PFOA,
                            "H2O(3)": 7.50e-4, #750 ppm of H2O(g) 
                            "O2(2)": 0.20975808, #Total of trace species: 4.02e-4 + 7.50e-4 = 1.152e-3, Remaining fraction for air: 1 - 1.152e-3 = 0.998848, 21% O2, 79% N2
                            "N2": 0.78908992,
                            }


print(gas_RMG.forward_rate_constants[4])
print(gas_RMG.net_production_rates[4])

# %%
gas_Weber.TPX = 650, ct.one_atm, {"PFOA": 4.02e-4, #402 ppm of PFOA,
                            "H2O": 7.50e-4, #750 ppm of H2O(g) 
                            "O2": 0.20975808, #Total of trace species: 4.02e-4 + 7.50e-4 = 1.152e-3, Remaining fraction for air: 1 - 1.152e-3 = 0.998848, 21% O2, 79% N2
                            "N2": 0.78908992,
                            }


print(gas_Weber.forward_rate_constants[46])
print(gas_Weber.net_production_rates[46])

# %%
dir(gas_Weber.species()[spec_index_Web].thermo)

# %%

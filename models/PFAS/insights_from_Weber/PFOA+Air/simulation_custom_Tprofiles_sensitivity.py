import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
import csv
from scipy.interpolate import interp1d
import json

def custom_gaussian(x, min_val, max_val, mean, std_dev):
    """
    Generates a Gaussian function with specified minimum and maximum values.

    Args:
        x (numpy.ndarray or float): The input value(s).
        min_val (float): The minimum value (offset) of the function.
        max_val (float): The desired maximum value of the function.
        mean (float): The mean (center) of the Gaussian peak.
        std_dev (float): The standard deviation (width) of the Gaussian.

    Returns:
        numpy.ndarray or float: The calculated Gaussian value(s).
    """
    amplitude = max_val - min_val
    # The standard Gaussian formula
    y = amplitude * np.exp(-(x - mean)**2 / (2 * std_dev**2)) + min_val
    return y


temperature_profiles = pd.read_csv('temperature_profiles_from_literature.csv')


# using isothermal for entire range
isothermal_all_dict = {}
                 
temperatures_to_customize = range(400,1150,50)
for cust_temp in temperatures_to_customize: 
    
    ##### get custom temperature profile edges (distances 0-25, 35-60)
 
    # Parameters
    overshoot = 100
    min_val = 25 #deg C
    max_val = cust_temp + overshoot #overshoot a little higher so the edges match up with the isothermal middle
    mean = 30 #cm
    std_dev = 7

    # Generate x values
    x_values = np.linspace(0, 60, 30) # Range from -4*sigma to 4*sigma covers most of the curve

    # Calculate the corresponding y values
    y_values_edges = custom_gaussian(x_values, min_val, max_val, mean, std_dev)
    T_profile_edges = interp1d(x_values, y_values_edges, kind='linear', fill_value='extrapolate')

    
    ##### get custom temperature profile for middle (distances 25-35)
    # Parameters
    mean_1, mean_2 = 25, 35 #cm
    std_dev = 5


    # Calculate the corresponding y values
    y_values_1 = custom_gaussian(x_values, min_val, max_val, mean_1, std_dev)
    y_values_2 = custom_gaussian(x_values, min_val, max_val, mean_2, std_dev)
    
    T_profile_1 = interp1d(x_values, y_values_1, kind='linear', fill_value='extrapolate')
    T_profile_2 = interp1d(x_values, y_values_2, kind='linear', fill_value='extrapolate')



    ######### put the y values together and run through the interpolate function to get final temperature profile 
    desired_y_values = []
    for x in x_values: 
        if x<22.5: 
            y = T_profile_edges(x)
        if x>=22.5 and x<=37.5: 
            y = cust_temp #isothermal middle
        if x>37.5: 
            y = T_profile_edges(x)
        desired_y_values.append(y)
        
    final_T_profile = interp1d(x_values, desired_y_values, kind='linear', fill_value='extrapolate')
    
    
    
    # Verify the min and max values (should be close to the defined parameters)
    print(f"Calculated maximum value for {cust_temp}: {np.max([final_T_profile(x) for x in x_values]):.2f}, should be close to {cust_temp}")

    
    #add to master dictionary
    isothermal_all_dict[cust_temp] = [x_values, [final_T_profile(x) for x in x_values], final_T_profile]
    
    

       
full_path = '/projects/westgroup/nora/Code/'


mech_paths = {'after': './with_all_halocarbon_chemistry/maxCarbonAtoms/chemkin/copies/copy_chem0116.cti', 
              'main': './just_PFAS/main/chemkin/copies/copy_chem0108.cti', 
              'just_PFAS': './just_PFAS/maxCarbonAtoms/chemkin/copies/copy_chem0169.cti', 

              'cpn_maxCarbonatoms_seed_spec_in_core': full_path+'projects/rebasing_PFAS/models/PFAS/insights_from_Weber/PFOA+Air/fix_reg_spec_in_core/lower_tolerance/cpn_maxCarbonatoms_seed/chemkin/copies/copy_chem0199.cti', #lower tolerance
              'forbidden_group': full_path+'projects/rebasing_PFAS/models/PFAS/insights_from_Weber/PFOA+Air/fix_reg_spec_in_core/lower_tolerance/forbidden_group/chemkin/copies/copy_chem0226.cti',

              'Weber': full_path+'projects/rebasing_PFAS/models/PFAS/insights_from_Weber/chemkin_for_Weber_2025_12_12/PFAS_O2_H2O_N2O_publish.cti', 
             
"forbidden_group_birad_recomb":
full_path+'projects/rebasing_PFAS/models/PFAS/insights_from_Weber/PFOA+Air/fix_reg_spec_in_core/lower_tolerance/forbidden_group_birad_recomb/chemkin/copies/copy_chem0271.cti',             
            "editing_mechanism": full_path+'projects/rebasing_PFAS/models/PFAS/insights_from_Weber/PFOA+Air/fix_reg_spec_in_core/lower_tolerance/forbidden_group_birad_recomb/chemkin/copies/copy_chem_annotated_edited_test.cti'
             }


def convert_to_Kelvin(temp_C): 
    temp_K = temp_C + 273.15 
    return temp_K


def run_sim(temp, mech, calculate_sensitivities=False):
    
    z_data = np.linspace(0,60,500)
    
    
    gas=ct.Solution(mech_paths[mech])
    
    print(mech_paths[mech])
    
    if mech!='Weber': 
        initial_composition ={"PFOA(1)": 4.02e-4, #402 ppm of PFOA,
                            "H2O(3)": 7.50e-4, #750 ppm of H2O(g) 
                            "O2(2)": 0.20975808, #Total of trace species: 4.02e-4 + 7.50e-4 = 1.152e-3, Remaining fraction for air: 1 - 1.152e-3 = 0.998848, 21% O2, 79% N2
                            "N2": 0.78908992,
                            }
    if mech=='Weber':    
        initial_composition ={"PFOA": 4.02e-4, #402 ppm of PFOA,
                        "H2O": 7.50e-4, #750 ppm of H2O(g) 
                        "O2": 0.20975808, #Total of trace species: 4.02e-4 + 7.50e-4 = 1.152e-3, Remaining fraction for air: 1 - 1.152e-3 = 0.998848, 21% O2, 79% N2
                        "N2": 0.78908992,
                        }

    
    #get the temperature profile 

    T_profile = isothermal_all_dict[temp][2]

        
    gas.TPX = 300, ct.one_atm, initial_composition
    
    int_diam = 0.007 #in meters, experiments have a 7 mm ID in reactor
    radius = int_diam/2
    cross_area = np.pi*radius*radius #cross sectional area [m**2]
    Vdot = 150e-6 / 60.0 #volumetric flow rate, 150 mL/min → m³/s
    v0 = Vdot/cross_area #m/s, velocity in x direction 
    length = 0.6 #in m, 60 cms
    

    n_steps = len(z_data)
    dz = length / n_steps
        
    states = ct.SolutionArray(gas, extra=['z', 'time', 'temp']) #for storing state information as we march through the reactor
    
    sensitivity_data = []        
    time = 0 
    for ind_z, z in enumerate(z_data):
        

        T_C = float(T_profile(z))
        T_K = convert_to_Kelvin(T_C)
        gas.TP = T_K, ct.one_atm
        tau = dz / v0  # residence time in this slice
        #r = ct.IdealGasConstPressureReactor(gas, energy='off')  #'off' = fixed T, because we are manually enforcing T
        r = ct.IdealGasReactor(gas, energy='off')  #'off' = fixed T, because we are manually enforcing T
        sim = ct.ReactorNet([r])
        
        if calculate_sensitivities==True: 
            if ind_z==250: #in center of the reactor, at 30.06012024048096 cm 
            #if z==z_data[-1]: #if we're on the last dz 
                for i in range(gas.n_reactions):
                    r.add_sensitivity_reaction(i) #enable sensitivity with respect to the rates
                    # set the tolerances for the solution and for the sensitivity coefficients
                #print('finished adding sensitivities')
        sim.rtol = 1.0e-6
        sim.atol = 1.0e-15
        sim.rtol_sensitivity = 1.0e-6
        sim.atol_sensitivity = 1.0e-6
        
        try:
            time+=tau
            sim.advance(tau)
            if calculate_sensitivities==True:
                if ind_z == 250: 
                #if z==z_data[-1]: #if we're on the last dz 
                    for i in range(gas.n_reactions):

                        reaction_for_sens = gas.reactions()[i]
                        try: 
                            # sensitivity_to_HF = sim.sensitivity('HF', i) #sensitivity of HF to reaction i 
                            # sensitivity_to_C2F4 = sim.sensitivity('C2F4', i) #C2F4
                            # sensitivity_to_CO2 =  sim.sensitivity('CO2', i) #CO2
                            # sensitivity_to_C2F6 =  sim.sensitivity('C2F6', i)   #C2F6, 
                            # sensitivity_to_CO =  sim.sensitivity('CO', i) #CO
                            # sensitivity_to_CF4 = sim.sensitivity('CF4', i) #CF4
                            # sensitivity_to_COF2 = sim.sensitivity('COF2', i)
                            
                            
                            
                            
                            sensitivity_to_HF = sim.sensitivity('HF(4)', i) #sensitivity of HF to reaction i 
                            sensitivity_to_C2F4 = sim.sensitivity('C2F4(5)', i) #C2F4
                            sensitivity_to_CO2 =  sim.sensitivity('CO2(6)', i) #CO2
                            sensitivity_to_C2F6 =  sim.sensitivity('C2F6(7)', i)   #C2F6, 
                            sensitivity_to_CO =  sim.sensitivity('CO(8)', i) #CO
                            sensitivity_to_CF4 = sim.sensitivity('CF4(9)', i) #CF4
                            sensitivity_to_COF2 = sim.sensitivity('COF2(10)', i)
                            sensitivity_data.append([sensitivity_to_HF, sensitivity_to_C2F4, sensitivity_to_CO2, sensitivity_to_C2F6, sensitivity_to_CO, sensitivity_to_CF4, sensitivity_to_COF2, reaction_for_sens.equation])
                        except Exception as e: 
                            print(e)
                            sensitivity_data.append([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, reaction_for_sens.equation])
                        
       
        except Exception as e:
            print(f"\nIntegration failed at z={z:.3f} m, T={T_K:.1f} K")
            print(f"Species with non-finite values: "
                  f"{[s for s, X in zip(gas.species_names, r.thermo.X) if not np.isfinite(X)]}")
            raise e
 
        
        states.append(r.thermo.state, z=z, time=time, temp=gas.T)
        #gas.TPX = T_K, ct.one_atm, r.thermo.X  # propagate composition only
        gas.X = r.thermo.X # propagate composition only
    return states, gas, sensitivity_data 
    
    
   
    
def run_edited_mech_sim(temp, mech_path, Weber=False, calculate_sensitivities=False):
    
    z_data = np.linspace(0,60,500)
    
    
    gas=ct.Solution(mech_path)
    
    if Weber:    
        initial_composition ={"PFOA": 4.02e-4, #402 ppm of PFOA,
                            "H2O": 7.50e-4, #750 ppm of H2O(g) 
                            "O2": 0.20975808, #Total of trace species: 4.02e-4 + 7.50e-4 = 1.152e-3, Remaining fraction for air: 1 - 1.152e-3 = 0.998848, 21% O2, 79% N2
                            "N2": 0.78908992,
                            }
    
    else: 
        initial_composition ={"PFOA(1)": 4.02e-4, #402 ppm of PFOA,
                            "H2O(3)": 7.50e-4, #750 ppm of H2O(g) 
                            "O2(2)": 0.20975808, #Total of trace species: 4.02e-4 + 7.50e-4 = 1.152e-3, Remaining fraction for air: 1 - 1.152e-3 = 0.998848, 21% O2, 79% N2
                            "N2": 0.78908992,
                            }


    
    #get the temperature profile 
    T_profile = isothermal_all_dict[temp][2]

        
    gas.TPX = 300, ct.one_atm, initial_composition
    
    int_diam = 0.007 #in meters, experiments have a 7 mm ID in reactor
    radius = int_diam/2
    cross_area = np.pi*radius*radius #cross sectional area [m**2]
    Vdot = 150e-6 / 60.0 #volumetric flow rate, 150 mL/min → m³/s
    v0 = Vdot/cross_area #m/s, velocity in x direction 
    length = 0.6 #in m, 60 cms
    

    n_steps = len(z_data)
    dz = length / n_steps
        
    states = ct.SolutionArray(gas, extra=['z', 'time', 'temp']) #for storing state information as we march through the reactor
    
    sensitivity_data = []        
    time = 0 
    for i, z in enumerate(z_data):
        
        T_C = float(T_profile(z))
        T_K = convert_to_Kelvin(T_C)
        

        gas.TP = T_K, ct.one_atm
        tau = dz / v0  # residence time in this slice
        #r = ct.IdealGasConstPressureReactor(gas, energy='off')  #'off' = fixed T, because we are manually enforcing T
        r = ct.IdealGasReactor(gas, energy='off')  #'off' = fixed T, because we are manually enforcing T
        sim = ct.ReactorNet([r])

        if calculate_sensitivities==True: 
            for i in range(gas.n_reactions):
                r.add_sensitivity_reaction(i) #enable sensitivity with respect to the rates
                # set the tolerances for the solution and for the sensitivity coefficients
            #print('finished adding sensitivities')
            sim.rtol = 1.0e-6
            sim.atol = 1.0e-15
            sim.rtol_sensitivity = 1.0e-6
            sim.atol_sensitivity = 1.0e-6
            
        try:
            time+=tau
            sim.advance(tau)
            if calculate_sensitivities==True:
                if z==z_data[-1]: #if we're on the last dz 
                    for i in range(gas.n_reactions):
                        print(i)
                        sensitivity_to_HF = sim.sensitivity('HF(4)', i) #sensitivity of HF to reaction i 
                        reaction_for_sens = gas.reactions()[i]
                        sensitivity_data.append([sensitivity_to_HF,reaction_for_sens.equation])
                        
                        
        except Exception as e:
            print(f"\nIntegration failed at z={z:.3f} m, T={T_K:.1f} K")
            print(f"Species with non-finite values: "
                  f"{[s for s, X in zip(gas.species_names, r.thermo.X) if not np.isfinite(X)]}")
            raise e
 
        
        states.append(r.thermo.state, z=z, time=time, temp=gas.T)
        #gas.TPX = T_K, ct.one_atm, r.thermo.X  # propagate composition only
        gas.X = r.thermo.X # propagate composition only
    return states, gas, sensitivity_data     
    
    
temp = 900
#model_path = '/projects/westgroup/nora/Code/projects/rebasing_PFAS/models/PFAS/insights_from_Weber/PFOA+Air/fix_reg_spec_in_core/lower_tolerance/forbidden_group/chemkin/copies/copy_chem0226.cti'
#saved in file: sens_data_forbidden_group.txt


#model_path = '/projects/westgroup/nora/Code/projects/rebasing_PFAS/models/PFAS/insights_from_Weber/chemkin_for_Weber_2025_12_12/PFAS_O2_H2O_N2O_publish.cti'
model_path = '/projects/westgroup/nora/Code/projects/rebasing_PFAS/models/PFAS/insights_from_Weber/PFOA+Air/fix_reg_spec_in_core/lower_tolerance/forbidden_group_birad_recomb/chemkin/copies/copy_chem_annotated_edited_test.cti'

#saved in file: sens_data_Weber_chatGPT.txt

# projects/rebasing_PFAS/models/PFAS/insights_from_Weber/chemkin_for_Weber_2025_12_12/PFAS_O2_H2O_N2O_publish.cti


#states, gas, sens_data = run_edited_mech_sim(temp, model_path, Weber=True, calculate_sensitivities=True)

states, gas, sens_data = run_sim(temp, "editing_mechanism", calculate_sensitivities=True)

with open('sens_data_fbr_edited_mech_900K_center.txt', 'w') as f: 
    f.write(str(sens_data))

# sens_data.sort(key=lambda x: abs(float(x[0])), reverse=True)

# for (sens, equation) in sens_data: 
#     print(sens, equation)
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
#     display_name: rmg_env
#     language: python
#     name: rmg_env
# ---

# %%
import time 
import datetime
from subprocess import PIPE, run
import cantera as ct
import os
import subprocess
import numpy as np
import pandas as pd
import re
import time

# %%
current_env = os.environ.copy()


#change this folder to where you want to save the model
folder_to_models = '/path/to/where/youwouldliketosavethemodels'

#change to RMG database location
RMG_database = '/path/to/your/RMG/database'


# %% [markdown]
# Functions

# %%
def convert_to_cti(description, chemkin_indices, format_desc=True): 
    #copy the necessary files to the folder
    
    if format_desc==False: 
        new_folder = f'./{description}/chemkin'
        RMG_model_directory = f'{folder_to_models}/{description}/'

    else: 
        new_folder = f'./{description}_rebased/chemkin'
        RMG_model_directory = f'{folder_to_models}/{description}_rebased/'

    # if not os.path.exists(f'{new_folder}/chemkin/converter.py'):
    process = subprocess.run(['scp', './converter.py', new_folder], stdout=PIPE, stderr=PIPE, cwd=folder)
    print(process.stdout)
    print(process.stderr)

    #run the converter.py script with the specific chemkin indices
    print(RMG_model_directory)
    process = subprocess.run(['python', 'converter.py', str(chemkin_indices)], stdout=PIPE, stderr=PIPE, cwd=f'{RMG_model_directory}chemkin/')
    print(process.stdout)
    print(process.stderr)
    print('finished')
    

def test_flamespeed(description, chemkin_indices, loglevel, format_desc=True, save=False): 
    
    if format_desc==False: 
        new_folder = f'./{description}/chemkin'
        RMG_model_directory = f'{folder_to_models}/{description}/'

    else: 
        new_folder = f'./{description}_rebased/chemkin'
        RMG_model_directory = f'{folder_to_models}/{description}_rebased/'

    copies_folder = f'{RMG_model_directory}chemkin/copies'
    
    beginning_chemkin = chemkin_indices[0]
    ending_chemkin = chemkin_indices[-1]
    
    if len(chemkin_indices)==1:
        list_=chemkin_indices
    else: 
        list_ = list(range(beginning_chemkin, ending_chemkin))

    ctis = []
    for chemkin_index in list_: 
        if chemkin_index<10:
            cti = f'copy_chem000{chemkin_index}.cti'
        if 100>chemkin_index>=10:
            cti = f'copy_chem00{chemkin_index}.cti'        
        if chemkin_index>=100:
            cti = f'copy_chem0{chemkin_index}.cti'
        ctis.append(cti)

    for cti in ctis: 
        
        cti_path = f'{copies_folder}/{cti}'

        print(f'*******************Starting cti: {cti} of {description} ***************')
        gas = ct.Solution(cti_path)
        halocarbon = 'CH3F(1)'

        To = 298
        Po = 1e5 # ct.one_atm


        mole_frac_list = list(np.linspace(0.025, 0.25, 50))
        mole_frac_list=[0.125]

        results = {}

        for i in range(len(mole_frac_list)): 
            try: 
                x = mole_frac_list[i]
                string = f'****************************starting new volume fraction: {x} **************************'
                print(string)

                norm_ox = (1-x)*.21
                mole_frac_dict = {halocarbon: x, 'O2(2)':((1-x)*.21), 'N2':((1-x)*0.79)} 
                #print(f'Unnormalized composition dictionary: {mole_frac_dict}')

                #normalize it to O2 
                mole_frac_dict = {halocarbon: (x/norm_ox), 'O2(2)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox } 
                #print(f'Normalizing: {mole_frac_dict}')

                gas.TPX = To, Po, mole_frac_dict
                width = 0.08
                flame = ct.FreeFlame(gas, width=width)
                flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
                flame.max_time_step_count = 3000
                flame.solve(loglevel=loglevel, auto=False)
                Su = flame.velocity[0]
                results[x] = Su
                if save==True: 
                    sltn = flame.to_solution_array()
                    pd = sltn.to_pandas()
                    pd.to_csv(f'data/{x}.csv', index=False)

            except Exception as e: 
                print(f'********************passed volume fraction:{mole_frac_list[i]}, error: {e}*************************************')
                pass

        vol_fracs = list(results.keys())
        flame_speeds = list(results.values())

        print(f"flamespeed: {flame_speeds}")




# %% [markdown]
#

# %%

# %% [markdown]
# # Testing Flamespeeds

# %%
#testing CH3F main

loglevel=0
chemkin_indices=[124] #if you know the numbers of the chemkin, write as integers
description = 'CH3F_main/CH3F_main' #this is the path to the model
test_flamespeed(description, chemkin_indices, loglevel, format_desc=False) #only format_desc==False if '_rebased' not in model name

# %%
#testing CH3F main

loglevel=0
chemkin_indices=[12, 30] #you can also provide a range
description = 'CH3F_main/CH3F_main'
convert_to_cti(description, chemkin_indices, format_desc=False)
test_flamespeed(description, chemkin_indices, loglevel, format_desc=False)

# %%
#can test many RMG models at once
loglevel=0
chemkin_indices=[21] #if you want to do the chem_annotated cti, write at chemkin_indices=['chem_annotated.inp'], give as string
for fam in families_ive_rebased:
    description = f'all_families_on_rebase_edited_by_Nora/one_kinfam_at_a_time/{fam}_rebased'
    convert_to_cti(description, chemkin_indices, format_desc=False)
    test_flamespeed(description, chemkin_indices, loglevel, format_desc=False)

# %%

# %%

# %%

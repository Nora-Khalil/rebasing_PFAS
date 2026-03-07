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
def clean_slate(): 
    "restore database to how it was (no edits to main's kinetics database)"
    process = subprocess.run(['git', 'reset', '--hard'], stdout=PIPE, stderr=PIPE, cwd=RMG_database)
    print('Restoring database to main.')
    print(process.stdout)
    print(process.stderr)
    
def checkout_family(family): 
    """
    Checkout the retrained family from commit 
    908a25d5f354e34ba46d44c1bd47b76cad14b5fc from RMG-database branch PFAS_kinetics_rebased_Jan2025.
    """
    #git_commit = '908a25d5f354e34ba46d44c1bd47b76cad14b5fc'
    git_commit = '4fee96d55dc156935f56a4cd4db68c7ab36cfe85' #one commit up from 908a25d5f35, with updated thermo for CS resonance structures
    family_folder_path = f'input/kinetics/families/{family}/*'
    process = subprocess.run(['git', 'checkout', git_commit, family_folder_path], stdout=PIPE, stderr=PIPE, cwd=RMG_database)
    print(f'Checking out the following family: {family}.')
    print(process.stdout) 
    print(process.stderr)
    
def checkout_family_from_before_rebase(family): 
    """
    Checkout the retrained family from commit 
    2d5ae5a4fe49390214d156d40839a5895aeccd48from RMG-database branch USNCM_blends. 
    Since that was a dropped commit, using edc0f17edb3c05354313c8e6e2e8d5c44ad85ddc from https://github.com/Nora-Khalil/RMG-database/commits/reinstating_dropped_commit_for_good_CH3F_model, where I reapplied the dropped commit. 
    """
    git_commit = 'edc0f17edb3c05354313c8e6e2e8d5c44ad85ddc' #this is the dropped commit from USNCM blends, before Matt's average fix was applied
    family_folder_path = f'input/kinetics/families/{family}/*'
    process = subprocess.run(['git', 'checkout', git_commit, family_folder_path], stdout=PIPE, stderr=PIPE, cwd=RMG_database)
    print(f'Checking out the following family: {family}.')
    print(process.stdout) 
    print(process.stderr)    
    
    
def checkout_multiple_families_from_before_rebase(families): 
    """
    Checkout the retrained family from commit 
    2d5ae5a4fe49390214d156d40839a5895aeccd48from RMG-database branch USNCM_blends.
    Since that was a dropped commit, using edc0f17edb3c05354313c8e6e2e8d5c44ad85ddc from https://github.com/Nora-Khalil/RMG-database/commits/reinstating_dropped_commit_for_good_CH3F_model, where I reapplied the dropped commit. 
    """
    git_commit = 'edc0f17edb3c05354313c8e6e2e8d5c44ad85ddc'
    for family in families: 
        family_folder_path = f'./input/kinetics/families/{family}/*'
        process = subprocess.run(['git', 'checkout', git_commit, family_folder_path], stdout=PIPE, stderr=PIPE, cwd=RMG_database)
    print(f'Checking out the following families: {families}.')
    print(process.stdout) 
    print(process.stderr)    

    
    
def checkout_multiple_families(families): 
    """
    Checkout the retrained family from commit 
    908a25d5f354e34ba46d44c1bd47b76cad14b5fc from RMG-database branch PFAS_kinetics_rebased_Jan2025.
    """
    #git_commit = '908a25d5f354e34ba46d44c1bd47b76cad14b5fc'
    git_commit = '4fee96d55dc156935f56a4cd4db68c7ab36cfe85' #one commit up from 908a25d5f35, with updated thermo for CS resonance structures
    for family in families: 
        family_folder_path = f'./input/kinetics/families/{family}/*'
        process = subprocess.run(['git', 'checkout', git_commit, family_folder_path], stdout=PIPE, stderr=PIPE, cwd=RMG_database)
    print(f'Checking out the following families: {families}.')
    print(process.stdout) 
    print(process.stderr)

def generate_RMG_model(description, input_type): 
    """
    Generates a new RMG model in a new folder. 
    """
    new_folder = f'./{description}_rebased'
    RMG_model_directory = f'{folder_to_models}/{description}_rebased/'
    if not os.path.exists(RMG_model_directory):
        os.mkdir(RMG_model_directory)
    
    #select the input.py needed 
    input_dict = {'no_PFAS_fams_or_new_libs': 'input_no_PFAS_families_or_new_libs.py', #does not use new PFAS fams or new CH2F2/PFAS libs
                  'no_PFAS_fams_or_PFAS_libs': 'input_no_PFAS_families_or_PFAS_libs.py', #does not use new PFAS fams or new PFAS libs (does use new CH2F2 lib)
                  'no_PFAS_fams': 'input_no_PFAS_families.py', #does not use new PFAS fams, but includes new CH2F2/PFAS libraries
                  'PFAS_fams': 'input.py', 
                  
                }
    input_file = input_dict[input_type]
    
    #copy the necessary files to the folder
    process = subprocess.run(['scp', f'./{input_file}', './run.sh', RMG_model_directory], stdout=PIPE, stderr=PIPE, encoding='utf-8', cwd=folder_to_models)
    print(process.stdout)
    print(process.stderr)
    
    #rename the input file so it matches what's in the run.sh
    process = subprocess.run(['mv', input_file, 'input.py'], stdout=PIPE, stderr=PIPE, encoding='utf-8', cwd=RMG_model_directory)
    print(process.stdout)   
    print(process.stderr)
    
    
    #run the RMG job
    process = subprocess.run(['sbatch', 'run.sh'], stdout=PIPE, stderr=PIPE, env=current_env, encoding='utf-8', cwd = RMG_model_directory)
    print(process.stdout)
    print(process.stderr)
    


# %%

#these are the families that I personally added PFAS reactions to
combined_families = ['CO_CF_bond_dissociation',
                 'PFAS_Hydrolysis', 
                 'Lactone_Formation',
                 'CO2_Elimination_From_Carboxylic_Acid',
                 'Perfluoroalkene_Formation',
                 'Enol_Ether_Formation',
    
                 '1+2_Cycloaddition',
                 'XY_Addition_MultipleBond',
                 'R_Addition_MultipleBond',
                 '1,2_Insertion_CO', 
                 '1,3_Insertion_CO2',
                 '1,3_sigmatropic_rearrangement',
                 'R_Addition_COm',
                 'F_Abstraction', 
                  'halocarbene_CO_dimerization'   #need to add in because halogens uses it

]

#the families that me and DF have added reactions to
families_ive_rebased = ['CO_CF_bond_dissociation',
                 'PFAS_Hydrolysis', 
                 'Lactone_Formation',
                 'CO2_Elimination_From_Carboxylic_Acid',
                 'Perfluoroalkene_Formation',
                 'Enol_Ether_Formation',
    
                '1+2_Cycloaddition',
                 'XY_Addition_MultipleBond',
                 'R_Addition_MultipleBond',
                 '1,2_Insertion_CO', 
                 '1,3_Insertion_CO2',
                 '1,3_sigmatropic_rearrangement',
                 'R_Addition_COm',
                 'F_Abstraction',
                                
                 '1,2_Insertion_carbene', #these were ones that David's reactions were in (but not PFAS)
                 'CO_Disproportionation',
                 'Disproportionation', 
                 'H_Abstraction',
                 'R_Recombination',
                 'Singlet_Carbene_Intra_Disproportionation', 
                  'halocarbene_CO_dimerization'   #need to add in because halogens uses it
 
                       ]


#these are families that either (1) I added reactions to, (2) DF added reactions to, or (3) I found there to be a git diff in the reactions.py when compared to main
halocarbon_fams = ['1+2_Cycloaddition',
                 'XY_Addition_MultipleBond',
                 'R_Addition_MultipleBond',
                 '1,2_Insertion_CO', 
                 '1,3_Insertion_CO2',
                 '1,3_sigmatropic_rearrangement',
                 'R_Addition_COm',
                 'F_Abstraction', 
                   
                   
                 '1,2_Insertion_carbene', #these were ones that David's reactions were in (but not PFAS)
                 'CO_Disproportionation',
                 'Disproportionation', 
                 'H_Abstraction',
                 'R_Recombination',
                 'Singlet_Carbene_Intra_Disproportionation',  
                  
                #other git diffs I found
                 '1,2-Birad_to_alkene',
                 '1,2_XY_interchange',
                 'Birad_R_Recombination',
                 'Br_Abstraction',
                   'Cl_Abstraction',
                   'Concerted_Intra_Diels_alder_monocyclic_1,2_shiftH',
                   'Disproportionation-Y',
                   'Disproportionation',
                   'Intra_RH_Add_Endocyclic',
                   'Intra_ene_reaction',
                   'XY_elimination_hydroxyl',
                   'halocarbene_CO_dimerization',
                   'halocarbene_recombination',
                   'intra_H_migration',
                   'intra_halogen_migration',

]

# %%

# %% [markdown]
# # Making the RMG models

# %% [markdown]
# General instructions: 
#
# Build models by putting together different code blocks. 
#
# - **clean_slate()**: this git resets the current RMG database commit. If you want to start on a branch that is essential an up-to-date main, make sure you're on commit eb6010936b0d0b5b6ea4445db2238643bde851d5 from https://github.com/Nora-Khalil/RMG-database/commits/rebased_PFAS_libs before using this function. 
#
# - **checkout_family(family)**: check out an individual family from my rebased (and retrained) RMG database. This has all of my PFAS kinetics changes on top of an up-to-date main, and is pulling from commit 4fee96d55dc156935f56a4cd4db68c7ab36cfe85 from https://github.com/Nora-Khalil/RMG-database/commits/rebased_Feb2025_and_including_PFLactone_carbene_ether_thermo. 
#
# - **checkout_family_from_before_rebase(family)**: checkout an individual family from the same non-rebased RMG database that I used to produce my good model (on my outdated development branch). Pulling from commit edc0f17edb3c05354313c8e6e2e8d5c44ad85ddc from https://github.com/Nora-Khalil/RMG-database/commits/reinstating_dropped_commit_for_good_CH3F_model. This is a reapplication of the original dropped commit used to make the good model, 2d5ae5a4fe49390214d156d40839a5895aeccd48from RMG-database branch USNCM_blends.
#
#
# - **checkout_multiple_families(families)**: check out multiple families from my rebased (and retrained) RMG database. This has all of my PFAS kinetics changes on top of an up-to-date main, and is pulling from commit 4fee96d55dc156935f56a4cd4db68c7ab36cfe85 from https://github.com/Nora-Khalil/RMG-database/commits/rebased_Feb2025_and_including_PFLactone_carbene_ether_thermo. 
#
# - **checkout_multiple_families_from_before_rebase(families)**: checkout multiple families from the same non-rebased RMG database that I used to produce my good model (on my outdated development branch). Pulling from commit edc0f17edb3c05354313c8e6e2e8d5c44ad85ddc from https://github.com/Nora-Khalil/RMG-database/commits/reinstating_dropped_commit_for_good_CH3F_model. This is a reapplication of the original dropped commit used to make the good model, 2d5ae5a4fe49390214d156d40839a5895aeccd48from RMG-database branch USNCM_blends.
#
# - **generate_RMG_model(description, input_type)**: copy the input.py, run.sh to the necessary folder, then submit it via sbatch. 

# %% [markdown]
#

# %% [markdown]
# Recreate Bad Model Database

# %%
#Recreate Bad Model Database

#make sure you're starting on eb6010936b0d0b5b6ea4445db2238643bde851d5, on https://github.com/Nora-Khalil/RMG-database/commits/rebased_PFAS_libs
clean_slate() #reset the commit in case you made any changes. This starts you on a database that looks similar to main 
checkout_multiple_families(families_ive_rebased) #checking out multiple from the RMG-database branch where I rebased and retrained my families, https://github.com/Nora-Khalil/RMG-database/commits/rebased_Feb2025_and_including_PFLactone_carbene_ether_thermo
description = './model_folder_name' #name the model
input_type = 'PFAS_fams' #the input.py will specify that we want to use PFAS families
generate_RMG_model(description, input_type)

# %%

# %%
#if you want to make multiple jobs that have one family checked out from my rebased database. 
input_type = 'PFAS_fams'

directory_for_multiple_models = './changing_one_fam_at_a_time' #make this directory if it doesn't exist
count=0
for fam in families_ive_rebased:
    
    if not os.path.exists(f'{directory_for_multiple_models}/{fam}_rebased/chemkin'):
        print(f'############## Starting family: {fam} #################')

        #reset the database
        clean_slate()

        #checkout out the new training data
        checkout_family(fam) #coming from rebased database, https://github.com/Nora-Khalil/RMG-database/commits/rebased_Feb2025_and_including_PFLactone_carbene_ether_thermo

        #generate the RMG model
        generate_RMG_model(f'{directory_for_multiple_models}/{fam}', input_type)
        count +=1
print(f'Count: {count}')

# %%

# %%
#for deleting folders if you need it
for fam in _: #iterate over a list of families
    print(f'Deleting useless folders for {fam}.')
    
    process = subprocess.run(['rm', '-rf', f'./{fam}_rebased/rms/',  f'./{fam}_rebased/pdep/',  f'./{fam}_rebased/species', f'./{fam}_rebased/rms/solver/'], stdout=PIPE, stderr=PIPE)


# %%

# %%

# %%
#regenerating_model_at_different_commits/checking_out_families_from_before_rebase/git_checkout_one_family_at_a_time_plus_all_ive_rebased

input_type = 'PFAS_fams'

RMG_families_path = '/home/khalil.nor/Code/RMG-database/input/kinetics/families'

families_from_good_model = [#'ANL_Brown_pdep',
 #'FFCM1(-)',
 #'halogens_pdep',
 'Birad_R_Recombination',
 'Birad_recombination',
 'R_Addition_MultipleBond',
 '1+2_Cycloaddition',
 'R_Recombination',
 'Intra_Disproportionation',
 'intra_H_migration',
 '1,2_Insertion_CO',
 '1,3_Insertion_CO2',
 'H_Abstraction',
 'Disproportionation',
 'Intra_R_Add_Endocyclic',
 'R_Addition_COm',
 'intra_halogen_migration',
 '1,4_Linear_birad_scission',
 'halocarbene_recombination',
 'XY_Addition_MultipleBond',
 'CO_Disproportionation',
 'Singlet_Val6_to_triplet',
 'Cyclic_Ether_Formation',
 'Intra_R_Add_Exocyclic',
 'F_Abstraction',
 'HO2_Elimination_from_PeroxyRadical',
 '1,2_Insertion_carbene',
 'XY_elimination_hydroxyl',
 'Retroene',
 '1,4_Cyclic_birad_scission',
 'Singlet_Carbene_Intra_Disproportionation',
 'PFAS_Hydrolysis',
 '1,3_Insertion_ROR',
 'Ketoenol',
 'Disproportionation-Y',
 '1,3_sigmatropic_rearrangement',
 '1,2-Birad_to_alkene',
 'Enol_Ether_Formation_charge_separated',
 '1,2_shiftC',
 '2+2_cycloaddition',
 'intra_OH_migration',
 '1,2_XY_interchange']

pdep_rxn_families = ["Birad_R_Recombination",

#"Birad_Recombination",

"R_Addition_MultipleBond",

"1+2_Cycloaddition",

"R_Recombination",

"Intra_Disproportionation",

"Intra_H_migration",

"1,2_Insertion_CO",

"1,3_Insertion_CO2"
                    ]

#other_fams_in_the_model = [fam for fam in families_from_USNCM_model if fam not in families_ive_rebased]


for fam in pdep_rxn_families:
    
        if os.path.exists(f'./checking_out_families_from_before_rebase/checkout_Birad_recomb_plus_one_other_family/{fam}_rebased/chemkin'):
    
            print(f'skipping {fam}')
            continue
    
        print(f'############## Starting family: {fam} #################')

        #reset the database
        clean_slate()
        
        #make sure to git checkout the recommended.py
        process = subprocess.run(['git', 'checkout', '908a25d5f354e34ba46d44c1bd47b76cad14b5fc', 'input/kinetics/families/recommended.py'], stdout=PIPE, stderr=PIPE, cwd='/home/khalil.nor/Code/RMG-database')
        print(process.stdout) 
        print(process.stderr)
        
        #checkout all of the families that I've rebased
        #checkout_multiple_families_from_before_rebase(families_ive_rebased)
        checkout_multiple_families(families_ive_rebased)
        
        #check out the birad recombination from before 
        checkout_family_from_before_rebase('Birad_recombination')
        
        #checkout out the new training data
        checkout_family_from_before_rebase(fam) #one of the families of the pdep networks

        #generate the RMG model
        generate_RMG_model(f'./checking_out_families_from_before_rebase/checkout_Birad_recomb_plus_one_other_family/{fam}', input_type)
        
        #wait until a chemkin file is generated 
        minutes = 5
        time.sleep(minutes*60)
        if os.path.exists(f'./checking_out_families_from_before_rebase/checkout_Birad_recomb_plus_one_other_family/{fam}_rebased/chemkin'):
            print('Job has started!')
            continue
        else: 
            print('Have to wait a little longer')
            minutes = 3
            time.sleep(minutes*60) 
            if os.path.exists(f'./checking_out_families_from_before_rebase/checkout_Birad_recomb_plus_one_other_family/{fam}_rebased/chemkin'):
                print('Ok lets go!')
                continue
            else: 
                print('Theres an issue...')
                break

        
        

# %%

# %% [markdown]
# # Scancel Jobs 

# %%
#for canceling jobs if you screw up   
jobs = [ #will need to modify this list, I get this from printing "squeue .."
          "48725319     short   S_CH3F khalil.n  R    1:23:08      1 c3018",
          "48725241     short   S_CH3F khalil.n  R    1:28:25      1 c0281",
          "48725193     short   S_CH3F khalil.n  R    1:33:37      1 c0171",
          "48725137     short   S_CH3F khalil.n  R    1:38:49      1 c0225",
          "48725061     short   S_CH3F khalil.n  R    1:43:49      1 c0229",
          "48724998     short   S_CH3F khalil.n  R    1:48:48      1 c0172",
          "48724925     short   S_CH3F khalil.n  R    1:53:49      1 c0232",
          "48724842     short   S_CH3F khalil.n  R    1:58:49      1 c0231",
          "48724764     short   S_CH3F khalil.n  R    2:04:00      1 c0314",
          "48724674     short   S_CH3F khalil.n  R    2:09:17      1 c0223",
          "48724555     short   S_CH3F khalil.n  R    2:14:35      1 c0206",
          "48724427     short   S_CH3F khalil.n  R    2:21:59      1 c0216",
          "48724356     short   S_CH3F khalil.n  R    2:27:52      1 c0222",
          "48724181     short   S_CH3F khalil.n  R    2:32:53      1 c2203",
          "48724085     short   S_CH3F khalil.n  R    2:37:53      1 c3022",
          "48724018     short   S_CH3F khalil.n  R    2:43:10      1 c0319",
          "48722782     short   S_CH3F khalil.n  R    2:53:55      1 c0315",
          "48722033     short   S_CH3F khalil.n  R    2:59:00      1 c0279",
          "48721518     short   S_CH3F khalil.n  R    3:04:56      1 c0243",
          "48721453     short   S_CH3F khalil.n  R    3:09:57      1 c0230"

]

job_ids = []
for line in jobs: 
    match = re.search('^([0-9]+)     short   S_CH3F', line).group(1)
    print(match)
    job_ids.append(match)


# %%
# #copy the necessary files to the folder
for job_id in job_ids: 
    process = subprocess.run(['scancel', job_id], stdout=PIPE, stderr=PIPE, encoding='utf-8')
    print(process.stdout)
    print(process.stderr)

# %%

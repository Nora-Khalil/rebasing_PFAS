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
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import os
import re

# %%
logs = [log for log in os.listdir('.') if 'fixing_reg_' in log and 'error' not in log]
error_logs = [log for log in os.listdir('.') if 'fixing_reg_' in log and 'error' in log]

# %%
PFAS_relevant_families = ['1+2_Cycloaddition',
'1,2_Insertion_CO',
'1,2_Insertion_carbene',
'1,3_Insertion_CO2',
'1,3_sigmatropic_rearrangement',
'CO2_Elimination_From_Carboxylic_Acid',
'CO_CF_bond_dissociation',
'CO_Disproportionation',
'Disproportionation',
'Enol_Ether_Formation',
'F_Abstraction',

'H_Abstraction',
'Lactone_Formation',
'PFAS_Hydrolysis',
'Perfluoroalkene_Formation',
'R_Addition_COm',
'R_Addition_MultipleBond',
'R_Recombination',
'Singlet_Carbene_Intra_Disproportionation',
'XY_Addition_MultipleBond',
'halocarbene_recombination'
]

# %%
for ind, fam in enumerate(PFAS_relevant_families):
    print(ind+1, fam)

# %%
for log in error_logs:
    cascade = False
    # print('\n')
    try: 
        with open(log, 'r') as f: 
            #family_index = int(re.search('fixing_reg_([0-9]+)', log).group(1))
            family_index = int(re.search('fixing_reg_error_([0-9]+)', log).group(1))

            #print(log, family_index)
            
            lines = f.readlines()
            
            if lines==[]:
                print(f'not completed yet: {log}')
                continue

            #last_line = lines[-1]
            #print(last_line)
            
            # if 'all done!' not in last_line:
            #     print(log)
            #     print(PFAS_relevant_families[family_index-1])
            
            for line in lines:
                if 'pruning tree' in line: 
                    cascade = True
            if cascade == True: 
                print('used cascading algo')
                print(PFAS_relevant_families[family_index-1])
    except IndexError: 
        print(f'passing {log}')


# %%
family_folders = [folder for folder in os.listdir('.') if '_retraining' in folder]

for folder in family_folders:
    
    try: 
        with open(f'./{folder}/groups.py', 'r') as f: 
            group_info = f.read()  

        #is there a group that include Si?
        if any([x in group_info for x in [',Si','Si,']]):

            #check if Si is in the training reactions 
            with open(f'./{folder}/training/dictionary.txt', 'r') as f: 
                dictionary_contents = f.read()
                if ' Si ' not in dictionary_contents: 
                    print(folder)
                    print('Group has Si but no Si in training reactions.')
    except FileNotFoundError: 
        continue
            
            
    

# %%

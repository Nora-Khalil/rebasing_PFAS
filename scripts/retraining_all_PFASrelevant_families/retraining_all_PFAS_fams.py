import rmgpy
import numpy as np
from rmgpy.molecule.molecule import *
from rmgpy.species import *
from rmgpy.chemkin import *
from rmgpy.data.rmg import RMGDatabase
from IPython.display import display
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.rmg.react import react
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.data.rmg import get_db
from rmgpy.molecule.group import Group
from rmgpy.kinetics.arrhenius import ArrheniusBM
from rmgpy import settings
import time
from arkane.main import get_git_commit
import matplotlib.pyplot as plt
import matplotlib
import os
import sys



print(f"You are on the following\nRMG-database: {get_git_commit(settings['database.directory'])}\nRMG-Py: {get_git_commit('/home/khalil.nor/Code/RMG-Py/rmgpy')}")


thermo_libs = [
'C1_C2_Fluorine', #adding Siddha's as first most trusted because this is the thermo library that Franklin used
'PFCA_thermo',
#'NCSU_C2_C8_PFAS', #adding Westmoreland's thermo as the second most trusted
'primaryThermoLibrary',
'Fluorine',
'FFCM1(-)',
'halogens',
'CHOF_G4',
'CHOCl_G4',
'CHOBr_G4',
'CHOFCl_G4',
'CHOFBr_G4',
'CHOFClBr_G4',
'DFT_QCI_thermo',
'2-BTP_G4',
'thermo_DFT_CCSDTF12_BAC',
'SulfurHaynes'
]

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
'halocarbene_recombination',
'Retroene'
]


family_index = int(sys.argv[-1])

kin_families = [PFAS_relevant_families[family_index]]

#make the folder to save to 
os.makedirs(f'./{PFAS_relevant_families[family_index]}_retraining', exist_ok=True)

database = RMGDatabase()
database.load(
    path = settings['database.directory'],
    thermo_libraries = thermo_libs,  # Can add others if necessary
    kinetics_families = kin_families,
    reaction_libraries = [],
    kinetics_depositories = ['training'],
)


family_to_train = PFAS_relevant_families[family_index]
family = database.kinetics.families[family_to_train]

family.clean_tree()
print('cleaned tree')

start = time.time()
family.generate_tree(thermo_database=database.thermo,
                     nprocs=1,
                     new_fraction_threshold_to_reopt_node=0.25,
                     max_batch_size=800,
                     extension_iter_max=2,
                     extension_iter_item_cap=100)

end = time.time()
print(end-start)
print('generated tree')

print(f'Length of family groups: {len(family.groups.entries)}')

start = time.time()
family.check_tree()
end = time.time()
print(end-start)

print('checked tree')

#let's look at the nodes that need regularization 
uncommon_atoms = ['S', 'Li', 'N', 'P', 'F', 'I', 'Cl', 'Br', 'O', 'Si']
nodes_to_watch = {}
print('collecting nodes to watch')
for node_label, node in family.groups.entries.items():
    if all(ua in node.item.to_adjacency_list() for ua in uncommon_atoms):
        nodes_to_watch[node_label] = [node.item.to_adjacency_list()]
        

start = time.time()
family.regularize(thermo_database=database.thermo)
end = time.time()
print(end-start)
print('regularized')


start = time.time()
family.check_tree()
end = time.time()
print(end-start)
print('checking tree')


print('After regularization is performed:')
for node_label, node in family.groups.entries.items():
    if node_label in nodes_to_watch.keys():
        print(node_label)
        print('before regularization')
        print(nodes_to_watch[node_label][0])
        print('\n\n')
        print('after regularization')
        print(node.item.to_adjacency_list())


start = time.time()
templateRxnMap = family.get_reaction_matches(thermo_database=database.thermo,remove_degeneracy=True, get_reverse=True,exact_matches_only=False,fix_labels=True)
end = time.time()
print(end-start)
print('got matches')

print(len(templateRxnMap))


family.clean_tree_rules()

start = time.time()
family.make_bm_rules_from_template_rxn_map(templateRxnMap)#,nprocs=6)
end = time.time()
print(end-start)
print('made bm rules')

start = time.time()
family.check_tree()
end = time.time()
print(end-start)
print('checking tree')


#only save to this folder 
save_here = f'./cascading_algo_fix/{PFAS_relevant_families[family_index]}_retraining/'
family.save(save_here)
print('all done!')



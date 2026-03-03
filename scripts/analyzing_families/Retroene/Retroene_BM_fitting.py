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
#     display_name: rmg_env_Jul2025
#     language: python
#     name: python3
# ---

# %%
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
import matplotlib.pyplot as plt
import matplotlib

# %%

# %%
import os

# %%

# %%
os.system('hostname')

# %%
thermo_libs = [
# 'C1_C2_Fluorine', #adding Siddha's as first most trusted because this is the thermo library that Franklin used
# 'PFCA_thermo',
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

kin_families = ['Retroene']


# %%
database = RMGDatabase()
database.load(
            path = settings['database.directory'],
            thermo_libraries = thermo_libs,
            transport_libraries = [],
            reaction_libraries = [],
            seed_mechanisms = [],#['BurkeH2O2inN2','ERC-FoundationFuelv0.9'],
            kinetics_families = kin_families,
            kinetics_depositories = ['training'],
            #frequenciesLibraries = self.statmechLibraries,
            depository = False, # Don't bother loading the depository information, as we don't use it
        )


# %%
database.kinetics.families

# %%
family_to_train = "Retroene"
family = database.kinetics.families[family_to_train]

# %%
start = time.time()
templateRxnMap = family.get_reaction_matches(thermo_database=database.thermo,remove_degeneracy=True,
                                             get_reverse=True,exact_matches_only=False,fix_labels=True)
end = time.time()
print(end-start)

# %%
len(templateRxnMap)

# %%
#"Root_N-4R!H->C",
for node_label, node in family.groups.entries.items():
    if node_label == "Root_N-4R!H->C":
        #family.simple_regularization(node, templateRxnMap, test=True)
        node_to_look_at=node

# %%
node_to_look_at.data

# %%
depo = family.get_training_depository()

# %%
for rxn in templateRxnMap["Root_N-4R!H->C"]:
    display(rxn)
    print(rxn.reactants[0].to_adjacency_list())

# %%
print(node_to_look_at.parent.item.to_adjacency_list())

# %%
print(node_to_look_at.item.to_adjacency_list())

# %%
for node_label, node in family.groups.entries.items():
    if str(node.parent)=='Root':
        print(node_label)

# %%
for node_label, node in family.groups.entries.items():
    if node_label == "Root_4R!H->C":
        print(node.item.to_adjacency_list())

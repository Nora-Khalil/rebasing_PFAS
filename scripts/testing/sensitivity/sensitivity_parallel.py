import cantera as ct
from PIL import Image
from subprocess import run
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import re 
import os
import numpy as np
import csv
import pandas as pd
import sys


###################### get the model #######################################


blend_name = sys.argv[1]


directory = './chemkin/copies/copy_chem_annotated_17.cti'


#create the species concentrations
halocarbon = 'CH3F(1)'
x = 0.13520408163265307 #this is the vf where CH3F reached a max

norm_ox = (1-x)*.21
species_con =  {halocarbon: (x/norm_ox), 'O2(2)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox } 



########### make the gas, specify the TPX

gas = ct.Solution(directory)
gas.TPX = 298, ct.one_atm, species_con

########### create the FreeFlame object in Cantera, solve the flame

f = ct.FreeFlame(gas)
f.solve()
    
########## save a flux diagram 
# diagram = ct.ReactionPathDiagram(f.gas,'H')
# diagram.show_details = True
# diagram.flow_type = 'NetFlow'
# diagram.scale=-1

# diagram.title = 'Reaction path diagram following H'
# diagram.label_threshold = 0.05

# dot_file = f'rxnpath_{blend_name}_rmg_details_2nd_time.dot'
# img_file = f'rxnpath_{blend_name}_rmg_details_2nd_time.png'
# img_path = Path.cwd().joinpath(img_file)

# diagram.write_dot(dot_file)
# #print(diagram.get_data())

# print("Wrote graphviz input file to '{0}'.".format(Path.cwd().joinpath(dot_file)))

# run('dot {0} -Tpng -o{1} -Gdpi=200'.format(dot_file, img_file).split())
# print("Wrote graphviz output file to '{0}'.".format(img_path))


########## get sensitivities of each reaction rate on flamespeed
'''
.get_flame_speed_reaction_sensitivities(): Compute the normalized sensitivities of the laminar flame speed S with respect to the reaction rate constants k:

                    s_i =   k_i    dS_u
                           -----  -------
                            S_U    dk_i

'''
# ######## sort the sensitivities by magnitude #######################################

sens = f.get_flame_speed_reaction_sensitivities()

sensitivity = {}
for m in range(gas.n_reactions):
    sensitivity[m] = abs(sens[m])

sorted_sensitivity = dict(sorted(sensitivity.items(), key=lambda item: item[1], reverse=True)) #sort with highest magnitude first

######### revert the sensitivity values back to original sign  ####################

#sorted_sensitivity_list = [[k,sens[k],gas.reaction(k)] for k,v in sorted_sensitivity.items() ]

data = {
    'k_s': [k for k,v in sorted_sensitivity.items()], #this is number of reaction in gas.reactions list
    'sensitivity': [sens[k] for k,v in sorted_sensitivity.items()], #sensitivity
    'cantera equation': [gas.reaction(k).equation for k,v in sorted_sensitivity.items()],
    'cantera products': [gas.reaction(k).products for k,v in sorted_sensitivity.items()],
    'cantera reactants': [gas.reaction(k).reactants for k,v in sorted_sensitivity.items()],
}

df = pd.DataFrame(data)
df.to_csv(f'{blend_name}_RMG_sensitivities.csv', index=False)





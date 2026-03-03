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
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
import csv
from scipy.interpolate import interp1d
import json

# %matplotlib inline

# %%
file = './sens_data_fbr_650K_center.txt'

print('650 K')

with open(file, 'r') as f: 
    data = f.read()

data = eval(data)


HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1][0])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2']


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec=='C2F4':
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break


file = './sens_data_fbr_edited_700K_center.txt'

print('700 K')
with open(file, 'r') as f: 
    data = f.read()

data = eval(data)


HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1][0])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2'
]


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec=='C2F4':
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break

                        

file = './sens_data_fbr_edited_750K_center.txt'

print('750 K')
with open(file, 'r') as f: 
    data = f.read()

data = eval(data)


HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1][0])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2'
]


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec=='C2F4':
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break


# %%

# %%

# %%
#tested mech
file = './sens_data_fbr_edited_mech_650K_center.txt'

lets_see = 'HF'

print('650 K')

with open(file, 'r') as f: 
    data = f.read()

data = eval(data)


HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2']


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec==lets_see:
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break


#tested mech
file = './sens_data_fbr_edited_mech_700K_center.txt'

print('700 K')

with open(file, 'r') as f: 
    data = f.read()

data = eval(data)


HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2']


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec==lets_see:
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break


#tested mech
file = './sens_data_fbr_edited_mech_750K_center.txt'

print('750 K')

with open(file, 'r') as f: 
    data = f.read()

data = eval(data)


HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2']


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec==lets_see:
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break
                
#tested mech
file = './sens_data_fbr_edited_mech_850K_center.txt'

print('850 K')

with open(file, 'r') as f: 
    data = f.read()

data = eval(data)


HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2']


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec==lets_see:
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break                

                
                
#tested mech
file = './sens_data_fbr_edited_mech_900K_center.txt'

print('900 K')

with open(file, 'r') as f: 
    data = f.read()

data = eval(data)


HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2']


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec==lets_see:
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break                 

# %%
#tested mech
file = './sens_data_fbr_edited_mech_650K_center.txt'

lets_see = 'C2F4'

print('650 K')

with open(file, 'r') as f: 
    data = f.read()

data = eval(data)


HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2']


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec==lets_see:
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break


#tested mech
file = './sens_data_fbr_edited_mech_700K_center.txt'

print('700 K')

with open(file, 'r') as f: 
    data = f.read()

data = eval(data)


HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2']


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec==lets_see:
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break


#tested mech
file = './sens_data_fbr_edited_mech_750K_center.txt'

print('750 K')

with open(file, 'r') as f: 
    data = f.read()

data = eval(data)


HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2']


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec==lets_see:
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break
                
#tested mech
file = './sens_data_fbr_edited_mech_850K_center.txt'

print('850 K')

with open(file, 'r') as f: 
    data = f.read()

data = eval(data)


HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2']


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec==lets_see:
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break                

                
                
#tested mech
file = './sens_data_fbr_edited_mech_900K_center.txt'

print('900 K')

with open(file, 'r') as f: 
    data = f.read()

data = eval(data)


HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2']


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec==lets_see:
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break                 

# %%

# %%

# %%
#file = './sens_data_Weber_650K_all_quick.txt'
file = './sens_data_Weber_650K_center.txt'
print('650 K')
with open(file, 'r') as f: 
    data = f.read()
    
data = eval(data)

HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1][0])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)


rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2'
]


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec=='HF':
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break

# %%
file = './sens_data_Weber_700K_center.txt'
print('700 K')
with open(file, 'r') as f: 
    data = f.read()
    
data = eval(data)

HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1][0])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)


rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2'
]


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec=='C2F4':
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break

# %%
file = './sens_data_Weber_750K_center.txt'
print('750 K')
with open(file, 'r') as f: 
    data = f.read()
    
data = eval(data)

HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1][0])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)


rankings = [HF_sens_ranks,
C2F4_sens_ranks,
CO2_sens_ranks,
C2F6_sens_ranks,
CO_sens_ranks ,
CF4_sens_ranks,
COF2_sens_ranks]

species_of_focus = ['HF',
'C2F4',
'CO2',
'C2F6',
'CO',
'CF4',
'COF2'
]


for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
    if spec=='C2F4':
        n=0
        print(f'*******************{spec}******************')
        for sens_data in ranking:
            if n<15: 
                print(sens_data[ind], sens_data[-1])
                n+=1
            else:
                break

# %%

# %%

# %%
#westmoreland thermo
temps = [600, 650, 700, 750]
lets_see = 'HF'

for temp in temps: 
    file = f'sens_data_westmoreland_thermo_{temp}K_238species.txt'


    print(f'{temp} K')

    with open(file, 'r') as f: 
        data = f.read()

    data = eval(data)


    HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
    C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
    CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
    C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
    CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
    CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
    COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

    rankings = [HF_sens_ranks,
    C2F4_sens_ranks,
    CO2_sens_ranks,
    C2F6_sens_ranks,
    CO_sens_ranks ,
    CF4_sens_ranks,
    COF2_sens_ranks]

    species_of_focus = ['HF',
    'C2F4',
    'CO2',
    'C2F6',
    'CO',
    'CF4',
    'COF2']


    for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
        if spec==lets_see:
            n=0
            print(f'*******************{spec}******************')
            for sens_data in ranking:
                if n<15: 
                    print(sens_data[ind], sens_data[-1])
                    n+=1
                else:
                    break

# %%
#westmoreland thermo edited 
temps = [600, 650, 700]
lets_see = 'C2F4'

for temp in temps: 
    file = f'sens_data_westmoreland_thermo_{temp}K_edited.txt'


    print(f'{temp} K')

    with open(file, 'r') as f: 
        data = f.read()

    data = eval(data)


    HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
    C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
    CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
    C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
    CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
    CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
    COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

    rankings = [HF_sens_ranks,
    C2F4_sens_ranks,
    CO2_sens_ranks,
    C2F6_sens_ranks,
    CO_sens_ranks ,
    CF4_sens_ranks,
    COF2_sens_ranks]

    species_of_focus = ['HF',
    'C2F4',
    'CO2',
    'C2F6',
    'CO',
    'CF4',
    'COF2']


    for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
        if spec==lets_see:
            n=0
            print(f'*******************{spec}******************')
            for sens_data in ranking:
                if n<15: 
                    print(sens_data[ind], sens_data[-1])
                    n+=1
                else:
                    break

# %%
#weber
temps = [650, 700, 750, 800, 850, 900, 950]
lets_see = 'HF'

for temp in temps: 
    file = f'./sens_data_Weber_{temp}K_center.txt'


    print(f'{temp} K')

    with open(file, 'r') as f: 
        data = f.read()

    data = eval(data)

    try: 
        HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
        C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1][0])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
        CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
        C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
        CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
        CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
        COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)
    except TypeError as e: 
        HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
        C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
        CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
        C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
        CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
        CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
        COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)
        
    rankings = [HF_sens_ranks,
    C2F4_sens_ranks,
    CO2_sens_ranks,
    C2F6_sens_ranks,
    CO_sens_ranks ,
    CF4_sens_ranks,
    COF2_sens_ranks]

    species_of_focus = ['HF',
    'C2F4',
    'CO2',
    'C2F6',
    'CO',
    'CF4',
    'COF2']


    for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
        if spec==lets_see:
            n=0
            print(f'*******************{spec}******************')
            for sens_data in ranking:
                if n<15: 
                    print(sens_data[ind], sens_data[-1])
                    n+=1
                else:
                    break

# %%

# %%

# %%

# %%
#"pdep_fix_RAMBlib"
temps = [650, 700, 750, 800, 850, 900, 950]
lets_see = 'C2F4'

for temp in temps: 
    file = f'./sens_data_pdep_fix_RAMBlib_{temp}K_edited.txt'

    print(f'{temp} K')

    with open(file, 'r') as f: 
        data = f.read()

    data = eval(data)


    HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
    C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
    CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
    C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
    CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
    CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
    COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)

    rankings = [HF_sens_ranks,
    C2F4_sens_ranks,
    CO2_sens_ranks,
    C2F6_sens_ranks,
    CO_sens_ranks ,
    CF4_sens_ranks,
    COF2_sens_ranks]

    species_of_focus = ['HF',
    'C2F4',
    'CO2',
    'C2F6',
    'CO',
    'CF4',
    'COF2']


    for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
        if spec==lets_see:
            n=0
            print(f'*******************{spec}******************')
            for sens_data in ranking:
                if n<15: 
                    print(sens_data[ind], sens_data[-1])
                    n+=1
                else:
                    break

# %%
#Brown Diflouromethane Seed
temps = [600, 650, 700, 750, 800, 850, 900, 950]

lets_see = 'C2F4'

for temp in temps: 
    file = f'./sens_data_FM_{temp}K.txt'

    print(f'{temp} K')

    with open(file, 'r') as f: 
        data = f.read()

    data = eval(data)


    HF_sens_ranks = sorted(data, key=lambda x: abs(float(x[0])), reverse=True)
    C2F4_sens_ranks = sorted(data, key=lambda x: abs(float(x[1])), reverse=True) #this was mistakenly appended to the sensitivity data list as a tuple
    CO2_sens_ranks = sorted(data, key=lambda x: abs(float(x[2])), reverse=True)
    C2F6_sens_ranks = sorted(data, key=lambda x: abs(float(x[3])), reverse=True)
    CO_sens_ranks = sorted(data, key=lambda x: abs(float(x[4])), reverse=True)
    CF4_sens_ranks = sorted(data, key=lambda x: abs(float(x[5])), reverse=True)
    COF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[6])), reverse=True)
    #CF2_sens_ranks = sorted(data, key=lambda x: abs(float(x[7])), reverse=True)

    rankings = [HF_sens_ranks,
    C2F4_sens_ranks,
    CO2_sens_ranks,
    C2F6_sens_ranks,
    CO_sens_ranks ,
    CF4_sens_ranks,
    COF2_sens_ranks, 
    #CF2_sens_ranks]
               ]
    species_of_focus = ['HF',
    'C2F4',
    'CO2',
    'C2F6',
    'CO',
    'CF4',
    'COF2', 
    #'CF2'
                       ]

    #looking at sensitivities one at a time

    # for ind, (spec, ranking) in enumerate(zip(species_of_focus, rankings)):
        
        #looking at sensitivities one at a time
        # if spec==lets_see:
        #     n=0
        #     print(f'*******************{spec}******************')
        #     for sens_data in ranking:
        #         if n<15: 
        #             print(sens_data[ind], sens_data[-1])
        #             n+=1
        #         else:
        #             break
    combined_sens = []    
            
    for sens_data in HF_sens_ranks: #just using HF because it species of interest only dictates what its sorted by, but info is same in all
        final_sens = 0
        for ind, spec in enumerate(species_of_focus):
            final_sens+=abs(sens_data[ind])
        combined_sens.append((final_sens, sens_data[-1]))
        
    n=0
    print(f'*******************{temp} K******************')
    for (final_sensitivity, eq) in combined_sens:
        if n<15: 
            print(final_sensitivity, eq)
            n+=1
        else:
            break        


# %%

import cantera as ct
import numpy as np
import pandas as pd
import os 
import re
from subprocess import getoutput
import sys
import csv


print("Running Cantera Version: " + str(ct.__version__))

indices_to_convert = eval(sys.argv[-1]) #turns string into a list

beginning_chemkin = indices_to_convert[0]
ending_chemkin = indices_to_convert[-1]

if type(beginning_chemkin)==int: 
    if len(indices_to_convert)==1 and beginning_chemkin<10: 
        list_of_blends = [f'chem000{beginning_chemkin}.inp']
    if len(indices_to_convert)==1 and beginning_chemkin>=10 and beginning_chemkin<100: 
        print(beginning_chemkin)
        list_of_blends = [f'chem00{beginning_chemkin}.inp'] 
    if len(indices_to_convert)==1 and beginning_chemkin>=100: 
        list_of_blends = [f'chem0{beginning_chemkin}.inp']
    if len(indices_to_convert)>1: #there's a list and we want a range 
        list_ = list(range(beginning_chemkin, ending_chemkin))
        list_of_blends = [f'chem000{num}.inp' for num in list_ if num<10] + [f'chem00{num}.inp' for num in list_ if 100>num>=10] + [f'chem0{num}.inp' for num in list_ if num>=100]
if type(beginning_chemkin)==str: 
    list_of_blends = indices_to_convert

    

fullpath = '.'

#file_name = 'chem0036.inp'

def convert(full_path_to_chemkin_folder, file_name):
    ############### copies chemkin files to dups folder #############################
        
    os.chdir(full_path_to_chemkin_folder)


    #os.system('source activate ct_env') #if on local, this is cantera 2.6 beta
    os.system('source activate cantera_env') #if on discovery, will be cantera 2.5


    #copy folders so i dont screw up the original, and change into the new directory with copies
    os.makedirs('copies', exist_ok=True)

    #copy chem.inp file into dups folder, and will then convert this copy into .cti
    command = f'scp {file_name} copies/copy_{file_name}'
    os.system(command)
    os.system('scp tran.dat copies/tran.dat')

    #now look in the dups folder
    os.chdir('./copies')


    ############################# converts the dup_chem.inp files to .cti files #############################

    x = 0
    while x == 0: 
    #this is a string of the output when I try to convert this file to a cti file. Will probably produce an error
        output = getoutput( f'ck2cti --input=copy_{file_name} --transport=tran.dat') 
        #if command passed without an error
        if re.search('PASSED',output):
            print('**************************command passed, converting to cti***************************')
            x += 1     
            
        #if command generated the Duplicate error
        else:
            print('*******************************needs some work****************************************')
            if re.search('Encountered\sunmarked\sduplicate\sreaction',output):
                match = re.search('See\slines\s([0-9]+)\sand\s([0-9]+)\sof\sthe\sinput\sfile',output)
                #capture the line numbers with the duplicate reactions
                line_numbers = [int(match.group(1)), int(match.group(2))]
                print(f'Unmarked duplicates on lines {match.group(1)} and {match.group(2)}')
                print('Editing chemkin file to allow conversion to .cti')
                #write the lines of the chemkin input file to a list so that I can insert the "DUPLICATE" statement
                with open(f'./copy_{file_name}','r') as f:
                    data = f.readlines()
                    print('reading in files')
                    print(data[line_numbers[0]-1], data[line_numbers[1]-1])

                #start editing the .inp file below

                #'adjustments' will make sure that, even when I add an element in 'data', my index will still be correct
                adjustments = [0,1]
                for i,adjust in zip(line_numbers,adjustments): 
                    start = i+adjust-1
                    count = 0 
                    while count == 0: 
                        #if you don't see a blank line after the duplicated reaction line, keep going until you do
                        if not re.search('^\n', data[start]): 
                            print('no match')
                            print(start)
                            print(data[start])
                            start += 1
                        #when we get to the blank line after the reaction block, insert "DUPLICATE" and stop the loop for this line number
                        else: 
                            print('there is a match')
                            data.insert(start,'DUPLICATE')
                            print('inserted DUPLICATE')
                            count = 1 
                #now overwrite the input file with the change 
                with open(f'copy_{file_name}','w+') as f: 
                    for l in data: 
                        f.write(l)
                        x==0
                print('overwrote file')

                

            #if you get this second Duplicate error          
            elif re.search('Undeclared duplicate reactions detected',output):
                
                print(f'Undeclared duplicate reactions error. See: {output}')
                duplicate_reactions = re.findall('Undeclared duplicate reactions detected:\nReaction [0-9]+: (.+)\nReaction [0-9]+: (.+)\n', output)
                edited_duplicate_reactions=[] #manually add in the (+M)s, change numbers, and remove the spaces
                for dup_rxn in list(duplicate_reactions[0]): #duplicate_reactions is a list of one 2-item tuple
                    if '2 ' in dup_rxn: #manually edit so its A+A=B+B instead of 2A=2B
                        splitup = dup_rxn.split()
                        edited_splitup = []
                        for ind, string in enumerate(splitup):
                            if string=='2':
                                edited_splitup.append(splitup[ind+1])
                                edited_splitup.append('+')
                            else: 
                                edited_splitup.append(string)
                        dup_rxn = ' '.join(edited_splitup)
                    dup_rxn_edited = dup_rxn.replace(' <=>',' (+M)<=>', 1).replace(' ','')
                    dup_rxn_edited = dup_rxn_edited.replace(dup_rxn_edited[-11:], dup_rxn_edited[-11:]+'(+M)')
                    edited_duplicate_reactions.append(dup_rxn_edited)
                
                #now edit the chemkin file 
                with open(f'./copy_{file_name}','r') as f:
                    data = f.readlines()
                                        
                #edit the data lines
                new_data_for_chemkin = []
                for ind, data_line in enumerate(data): 
                    
                    
                    if data_line=='\n' and ind>11:
                        for unmarked_dup in edited_duplicate_reactions: 
                            if unmarked_dup in data[ind-10]: 
                                new_data_for_chemkin.append('DUPLICATE\n')
                                new_data_for_chemkin.append(data_line)
                                print('added in flag')
                            else: #add in the \n to the new data for chemkin
                                new_data_for_chemkin.append(data_line)
                                
                                
                    
                    else: 
                        new_data_for_chemkin.append(data_line)

                #now overwrite the input file with the change 
                with open(f'copy_{file_name}','w+') as f: 
                    for l in new_data_for_chemkin: 
                        f.write(l)
                        x==0
             
          
            #if the command generated an error that is not the Duplicate error
            else:
                #if code ever gets to here, just cry
                print('There is another error, see Output')
                print(re.search('Encountered\sunmarked\sduplicate\sreaction',output))
                print(output)
                x += 1
    os.chdir('../')

for file_name in list_of_blends: 
    convert(f'{fullpath}', file_name)

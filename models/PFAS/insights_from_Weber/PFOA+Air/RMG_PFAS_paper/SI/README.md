## Summary of SI Materials

- **RMG_main.yaml**: mechanism file for the RMG main model (plotted in Figure 17 and 18 in dashed blue).
- **RMG_present_work.yaml**: mechanism file for RMG model of the present work (PW) (plotted in Figures 13, 14, 17 and 18 in solid blue).
- **RMG_input_file.py**: example RMG input file for model generation.
- **NCSU_2025_dictionary.txt**: adjacency-list formatted dictionary defining the species reported to be relevant to PFCA gas-phase thermal oxidation by [Ram et al. (2024)](https://pubs.acs.org/doi/full/10.1021/acs.jpca.3c06937).
- **cantera_pfr_simulation.ipynb**: IPython notebook with Cantera plug flow reactor simulation for PFOA thermal degradation. Includes plotting scripts to recreate Figures 1, 13, 14, 16, 17, 18 (temperature profiles, concentration ratio graphs, PFOA decomposition along reactor).
- **RMG_present_work_results.json**: PFR simulation results generated using the RMG PW model.
- **weber_mechanism_results.json**: PFR simulation results generated using the Weber mechanism. 
- **weber_experiment_data** folder contains the following:
    - CSV files containing data points manually extracted by the authors of the present work from the concentration ratio graphs reported by [Weber et al.](https://pubs.acs.org/doi/full/10.1021/acs.jpca.4c01909)
    - **PFAS_O2_H2O_N2O_publish.yaml**: YAML-formatted mechanism file of the kinetic model for thermal decomposition of PFOA presented by [Weber et al.](https://pubs.acs.org/doi/full/10.1021/acs.jpca.4c01909)
      


## Where to find the database presented in this work? 

The modified RMG-database (and corresponding RMG-Py) presented in this work exists as branches on the official RMG github repository. 
- [RMG-database](https://github.com/ReactionMechanismGenerator/RMG-database/tree/PFAS_JPCA_paper)
- [RMG-Py](https://github.com/ReactionMechanismGenerator/RMG-Py/tree/PFAS_JPCA_paper)

## Which Github commits should be used to reproduce the RMG models included in this work? 

The PW model (included in Figures 13, 14, 17 and 18 in solid blue) can be generated using the following Github commits.
- RMG-database: commit [55117ddae9820a7d6142676100c27743f86b2e95](https://github.com/ReactionMechanismGenerator/RMG-database/tree/55117ddae9820a7d6142676100c27743f86b2e95)
- RMG-Py: commit [c30f89a5961eae5e530aef7480835accdda2f9ef](https://github.com/ReactionMechanismGenerator/RMG-Py/tree/c30f89a5961eae5e530aef7480835accdda2f9ef]

The main model (included in Figure 17 and 18 in dashed blue) can be generated with the following Github commits: 

- RMG-database: commit [0ce1c5004e66567190113ffa14e30e87e16d5415](https://github.com/ReactionMechanismGenerator/RMG-database/tree/0ce1c5004e66567190113ffa14e30e87e16d5415)

- RMG-Py: commit [c30f89a5961eae5e530aef7480835accdda2f9ef](https://github.com/ReactionMechanismGenerator/RMG-Py/tree/c30f89a5961eae5e530aef7480835accdda2f9ef]

## Instructions for generating RMG models 

RMG models should be generated according to the instructions on the official RMG [website](https://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/running.html). 

The **RMG_input_file.py** includes the RMG input parameters required to generate the models in this work. This file can be used as is to generate the RMG main model, or can be modified (by uncommenting specific blocks of code) to generate the PW model. 

## Instructions for generating concentration ratio plots

The **cantera_pfr_simulation.ipynb** should be used for PFR simulation. Within the second cell of section entitled "PFR Simulation", the yaml_path variable should be edited depending on which mechanism the user would like to simulate. Additionally, the Vdot variable (representing PFR inlet flow) should be edited depending on the desired initial conditions of the simulation. 

PFR simulation results can be generated within this script and used for plotting directly, or saved as .json files to be plotted later. 









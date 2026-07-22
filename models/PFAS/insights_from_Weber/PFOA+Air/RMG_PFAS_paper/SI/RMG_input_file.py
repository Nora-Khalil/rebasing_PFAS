thermolibs = [

###### uncomment to generate the RMG model of the present work ########
#'NCSU_2025_revised',
#'C1_C2_Fluorine',
#'PFCA_thermo',
    
'primaryThermoLibrary',
'FFCM1(-)',
'halogens',
'CHOF_G4',
'CHOCl_G4',
'CHOBr_G4',
'CHOFCl_G4',
'CHOFBr_G4',
'CHOFClBr_G4',
'DFT_QCI_thermo',
'Fluorine',
'2-BTP_G4',
'thermo_DFT_CCSDTF12_BAC',
'SulfurHaynes',
]


database(
thermoLibraries = thermolibs,
reactionLibraries = ['FFCM1(-)'], 
seedMechanisms = [],
kineticsDepositories = ['training'],
kineticsFamilies = ['default', 'halogens','Disproportionation-Y'],
frequenciesLibraries = ['halogens_G4'],
kineticsEstimator = 'rate rules',
)

###### uncomment to generate the RMG model of the present work 
#database(
#thermoLibraries = thermolibs,
#reactionLibraries = ['FFCM1(-)', 'halogens_pdep'],
#seedMechanisms = [ "ANL_Brown_pdep", "ANL_Brown_pdep_2", 'Brown_difluoromethane_seed', "PFOA_unimolecular_decomp"],
#kineticsDepositories = ['training'],
#kineticsFamilies = ['default','PFAS', 'halogens','Disproportionation-Y'],
#frequenciesLibraries = ['halogens_G4'],
#kineticsEstimator = 'rate rules',
#)


#in N2 bath gas 
species(
    label = 'N2',
    reactive = False,
    structure = SMILES('N#N')
)

coreSpeciesFile(
    species_dictionary_file="./NCSU_2025_dictionary.txt"
)

simpleReactor(
        temperature=[(400,'K'),(1600,'K')],
        pressure=(1.0,'bar'),
        nSims=12,
        initialMoleFractions={
        "PFOA": 4.02e-4, #402 ppm of PFOA
        "H2O": 7.50e-4, #750 ppm of H2O(g) 
        "O2": 0.20975808, #Total of trace species: 4.02e-4 + 7.50e-4 = 1.152e-3, Remaining fraction for air: 1 - 1.152e-3 = 0.998848, 21% O2, 79% N2
        "N2": 0.78908992,
        },
        terminationConversion={
        "PFOA": 0.9999,
        },
        terminationTime=(4,'s'), 
)
    
model(
    toleranceMoveToCore = 0.01,
    toleranceInterruptSimulation = 0.1,
    maximumEdgeSpecies = 3e5,
    filterReactions = True,
    filterThreshold = 5e8,
    minCoreSizeForPrune = 50,
    minSpeciesExistIterationsForPrune = 4,
)

pressureDependence(
    method='modified strong collision',
    #method = 'reservoir state',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,2500,'K',8),
    pressures=(0.01,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
    maximumAtoms=16,
    completedNetworks=['CH2O2', 'CO3'],
    
)


simulator(
    atol = 1e-16,
    rtol = 1e-08,
    sens_atol = 1e-06,
    sens_rtol = 0.0001,
)

generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=8, 
    maximumRadicalElectrons=1, 
    maximumSingletCarbenes=1,
    maximumCarbeneRadicals=0, 
    allowSingletO2 = True,
)

options(
    units = "si",
    generateSeedEachIteration = True,
    generateOutputHTML = True,
    generatePlots = True,
    saveSimulationProfiles = True,
    saveEdgeSpecies = True,
    keepIrreversible = True,
    verboseComments = False,
)



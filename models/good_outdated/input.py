
thermolibs = [
'C1_C2_Fluorine', #adding Siddha's as second most trusted 
'PFCA_thermo',
#'NCSU_C2_C8_PFAS', #adding Westmoreland's thermo as the first most trusted
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
reactionLibraries = ['FFCM1(-)', 'halogens_pdep', 'ANL_Brown_pdep', "PFAS_charged_reactions", "PFAS_HPL_multiple_matches"],
seedMechanisms = [],
kineticsDepositories = ['training'],
kineticsFamilies = ['default','PFAS', 'halogens','Disproportionation-Y'],
frequenciesLibraries = ['halogens_G4'],
kineticsEstimator = 'rate rules',
)


species(
    label = 'CH3F',
    reactive = True,
    structure = SMILES('CF')
)


species(
    label = 'O2',
    reactive = True,
    structure = SMILES('[O][O]')
)


species(
    label = 'N2',
    reactive = False,
    structure = SMILES('N#N')
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
    maximumOxygenAtoms=6,
    #maximumHeavyAtoms=24,
    maximumRadicalElectrons=2,
    maximumSingletCarbenes=1,
    maximumCarbeneRadicals=0,
    allowSingletO2 = True,
)

options(
    units = "si",
    generateSeedEachIteration = True,
    generateOutputHTML = False,
    generatePlots = False,
    saveSimulationProfiles = False,
    saveEdgeSpecies = True,
    keepIrreversible = True,
    verboseComments = False,
)
    
#stoichiometric combustion : 0.6666666666666666CH3F + 1.0(O2 + 3.76N2) => 0.6666666666666666CO2 + 0.6666666666666666HF + 0.0CF2O + 0.6666666666666666H2O + 3.76N2
#stoichiometric initial composition: {'CH3F': 0.6666666666666666, 'O2': 1.0, 'N2': 3.76}
    
simpleReactor(
        temperature=[(1000,'K'),(2000,'K')],
        pressure= [(1.0,'bar'),(10.0,'bar')],
        nSims=10,
        initialMoleFractions={
        "CH3F": [0.4666666666666666, 0.8666666666666667], #going a bit over and a bit under stoic conditions
        "O2": 1.0,
        "N2": 3.76,
        },
        # terminationConversion={
        # 'halogen': 0.999,
        # },
        #terminationRateRatio=1e-4,
        #terminationTime=(10,'s'),
        terminationTime=(3,'s'),
        #sensitivity=['halogen','OH'],
        #sensitivityThreshold=0.001,
        )
        
model(
    toleranceMoveToCore = 0.1,
    toleranceInterruptSimulation = 0.1,
    maximumEdgeSpecies = 3e5,
    filterReactions = True,
    filterThreshold = 5e8,
    minCoreSizeForPrune = 50,
    minSpeciesExistIterationsForPrune = 4,
)

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,2500,'K',8),
    pressures=(0.01,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
    maximumAtoms=16,
)



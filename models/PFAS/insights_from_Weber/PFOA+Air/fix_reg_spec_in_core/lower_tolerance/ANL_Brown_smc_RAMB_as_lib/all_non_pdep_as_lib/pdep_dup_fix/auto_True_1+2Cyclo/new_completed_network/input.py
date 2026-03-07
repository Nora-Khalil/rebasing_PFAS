thermolibs = [
'NCSU_2025_revised',
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
reactionLibraries = ['FFCM1(-)', 'halogens_pdep', "2-BTP/full", "PFAS_HPL_Test/R_Addition_MultipleBond", "PFAS_charged_reactions", "PFAS_HPL_multiple_matches"],
seedMechanisms = [ "ANL_Brown_pdep", 'Brown_difluoromethane_seed'],
kineticsDepositories = ['training'],
kineticsFamilies = ['default','PFAS', 'halogens','Disproportionation-Y'],
frequenciesLibraries = ['halogens_G4'],
kineticsEstimator = 'rate rules',
)


#Only put PFPA in core initially, want to see what RMG will place in the core on its own
species(
    label = 'PFOA',
    reactive = True,
    structure = SMILES('C(=O)(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)O')
)

#in N2 bath gas 
species(
    label = 'N2',
    reactive = False,
    structure = SMILES('N#N')
)


#air as O2 
species(
    label = 'O2',
    reactive = True,
    structure = SMILES('[O][O]')
)

#water in gaseous state is an impurity in the reactor
species(
    label = 'H2O',
    reactive = True,
    structure = SMILES('O')
)

#making sure all the species that are monitored by Weber are in the core

species(
    label = "HF",
    reactive = True,
    structure = SMILES('F[H]')
)


species(
    label = "C2F4",
    reactive = True,
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {6,S}
5 C u0 p0 c0 {1,S} {2,S} {6,D}
6 C u0 p0 c0 {3,S} {4,S} {5,D}""")
)


species(
    label = "CO2",
    reactive = True,
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}""")
)


species(
    label = "C2F6",
    reactive = True,
    structure = SMILES('C(C(F)(F)(F))(F)(F)(F)')
)


species(
    label = "CO",
    reactive = True,
    structure = adjacencyList("""1 O u0 p1 c+1 {2,T}
2 C u0 p1 c-1 {1,T}""")
)


species(
    label = "CF4",
    reactive = True,
    structure = SMILES('FC(F)(F)F')
)


species(
    label = "COF2",
    reactive = True,
    structure = SMILES('O=C(F)F')
)



species(
    label = "C7F14singletcarbene",
    reactive = True,
    structure = adjacencyList("""1  F u0 p3 c0 {18,S}
2  F u0 p3 c0 {18,S}
3  F u0 p3 c0 {16,S}
4  F u0 p3 c0 {16,S}
5  F u0 p3 c0 {15,S}
6  F u0 p3 c0 {15,S}
7  F u0 p3 c0 {17,S}
8  F u0 p3 c0 {17,S}
9  F u0 p3 c0 {19,S}
10 F u0 p3 c0 {19,S}
11 F u0 p3 c0 {20,S}
12 F u0 p3 c0 {20,S}
13 F u0 p3 c0 {20,S}
14 F u0 p3 c0 {21,S}
15 C u0 p0 c0 {5,S} {6,S} {16,S} {17,S}
16 C u0 p0 c0 {3,S} {4,S} {15,S} {18,S}
17 C u0 p0 c0 {7,S} {8,S} {15,S} {19,S}
18 C u0 p0 c0 {1,S} {2,S} {16,S} {20,S}
19 C u0 p0 c0 {9,S} {10,S} {17,S} {21,S}
20 C u0 p0 c0 {11,S} {12,S} {13,S} {18,S}
21 C u0 p1 c0 {14,S} {19,S}""")
    
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
        #terminationRateRatio=1e-6,
        terminationTime=(4,'s'), #source says its in PFR for 2 and 25 s, should i try and change this to 2 
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

#left the same as when David ran his refrigerants
generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=8, 
    # maximumOxygenAtoms=4,
    # maximumRadicalElectrons=2,
    # maximumSingletCarbenes=1,
    # maximumCarbeneRadicals=0,
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



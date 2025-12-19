
# File containing the model parameters and which acquire the experimental data for validation.
# This file feeds the RunmeDDE.py file
#
# Created by Luca Rossini on 11 March 2025
# Last update 28 August 2025
#
# E-mail: luca.rossini@ulb.be


# Import list

import pandas as pd
import numpy as np


# Absorb the temperature array needed by the model

InputFile = pd.read_excel('TemperatureInput.xlsx')
DailyTemp = InputFile['Average_Temperature']


# Absorb the experimental data

ExpDataFile_OldCode = pd.read_excel('DifferentialRep-AbsTime.xlsx', sheet_name='Temp_28')
ExpDataDay_OldCode = ExpDataFile_OldCode['Time']
ExpDataEgg_OldCode = ExpDataFile_OldCode['Egg']
ExpDataL1_OldCode = ExpDataFile_OldCode['L1']
ExpDataL2_OldCode = ExpDataFile_OldCode['L2']
ExpDataL3_OldCode = ExpDataFile_OldCode['L3']
ExpDataPupa_OldCode = ExpDataFile_OldCode['Pupa']
ExpDataAdMale_OldCode = ExpDataFile_OldCode['AdMale']
ExpDataAdFemale_OldCode = ExpDataFile_OldCode['AdFemale']


# List of the parameters from life tables analyses

    # Sex ratio

SR = 0.5

    # Fertility rate

alpha = 37000.000000
gamma = 17.761392
Lambda = 26.036039
tau = 22.944989
delta = 4.0000

FertPar = [alpha, gamma, Lambda, tau, delta]

    # Development Egg - Sharpe and De Michele

A_Egg = 4.782932
B_Egg = -113.997964
C_Egg = 20.442876
D_Egg = 272.1352
E_Egg = 6.756630
F_Egg = -137.529823

DevPar_Egg = [A_Egg, B_Egg, C_Egg, D_Egg, E_Egg, F_Egg]

    # Development L1 - Brière
    
a_L1 = 0.001134
TL_L1 = 3.148968
TM_L1 = 32.023358
m_L1 = 5.686684

DevPar_L1 = [a_L1, TL_L1, TM_L1, m_L1]

    # Development L2 - Logan
    
psi_L2 = 0.043368
rho_L2 = 0.116276
TM_L2 = 32.329243
DeltaT_L2 = 2.217094

DevPar_L2 = [psi_L2, rho_L2, TM_L2, DeltaT_L2]

    # Development L3 - Lactin-1

rho_L3 = 0.185309
TM_L3 = 32.973793
DeltaT_L3 = 5.374509

DevPar_L3 = [rho_L3, TM_L3, DeltaT_L3]
    
    # Development P - Brière

a_P = 0.000046
TL_P = 6.602687
TM_P = 35.271650
m_P = 1.00

DevPar_P = [a_P, TL_P, TM_P, m_P]
    
    # Survival Ad Males - Lactin-1

a_AM = 0.106797
TM_AM = 34.999382
DeltaT_AM = 9.333894
    
DevPar_AM = [a_AM, TM_AM, DeltaT_AM]
    
    # Survival Ad Females - Lactin-1

a_AF = 0.097196
TM_AF = 35.000000
DeltaT_AF = 10.243252

DevPar_AF = [a_AF, TM_AF, DeltaT_AF]

    # Mortality Egg - Exponential

a_Mort_Egg = 0.012533
b_Mort_Egg = -0.482365
c_Mort_Egg = 2.256795

MortPar_Egg = [a_Mort_Egg, b_Mort_Egg, c_Mort_Egg]
    
    # Mortality L1 - Exponential

a_Mort_L1 = 0.0143505
b_Mort_L1 = -0.5538988
c_Mort_L1 = 2.584712

MortPar_L1 = [a_Mort_L1, b_Mort_L1, c_Mort_L1]
    
    # Mortality L2 - Exponential
    
a_Mort_L2 = 0.017029
b_Mort_L2 = -0.658377
c_Mort_L2 = 3.308932

MortPar_L2 = [a_Mort_L2, b_Mort_L2, c_Mort_L2]
    
    # Mortality L3 - Exponential

a_Mort_L3 = 0.018840
b_Mort_L3 = -0.727484
c_Mort_L3 = 3.652471

MortPar_L3 = [a_Mort_L3, b_Mort_L3, c_Mort_L3]

    # Mortality Pupa - Constant value

MortPar_P = 0.1

# Lag functions parameters - for DDE model

         # Egg lag - 2nd order rational (ax^2+bx+c)/(x+d)

a_LagEgg = 0.4114855
b_LagEgg = -18.8779629
c_LagEgg = 241.8451923
d_LagEgg = 12.6277356

LagPar_Egg = [a_LagEgg, b_LagEgg, c_LagEgg, d_LagEgg]

        # L1 lag - 2nd order polynomial ax^2+bx+c
        
a_LagL1 = 0.0113729
b_LagL1 = -0.6306989
c_LagL1 = 9.7507676

LagPar_L1 = [a_LagL1, b_LagL1, c_LagL1]

        # L2 lag - 2nd order polynomial ax^2+bx+c
        
a_LagL2 = 0.0106494
b_LagL2 = -0.6290057
c_LagL2 = 10.422407

LagPar_L2 = [a_LagL2, b_LagL2, c_LagL2]

        # L3 lag - 2nd order polynomial ax^2+bx+c

a_LagL3 = 0.015178
b_LagL3 = -0.8782897
c_LagL3 = 14.0351341

LagPar_L3 = [a_LagL3, b_LagL3, c_LagL3]

        # Pupa Lag - Exponential a + b exp(cx)
        
a_LagP = 4.5004072
b_LagP = 2126.6427
c_LagP = -0.4302271

LagPar_P = [a_LagP, b_LagP, c_LagP]

        # AF2 Lag - PreOviposition - 2nd order polynomial ax^2+bx+c
        
a_LagAF2 = 0.0070255
b_LagAF2 = -0.2766729
c_LagAF2 = 3.7522346

LagPar_AF2 = [a_LagAF2, b_LagAF2, c_LagAF2]

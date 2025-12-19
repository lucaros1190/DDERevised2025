from matplotlib.pyplot import pause

# import list

from Parameters import *
from Functions import *
import numpy as np
import scipy.integrate as sc 
from ddeint import ddeint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import time


# Timer start

t_start = time.time()


# Acquisition of the time span: consider that the Delta t subdivision affects the precision of the algorithm. Make some attempts to find a good compromise. The Delta t is the variable IntegrationStep, to define according to the needs of the user.

IntegrationStep = 7500

t_span = np.linspace(0, len(DailyTemp), IntegrationStep)


# Define the initial history for the DDE: change the array inside the function to modify the initial condition

def history(t):
    s = 2 * np.heaviside(t, 0.5) * np.array([50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    return s


# Compute the normalisation coefficients for the development rates and adult survival

[NormCoeff_Egg, NormCoeff_L1, NormCoeff_L2, NormCoeff_L3, NormCoeff_P, NormCoeff_AF, NormCoeff_AM] = NormCoefficients(DevPar_Egg, DevPar_L1, DevPar_L2, DevPar_L3, DevPar_P, DevPar_AF, DevPar_AM)


# Solve the equation using 'ddeint'

print('\nPlease wait... Solving DDE model \n')

SolDDE_NormRevised = ddeint(DDE_FuncNormRevised, history, t_span, fargs = (FertPar, MortPar_Egg, MortPar_L1, MortPar_L2, MortPar_L3, MortPar_P, DevPar_Egg, DevPar_L1, DevPar_L2, DevPar_L3, DevPar_P, DevPar_AF, DevPar_AM, LagPar_Egg, LagPar_L1, LagPar_L2, LagPar_L3, LagPar_P, LagPar_AF2, SR, DailyTemp, NormCoeff_Egg, NormCoeff_L1, NormCoeff_L2, NormCoeff_L3, NormCoeff_P, NormCoeff_AF, NormCoeff_AM))


# Timer stop

t_f = time.time() - t_start

print('\nExecution time:', t_f, 's\n')

# Plot the results


    # Mosaic plot - First part
    
fig = plt.figure(8)

pltArray = fig.subplots(2, 2) # (rows, columns)

pltArray[0, 0].plot(t_span, SolDDE_NormRevised[:,0] + SolDDE_NormRevised[:,1], color = 'blue', label = 'Estimated')
pltArray[0, 0].scatter(ExpDataDay_OldCode, ExpDataEgg_OldCode, marker ='o', color='red', zorder=4, s=30, label='Experimental data')
pltArray[0, 0].set_title('Egg', fontsize = 17)
pltArray[0, 0].set_ylabel('N. of individuals', fontsize = 13)
pltArray[0, 0].legend(['Estimated', 'Experimental data'], loc='upper right', fontsize = 13)


pltArray[0, 1].plot(t_span, SolDDE_NormRevised[:,2] + SolDDE_NormRevised[:,3], color = 'blue', label = 'Estimated')
pltArray[0, 1].scatter(ExpDataDay_OldCode, ExpDataL1_OldCode, marker ='o', color='red', zorder=4, s=30, label='Experimental data')
pltArray[0, 1].set_title('L1', fontsize = 17)


pltArray[1, 0].plot(t_span, SolDDE_NormRevised[:,4] + SolDDE_NormRevised[:,5], color = 'blue', label = 'Estimated')
pltArray[1, 0].scatter(ExpDataDay_OldCode, ExpDataL2_OldCode, marker ='o', color='red', zorder=4, s=30, label='Experimental data')
pltArray[1, 0].set_ylabel('N. of individuals', fontsize = 13)
pltArray[1, 0].set_xlabel('Time (days)', fontsize = 13)
pltArray[1, 0].set_title('L2', fontsize = 17)


pltArray[1, 1].plot(t_span, SolDDE_NormRevised[:,6] + SolDDE_NormRevised[:,7], color = 'blue', label = 'Estimated')
pltArray[1, 1].scatter(ExpDataDay_OldCode, ExpDataL3_OldCode, marker ='o', color='red', zorder=4, s=30, label='Experimental data')
pltArray[1, 1].set_xlabel('Time (days)', fontsize = 13)
pltArray[1, 1].set_title('L3', fontsize = 17)


    # Mosaic plot - Second part

figBis = plt.figure(9)

pltArrayBis = figBis.subplots(2, 2)

pltArrayBis[0, 0].plot(t_span, SolDDE_NormRevised[:,8] + SolDDE_NormRevised[:,9], color = 'blue', label = 'Estimated')
pltArrayBis[0, 0].scatter(ExpDataDay_OldCode, ExpDataPupa_OldCode, marker ='o', color='red', zorder=4, s=30, label='Experimental data')
pltArrayBis[0, 0].set_ylabel('N. of individuals', fontsize = 13)
pltArrayBis[0, 0].set_title('Pupa', fontsize = 17)


pltArrayBis[0, 1].plot(t_span, SolDDE_NormRevised[:,10], color = 'blue', label = 'Estimated')
pltArrayBis[0, 1].scatter(ExpDataDay_OldCode, ExpDataAdFemale_OldCode, marker ='o', color='red', zorder=4, s=30, label='Experimental data')
pltArrayBis[0, 1].set_title('Adult males', fontsize = 17)
pltArrayBis[0, 1].set_xlabel('Time (days)', fontsize = 13)


pltArrayBis[1, 0].plot(t_span, SolDDE_NormRevised[:,11] + SolDDE_NormRevised[:,12], color = 'blue', label = 'Estimated')
pltArrayBis[1, 0].scatter(ExpDataDay_OldCode, ExpDataAdMale_OldCode, marker ='o', color='red', zorder=4, s=30, label='Experimental data')
pltArrayBis[1, 0].set_xlabel('Time (days)', fontsize = 13)
pltArrayBis[1, 0].set_ylabel('N. of individuals', fontsize = 13)
pltArrayBis[1, 0].set_title('Adult females', fontsize = 17)

pltArrayBis[1, 1].plot(1,1)
pltArrayBis[-1, -1].axis('off')



plt.show()

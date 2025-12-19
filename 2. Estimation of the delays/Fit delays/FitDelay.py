
# Python script to fit all the minimum delay functions in once.
# Input file: FitData.csv
#
# Created by Luca Rossini on 27 August 2025
# Last update 27 August 2025
#
# E-mail: luca.rossini@ulb.be


# Import list

import pandas as pd
from math import *
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats.distributions import chi2
import numpy as np
import matplotlib.pyplot as plt


# Read the data to fit and plot
  
Data = pd.read_csv('FitData.csv', sep='\t', header = 0)

# Set the header of the dataset

Data.columns = ["x", "y", "err_y"]

x = Data['x']
y = Data['y']

err_y = Data['err_y']


# Definition of the functions to fit

    # 1. Rational 2nd order
    
def ratpoly_2ord(x, a, b, c, d):
    return (a * x**2 + b * x + c) / (x + d)

    # 2. Exponential 1

def expofun_1(x, a, b):
    return a * np.exp(b * x)
    
    # 3. Second-order polynomial

def poly_2(x, a, b, c):
    return a * x**2 + b * x + c
    
    # 4. Rational 1st order
    
def ratpoly_1ord(x, a, b, c):
    return (a * x + b) / (x + c)
      
    # 5. Exponential 2

def expofun_2(x, a, b, c):
    return a * np.pow(x,b) + c

    # 6. Exponential 3

def expofun_3(x, a, b, c):
    return a + b * np.exp(c * x)

    
    
# Fitting the function 1
    
popt_1, pcov_1 = curve_fit(ratpoly_2ord, x, y, bounds = ([-100000, -100000, -100000, -100000], [50000, 50000, 50000, 50000]), p0 = (1, 1, 1, 1), method='dogbox', max_nfev=10000)

    # Best fit values
    
a_1 = popt_1[0]
b_1 = popt_1[1]
c_1 = popt_1[2]
d_1 = popt_1[3]

    # Perameters error

    # np.diag extracts the diagonal elements of the covariance matrix
    # np.sqrt computes the square root of each element in an array

perr_1 = np.sqrt(np.diag(pcov_1))

err_a_1 = perr_1[0]
err_b_1 = perr_1[1]
err_c_1 = perr_1[2]
err_d_1 = perr_1[3]

    # Goodness-of-fit indicators

x_fun = np.linspace(0, 40, 5000)
fitted_fun_1 = ratpoly_2ord(x, *popt_1)
fitted_plot_1 = ratpoly_2ord(x_fun, *popt_1)

    # Calculating R-squared

resid = y - fitted_fun_1
ss_res = np.sum(resid**2)
ss_tot = np.sum((y - np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)

    # Number of degrees of freedom (NDF)

ndf = len(x) - 4

    # Calculate the chi-squared (with error below the fraction)

chi_sq = 0
for i in range(len(x)):
    chi_sq = pow((y[i] - fitted_fun_1[i])/ fitted_fun_1[i], 2)

    # Calculate the P-value from chi-square

Pvalue = 1 - chi2.sf(chi_sq, ndf)

    # Calculate AIC and BIC

AIC = 2 * 4 - 2 * np.log(ss_res/len(x))
BIC = 4 * np.log(len(x)) - 2 * np.log(ss_res/len(x))

    # Print the results

print('\n Rational second order fit (a * x**2 + b * x + c) / (x + d) results: \n')

    # Define the parameters' name
    
parname = ('a = ', 'b = ', 'c = ', 'd = ')

for i in range(len(popt_1)):
    print(parname[i] + str(round(popt_1[i], 7)) + ' +/- ' + str(round(perr_1[i],7)))

print(' R-squared = ', round(r_squared, 5))
print(' Chi-squared = ', round(chi_sq, 5))
print(' P-value = ', round(Pvalue, 7))
print(' Number of degrees of freedom (NDF) =', ndf)
print(' Akaike Information Criterion (AIC):', round(AIC, 5))
print(' Bayesian Information Criterion (BIC)', round(BIC, 5))
print('\n')

print(' Covariance matrix: \n')
    # Define the row names to print
rowname = (' ', 'a', 'b', 'c', 'd')
print(rowname[0] + '     \t' + rowname[1] + '     \t' + rowname[2] + '     \t' + rowname[3] + rowname[4] + '     \t')

for i in range(len(pcov_1)):
    print(rowname[i+1] + ' ' + str(pcov_1[i]))

print(' ')



# Fitting the function 2
    
popt_2, pcov_2 = curve_fit(expofun_1, x, y, bounds = ([-100000, -100000], [50000, 50000]), p0 = (1, -1), method='dogbox', max_nfev=10000)

    # Best fit values
    
a_2 = popt_2[0]
b_2 = popt_2[1]

    # Perameters error

perr_2 = np.sqrt(np.diag(pcov_2))

err_a_2 = perr_2[0]
err_b_2 = perr_2[1]

    # Goodness-of-fit indicators

x_fun = np.linspace(0, 40, 5000)
fitted_fun_2 = expofun_1(x, *popt_2)
fitted_plot_2 = expofun_1(x_fun, *popt_2)

    # Calculating R-squared

resid = y - fitted_fun_2
ss_res = np.sum(resid**2)
ss_tot = np.sum((y - np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)

    # Number of degrees of freedom (NDF)

ndf = len(x) - 2

    # Calculate the chi-squared (with error below the fraction)

chi_sq = 0
for i in range(len(x)):
    chi_sq = pow((y[i] - fitted_fun_2[i])/ fitted_fun_2[i], 2)

    # Calculate the P-value from chi-square

Pvalue = 1 - chi2.sf(chi_sq, ndf)

    # Calculate AIC and BIC

AIC = 2 * 4 - 2 * np.log(ss_res/len(x))
BIC = 4 * np.log(len(x)) - 2 * np.log(ss_res/len(x))

    # Print the results

print('\n Exponential-1 fit a * exp(b * x) results: \n')

    # Define the parameters' name
    
parname = ('a = ', 'b = ')

for i in range(len(popt_2)):
    print(parname[i] + str(round(popt_2[i], 7)) + ' +/- ' + str(round(perr_2[i],7)))

print(' R-squared = ', round(r_squared, 5))
print(' Chi-squared = ', round(chi_sq, 5))
print(' P-value = ', round(Pvalue, 7))
print(' Number of degrees of freedom (NDF) =', ndf)
print(' Akaike Information Criterion (AIC):', round(AIC, 5))
print(' Bayesian Information Criterion (BIC)', round(BIC, 5))
print('\n')

print(' Covariance matrix: \n')
    # Define the row names to print
rowname = (' ', 'a', 'b')
print(rowname[0] + '     \t' + rowname[1] + '     \t' + rowname[2] + '     \t' + '     \t')

for i in range(len(pcov_2)):
    print(rowname[i+1] + ' ' + str(pcov_2[i]))

print(' ')



# Fitting the function 3
    
popt_3, pcov_3 = curve_fit(poly_2, x, y, bounds = ([-100000, -100000, -100000], [50000, 50000, 50000]), p0 = (1, 1, 1), method='dogbox', max_nfev=10000)

    # Best fit values
    
a_3 = popt_3[0]
b_3 = popt_3[1]
c_3 = popt_3[2]

    # Perameters error

perr_3 = np.sqrt(np.diag(pcov_3))

err_a_3 = perr_3[0]
err_b_3 = perr_3[1]
err_c_3 = perr_3[2]

    # Goodness-of-fit indicators

x_fun = np.linspace(0, 40, 5000)
fitted_fun_3 = poly_2(x, *popt_3)
fitted_plot_3 = poly_2(x_fun, *popt_3)

    # Calculating R-squared

resid = y - fitted_fun_3
ss_res = np.sum(resid**2)
ss_tot = np.sum((y - np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)

    # Number of degrees of freedom (NDF)

ndf = len(x) - 3

    # Calculate the chi-squared (with error below the fraction)

chi_sq = 0
for i in range(len(x)):
    chi_sq = pow((y[i] - fitted_fun_3[i])/ fitted_fun_3[i], 2)

    # Calculate the P-value from chi-square

Pvalue = 1 - chi2.sf(chi_sq, ndf)

    # Calculate AIC and BIC

AIC = 2 * 4 - 2 * np.log(ss_res/len(x))
BIC = 4 * np.log(len(x)) - 2 * np.log(ss_res/len(x))

    # Print the results

print('\n 2nd order polynomial fit a * x^2 + b * x + c results: \n')

    # Define the parameters' name
    
parname = ('a = ', 'b = ', 'c = ')

for i in range(len(popt_3)):
    print(parname[i] + str(round(popt_3[i], 7)) + ' +/- ' + str(round(perr_3[i],7)))

print(' R-squared = ', round(r_squared, 5))
print(' Chi-squared = ', round(chi_sq, 5))
print(' P-value = ', round(Pvalue, 7))
print(' Number of degrees of freedom (NDF) =', ndf)
print(' Akaike Information Criterion (AIC):', round(AIC, 5))
print(' Bayesian Information Criterion (BIC)', round(BIC, 5))
print('\n')

print(' Covariance matrix: \n')
    # Define the row names to print
rowname = (' ', 'a', 'b', 'c')
print(rowname[0] + '     \t' + rowname[1] + '     \t' + rowname[2] + '     \t' + rowname[3] + '     \t' + '     \t')

for i in range(len(pcov_3)):
    print(rowname[i+1] + ' ' + str(pcov_3[i]))

print(' ')



# Fitting the function 4
    
popt_4, pcov_4 = curve_fit(ratpoly_1ord, x, y, bounds = ([-100000, -100000, -1], [50000, 50000, 50000]), p0 = (1, 1, 1), method='dogbox', max_nfev=10000)

    # Best fit values
    
a_4 = popt_4[0]
b_4 = popt_4[1]
c_4 = popt_4[2]

    # Perameters error

perr_4 = np.sqrt(np.diag(pcov_4))

err_a_4 = perr_4[0]
err_b_4 = perr_4[1]
err_c_4 = perr_4[2]

    # Goodness-of-fit indicators

x_fun = np.linspace(0, 40, 5000)
fitted_fun_4 = ratpoly_1ord(x, *popt_4)
fitted_plot_4 = ratpoly_1ord(x_fun, *popt_4)

    # Calculating R-squared

resid = y - fitted_fun_4
ss_res = np.sum(resid**2)
ss_tot = np.sum((y - np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)

    # Number of degrees of freedom (NDF)

ndf = len(x) - 3

    # Calculate the chi-squared (with error below the fraction)

chi_sq = 0
for i in range(len(x)):
    chi_sq = pow((y[i] - fitted_fun_4[i])/ fitted_fun_4[i], 2)

    # Calculate the P-value from chi-square

Pvalue = 1 - chi2.sf(chi_sq, ndf)

    # Calculate AIC and BIC

AIC = 2 * 4 - 2 * np.log(ss_res/len(x))
BIC = 4 * np.log(len(x)) - 2 * np.log(ss_res/len(x))

    # Print the results

print('\n Rational first order fit (a * x + b) / (x + c) results: \n')

    # Define the parameters' name
    
parname = ('a = ', 'b = ', 'c = ')

for i in range(len(popt_4)):
    print(parname[i] + str(round(popt_4[i], 7)) + ' +/- ' + str(round(perr_4[i],7)))

print(' R-squared = ', round(r_squared, 5))
print(' Chi-squared = ', round(chi_sq, 5))
print(' P-value = ', round(Pvalue, 7))
print(' Number of degrees of freedom (NDF) =', ndf)
print(' Akaike Information Criterion (AIC):', round(AIC, 5))
print(' Bayesian Information Criterion (BIC)', round(BIC, 5))
print('\n')

print(' Covariance matrix: \n')
    # Define the row names to print
rowname = (' ', 'a', 'b', 'c')
print(rowname[0] + '     \t' + rowname[1] + '     \t' + rowname[2] + '     \t' + rowname[3] + '     \t' + '     \t')

for i in range(len(pcov_4)):
    print(rowname[i+1] + ' ' + str(pcov_4[i]))

print(' ')



# Fitting the function 5
    
popt_5, pcov_5 = curve_fit(expofun_2, x, y, bounds = ([-100000, -100000, -100000], [50000, 50000, 50000]), p0 = (1, 1, 1), method='dogbox', max_nfev=10000)

    # Best fit values
    
a_5 = popt_5[0]
b_5 = popt_5[1]
c_5 = popt_5[2]

    # Perameters error

perr_5 = np.sqrt(np.diag(pcov_5))

err_a_5 = perr_5[0]
err_b_5 = perr_5[1]
err_c_5 = perr_5[2]

    # Goodness-of-fit indicators

x_fun = np.linspace(0, 40, 5000)
fitted_fun_5 = expofun_2(x, *popt_5)
fitted_plot_5 = expofun_2(x_fun, *popt_5)

    # Calculating R-squared

resid = y - fitted_fun_5
ss_res = np.sum(resid**2)
ss_tot = np.sum((y - np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)

    # Number of degrees of freedom (NDF)

ndf = len(x) - 3

    # Calculate the chi-squared (with error below the fraction)

chi_sq = 0
for i in range(len(x)):
    chi_sq = pow((y[i] - fitted_fun_5[i])/ fitted_fun_5[i], 2)

    # Calculate the P-value from chi-square

Pvalue = 1 - chi2.sf(chi_sq, ndf)

    # Calculate AIC and BIC

AIC = 2 * 4 - 2 * np.log(ss_res/len(x))
BIC = 4 * np.log(len(x)) - 2 * np.log(ss_res/len(x))

    # Print the results

print('\n Exponential-2 fit a * x^b + c results: \n')

    # Define the parameters' name
    
parname = ('a = ', 'b = ', 'c = ')

for i in range(len(popt_5)):
    print(parname[i] + str(round(popt_5[i], 7)) + ' +/- ' + str(round(perr_5[i],7)))

print(' R-squared = ', round(r_squared, 5))
print(' Chi-squared = ', round(chi_sq, 5))
print(' P-value = ', round(Pvalue, 7))
print(' Number of degrees of freedom (NDF) =', ndf)
print(' Akaike Information Criterion (AIC):', round(AIC, 5))
print(' Bayesian Information Criterion (BIC)', round(BIC, 5))
print('\n')

print(' Covariance matrix: \n')
    # Define the row names to print
rowname = (' ', 'a', 'b', 'c')
print(rowname[0] + '     \t' + rowname[1] + '     \t' + rowname[2] + '     \t' + rowname[3] + '     \t' + '     \t')

for i in range(len(pcov_5)):
    print(rowname[i+1] + ' ' + str(pcov_5[i]))

print(' ')



# Fitting the function 6
    
popt_6, pcov_6 = curve_fit(expofun_3, x, y, bounds = ([0, 0, -100000], [50000, 50000, 0]), p0 = (1, 1, -1), method='dogbox', max_nfev=10000)

    # Best fit values
    
a_6 = popt_6[0]
b_6 = popt_6[1]
c_6 = popt_6[2]

    # Perameters error

perr_6 = np.sqrt(np.diag(pcov_6))

err_a_6 = perr_6[0]
err_b_6 = perr_6[1]
err_c_6 = perr_6[2]

    # Goodness-of-fit indicators

x_fun = np.linspace(0, 40, 5000)
fitted_fun_6 = expofun_3(x, *popt_6)
fitted_plot_6 = expofun_3(x_fun, *popt_6)

    # Calculating R-squared

resid = y - fitted_fun_6
ss_res = np.sum(resid**2)
ss_tot = np.sum((y - np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)

    # Number of degrees of freedom (NDF)

ndf = len(x) - 3

    # Calculate the chi-squared (with error below the fraction)

chi_sq = 0
for i in range(len(x)):
    chi_sq = pow((y[i] - fitted_fun_6[i])/ fitted_fun_6[i], 2)

    # Calculate the P-value from chi-square

Pvalue = 1 - chi2.sf(chi_sq, ndf)

    # Calculate AIC and BIC

AIC = 2 * 4 - 2 * np.log(ss_res/len(x))
BIC = 4 * np.log(len(x)) - 2 * np.log(ss_res/len(x))

    # Print the results

print('\n Exponential-3 fit a + exp(b * c) results: \n')

    # Define the parameters' name
    
parname = ('a = ', 'b = ', 'c = ')

for i in range(len(popt_6)):
    print(parname[i] + str(round(popt_6[i], 7)) + ' +/- ' + str(round(perr_6[i],7)))

print(' R-squared = ', round(r_squared, 5))
print(' Chi-squared = ', round(chi_sq, 5))
print(' P-value = ', round(Pvalue, 7))
print(' Number of degrees of freedom (NDF) =', ndf)
print(' Akaike Information Criterion (AIC):', round(AIC, 5))
print(' Bayesian Information Criterion (BIC)', round(BIC, 5))
print('\n')

print(' Covariance matrix: \n')
    # Define the row names to print
rowname = (' ', 'a', 'b', 'c')
print(rowname[0] + '     \t' + rowname[1] + '     \t' + rowname[2] + '     \t' + rowname[3] + '     \t' + '     \t')

for i in range(len(pcov_6)):
    print(rowname[i+1] + ' ' + str(pcov_6[i]))

print(' ')




# Plot the data

# Define the limits for x and y axes

x_min, x_max = 0, 40  # Example limits for x-axis
y_min, y_max = 0, 20  # Example limits for y-axis

plt.figure(1)

plt.scatter(x, y, color="C0", alpha=0.5, label=f'Experimental data')
plt.errorbar(x, y, yerr=err_y, fmt='o') #fmt specifies point as circles
plt.plot(x_fun, fitted_plot_1, color="C1", alpha=0.7, label=f'y = A(x)^2 / B(x)')
plt.plot(x_fun, fitted_plot_2, color="C2", alpha=0.7, label=f'y = a exp(b x)')
plt.plot(x_fun, fitted_plot_3, color="C3", alpha=0.7, label=f'y = ax^2 + bx + c)')
plt.plot(x_fun, fitted_plot_4, color="C4", alpha=0.7, label=f'y = (a * x + b) / (x + c)')
plt.plot(x_fun, fitted_plot_5, color="C5", alpha=0.7, label=f'y = a * x^b + c')
plt.plot(x_fun, fitted_plot_6, color="C6", alpha=0.7, label=f'y = a + b exp(c x)')

plt.xlabel('Temperature (°C)')
plt.ylabel('Delay (days)')
plt.legend()
plt.title('Preoviposition time')

# Set the limits for x and y axes

plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)

plt.show()

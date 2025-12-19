
# Python script that resamples the values from life tables data. Please, USE THE STANDARD EXCEL TEMPLATE, this script is built ad hoc for that.
#
# This script is helpful to compute the minimum development time of each life stage from life tables datasets written as differential representation.
#
# Created by Luca Rossini on 26 January 2025
# Last update: 27 January 2025
#
# E-mail: luca.rossini@ulb.be


# Import list

import random
import statistics
import multiprocessing
import pandas as pd
import time
import numpy as np


# Timer start

t_start = time.time()


# Step 0: Declarations

    # Number of samples per subgroup

    # Development
    
N = 3

    # Pre-oviposition
    
N_PreOvi = 3

    # Number of bootstrap cycles
    
BootstrapIterations = 10000000

    # Number of CPU cores
    
NumCores = multiprocessing.cpu_count()
multiprocessing.set_start_method('fork', force = True)


# Step 1: Load data from Excel

DevData = pd.read_excel('LifeTablesDataset.xlsx', sheet_name='Individual-LifeHistory')

# Define the temperature mappings - You can change the label names on the left of each line, but don't change the range, as it is ad hoc to select the columns properly!

RearingConditions = {'6': list(range(2, 10)),
                     '9': list(range(15, 23)),
                     '13': list(range(28, 36)),
                     '18': list(range(41, 49)),
                     '20': list(range(54, 62)),
                     '24': list(range(67, 75)),
                     '25': list(range(80, 88)),
                     '26': list(range(93, 101)),
                     '27': list(range(106, 114)),
                     '28': list(range(119, 127)),
                     '29': list(range(132, 140)),
                     '31': list(range(145, 153)),
                     '32': list(range(158, 166)),
                     '33': list(range(171, 179))}
                       

# Define the stage names assigned to the columns - You can change the name of the stages at your best convenience and based on the definitions of your Excel file

stages = ['Egg',
          'L1',
          'L2',
          'L3',
          'NotAvailable-1',
          'NotAvailable-2',
          'NotAvailable-3',
          'Pupa']


    # Store datasets in a dictionary

DataDictionary = {temp: DevData.iloc[21:1020, cols].values for temp, cols in RearingConditions.items()}


# Define some functions to run the parallel process

    # Function to clean data

def clean_data(data, col_idx):
    return [x for x in data[:, col_idx].tolist() if x != 0 and not np.isnan(x)]


    # Step 2: Define the function to sample and compute the minimum value

def Sampler(N, DataInput):
    if len(DataInput) < N-1:
        return None  # Return None if there are not enough elements to sample
    SubGroup = random.sample(DataInput, N)  # Sample N elements
    return min(SubGroup)  # Compute min value


    # Step 3: Define parallel function

def ParallelSampler(Args):
    return Sampler(N, Args)


    # Step 4: Bootstrap function

def BootstrapCalculator(DataInput, N):

    if len(DataInput) < N:
        return None, None  # Return None values if there are not enough elements
    
    if __name__ == '__main__':
        with multiprocessing.Pool(processes=NumCores) as pool:
            MinDelays = pool.map(ParallelSampler, [DataInput] * BootstrapIterations)
        
        # Filter out None values
        MinDelays = [x for x in MinDelays if x is not None]
        
        if not MinDelays:  # Check again if MinDelays is empty
            return None, None
        
        # Compute mean and standard deviation
        Mean_MinDelay = statistics.mean(MinDelays)
        StdDev_MinDelay = statistics.stdev(MinDelays)
        
        return Mean_MinDelay, StdDev_MinDelay


# Step 5: Iterate over all datasets and compute bootstrap results

if __name__ == '__main__':

    Results = []
    
    for temp, data in DataDictionary.items():
    
        for col_idx, stage in enumerate(stages):
        
            LifeTabData = clean_data(data, col_idx)
            Mean_MinDelay, StdDev_MinDelay = BootstrapCalculator(LifeTabData, N)
            Results.append([temp, stage, Mean_MinDelay if Mean_MinDelay is not None else "", StdDev_MinDelay if StdDev_MinDelay is not None else ""])
    
    # Save results to CSV
    
    ResultsDataFrame = pd.DataFrame(Results, columns=['RearingCondition', 'Stage', 'Mean_MinDelay', 'StdDev_MinDelay'])
    ResultsDataFrame = ResultsDataFrame.sort_values(by = 'Stage')
    ResultsDataFrame.to_csv('BootstrapDevelopment.csv', index = False, sep = ';')
    


## COMPUTE THE PREOVIPOSITION TIME


# Step 1: Load data from Excel

FertData = pd.read_excel('Fertility.xlsx', sheet_name='Fertility')

# Select the columns containing the data for bootstrap

SelectedColumns = [1, 4, 7, 10, 13, 16, 19, 22, 25]

# Assign the names of the columns

ColNamesFert = {1: '13',
                4: '18',
                7: '20',
                10: '24',
                13: '25',
                16: '26',
                19: '27',
                22: '28',
                25: '29'}


# Iterate over all datasets and compute bootstrap results

if __name__ == '__main__':

    ResultsFert = []
    
    for i in SelectedColumns:
    
        # Estract the data from the different columns
        
        FertCol = FertData.iloc[2:12, i].dropna().tolist()
        
        PreOviTime, DevStPreOvi = BootstrapCalculator(FertCol, N_PreOvi)
        
        ResultsFert.append({'Temperature': ColNamesFert.get(i),
                            'Pre-oviposition time': PreOviTime,
                            'Dev.st': DevStPreOvi})

    # Create a dataframe to store results
    
    ResultsPreOvi = pd.DataFrame(ResultsFert)

    # Print results on a CSV
    
    ResultsPreOvi.to_csv('BootstrapFertility.csv', index=False, sep = ';')
    

# Timer stop

    t_f = time.time() - t_start

    print('\nExecution time:', t_f, 's\n')

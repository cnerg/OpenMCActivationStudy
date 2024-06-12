# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 02:05:08 2024

@author: Anupama Rajendra
"""
import csv
import numpy as np
import pandas as pd

#The purpose of this code is to take the flux vs. energy data derived from OpenMC and
# restructure it into a form appropriate for ALARA.

#Loading the csv that contains the energy bounds of the Vitamin J structure 
Vit_J_Groups = pd.read_csv('VitaminJEnergyGroupStructure.csv', header=1)
x_lh = Vit_J_Groups.iloc[:, 1]

x = sorted(x_lh, reverse=True) #Energy bounds arranged from high energy to low energy
        
#Loading the csv file that contains the energy bins and neutron flux values (from OpenMC)
#This csv file was created in a Python code that reads the h5 output from OpenMC
Flux = 'Neutron_Flux.csv'
Flux_Data = pd.read_csv(Flux, header=1)
# Extract both columns from the CSV file
y_lh = Flux_Data.iloc[:, 0]
y = sorted(y_lh, reverse=True) #Energy bins arranged from highest to lowest energy
z_lh = Flux_Data.iloc[:, 1]
z = sorted(z_lh, reverse=True) #Fluxes arranged to highest to lowest corresponding energy
#This section adds up the neutron fluxes between each pair of energy bounds in
#the Vitamin J structure

#Initializing the array of fluxes that are summed up
flux_sums = [0] * (len(x) - 1)

#Iterating over each interval of energies in the Vit J structure
#Since there are 175 energy bins, there are 175 - 1 = 174 intervals

for i in range(len(x) - 1):
    lower_bound = x[i]
    upper_bound = x[i + 1]
    
    # If the OpenMC energy bin corresponding to the flux is equal to energy level of bin 175,
    # this 'if statement' will include the flux in the summation of the last interval
    # (This is done for completeness - it is unlikely that these two energy values will be
    # exactly equal to each other)
    
    if i == len(x) - 2:  # Checking if it's the last interval
        for j in range(len(y)):
            if lower_bound >= y[j] >= upper_bound:  # Including the upper energy bound
                flux_sums[i] += z[j]
                
    # For all other energy intervals, the upper bound is not included in the summation            
    else:
        for j in range(len(y)):
            if lower_bound >= y[j] > upper_bound:
                flux_sums[i] += z[j]
    
#Converting list of fluxes to array format    
flux_sums_array = np.array(flux_sums)

#Saving flux_sums as a 30x6 tab-delimited file:
    
Array_1 = flux_sums_array[:174] #The first 174 elements of the flux array
Array_2 = flux_sums_array[174] #The 175th element of the array
    
flux_sums_shaped = Array_1.reshape(29,6) #Accounts for the first 174 elements

# Open the file in append mode and write each row of the reshaped array
with open('Flux_Summation.txt', 'w') as f:
    for row in flux_sums_shaped:
        formatted_row = ' '.join([f'{value:.5e}' for value in row])
        f.write(formatted_row + '\n')

#Adding the 175th element to the array:
with open('Flux_Summation.txt', 'a') as f:
    f.write(f'{Array_2:.5e}\n')

#print(flux_sums_shaped)
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 02:05:08 2024

@author: Anupama Rajendra
"""
import csv
import numpy as np
import pandas as pd
import openmc

#The purpose of this code is to take the flux vs. energy data derived from OpenMC and
# restructure it into a form appropriate for ALARA.

#Loading the csv that contains the energy bounds of the Vitamin J structure 
Vit_J = pd.read_csv('VitaminJEnergyGroupStructure.csv')
ebounds_lh = Vit_J.iloc[:, 1]

ebounds = sorted(ebounds_lh, reverse=True) #Energy bounds arranged from high energy to low energy
        
#Loading the csv file that contains the energy bins and neutron flux values (from OpenMC)
#This csv file was created in a Python code that reads the h5 output from OpenMC
Flux_Data = pd.read_csv('Neutron_Flux.csv')
# Extract both columns from the CSV file
energy_lh = Flux_Data.iloc[:, 0]
energy = sorted(energy_lh, reverse=True) #Energy bins arranged from highest to lowest energy
flux_lh = Flux_Data.iloc[:, 1]
flux = sorted(flux_lh, reverse=True) #Fluxes arranged to highest to lowest corresponding energy
#This section adds up the neutron fluxes between each pair of energy bounds in
#the Vitamin J structure

#Initializing the array of fluxes that are summed up
flux_sums = np.zeros(175)

results = zip(energy, flux)

#Iterating over each interval of energies in the Vit J structure
#Since there are 175 energy bins, there are 175 - 1 = 174 intervals

groups = len(ebounds) - 1
bin = 0

for energy, flux in results
        while bin < groups and ebounds[bin + 1] >= energy
                bin += 1
        if bin < len(flux_sums):
                flux_sums[bin] += flux

print(flux_sums)

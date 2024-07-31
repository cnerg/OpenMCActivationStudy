# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 12:54:10 2024

@author: Anupama Rajendra
"""

import matplotlib.pyplot as plt
import csv

#Initializing arrays to store data from OpenMC results .csv file:
OpenMC_Nuclides = []
OpenMC_Day = []
OpenMC_Month = []
OpenMC_Shutdown = []

with open('Densities_CSV.csv', 'r') as OpenMC_Data:
    OpenMC_csv_reader = csv.reader(OpenMC_Data)
    for OpenMC_line in OpenMC_csv_reader:
        OpenMC_Nuclides.append(OpenMC_line[0])
        OpenMC_Day.append(OpenMC_line[2])
        OpenMC_Month.append(OpenMC_line[3])
        OpenMC_Shutdown.append(OpenMC_line[4])

# Reading content from ALARA Output file:
with open('ALARA_Data.txt', 'r') as Data_File:
    ALARA_Data = Data_File.readlines()

#Identifying the part of the ALARA output file that contains the relevant density data:
ALARA_FileBounds = []
for Density_Index, Density_Content in enumerate(ALARA_Data):
    if '==' in Density_Content:
        ALARA_FileBounds.append(Density_Index)
        #If '==' is found twice, end the loop  
        if len(ALARA_FileBounds) == 2:
            break
        
#Process the data in between the two lines that begin with '=='
ALARA_Nuclide_Data = ALARA_Data[ALARA_FileBounds[0] + 1:ALARA_FileBounds[1]]

# Stable W nuclides present at beginning of operation:
Stable_Nuc = ['w-180', 'w-182', 'w-183', 'w-184', 'w-186']

#Initializing arrays to store data from ALARA and data common to both ALARA & OpenMC
ALARA_Nuc_no_W = []
ALARA_no_W_Day = []
ALARA_no_W_Month = []
ALARA_no_W_Shutdown = []
ALARA_List = [ALARA_no_W_Day, ALARA_no_W_Month, ALARA_no_W_Shutdown]

Common_Nuclides = []
Common_Day_Diff = []
Common_Month_Diff = []
Common_Shutdown_Diff = []

#Nuclide data found in ALARA only
ALARA_Nuclides_Only = []
ALARA_no_W_Day_Only = []
ALARA_no_W_Month_Only = []
ALARA_no_W_Shutdown_Only = []

for ALARA_Filtered_Lines in ALARA_Nuclide_Data:
    #For any nuclide that is not one of the stable W nuclides:
    if not any(ALARA_Filtered_Lines.startswith(Nuclide) for Nuclide in Stable_Nuc):
        #Formatting change to make OpenMC and ALARA nuclide formats match up
        ALARA_Filtered_Lines = ALARA_Filtered_Lines.replace('-','').capitalize()
        #Storing all density information from ALARA output
        
        #The first part of each line is the nuclide
        Nuc_Each_Line = (ALARA_Filtered_Lines.strip().split()[0])
        ALARA_Nuc_no_W.append(Nuc_Each_Line)
        ALARA_no_W_Day.append(ALARA_Filtered_Lines.strip().split()[1])
        ALARA_no_W_Month.append(ALARA_Filtered_Lines.strip().split()[2])
        ALARA_no_W_Shutdown.append(ALARA_Filtered_Lines.strip().split()[3])
    #Identifying nuclides from ALARA that are also found from OpenMC data:
    if Nuc_Each_Line in OpenMC_Nuclides:
            Common_Nuclides.append(Nuc_Each_Line)
    else:
        #Add to lists that contain nuclides/densities only found in ALARA output
        ALARA_Nuclides_Only.append(Nuc_Each_Line)
        ALARA_no_W_Day_Only.append(ALARA_Filtered_Lines.split()[1])
        ALARA_no_W_Month_Only.append(ALARA_Filtered_Lines.split()[2])
        ALARA_no_W_Shutdown_Only.append(ALARA_Filtered_Lines.split()[3])

# Times after shutdown [s]
time_steps = [0, 86400, 2.6864e+06]

# Plot data for each isotope found in ALARA
for i, isotope in enumerate(ALARA_Nuc_no_W):
    densities = [ALARA_no_W_Day[i], ALARA_no_W_Month[i], ALARA_no_W_Shutdown[i]]
    plt.plot(time_steps, densities, marker='o', label=isotope)

plt.xlabel('Time (seconds)')
plt.ylabel('Number Density (atoms/cm^3)')
plt.title('Number Density [atoms/cm^3] vs. Time After Shutdown [s]')
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.subplots_adjust(right=0.7)
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.savefig('Nuclide_density_ALARA')
plt.show()

#Initializing arrays to store ratios of nuclide densities from both codes
Ratio_Day = []
Ratio_Month = []
Ratio_Shutdown = []

#print("All nuclides found in OpenMC", OpenMC_Nuclides)
print("Nuclides common to ALARA and OpenMC", Common_Nuclides)
#print("Nuclides found only in ALARA:", ALARA_Nuclides_Only)

#Iterating over all nuclides common to both sets of results
for Common_Name in Common_Nuclides:
    #Identifying the location of each common nuclide in the OpenMC nuclide list
    Index_OpenMC = OpenMC_Nuclides.index(Common_Name)
    Density_OpenMC_Day = float(OpenMC_Day[Index_OpenMC])
    Density_OpenMC_Month = float(OpenMC_Month[Index_OpenMC])
    Density_OpenMC_Shutdown = float(OpenMC_Shutdown[Index_OpenMC])
    
    #Identifying the location of each common nuclide in the ALARA nuclide list
    Index_ALARA  = ALARA_Nuc_no_W.index(Common_Name)
    Density_ALARA_Day = float(ALARA_no_W_Day[Index_ALARA])
    Density_ALARA_Month = float(ALARA_no_W_Month[Index_ALARA]) 
    Density_ALARA_Shutdown = float(ALARA_no_W_Shutdown[Index_ALARA])
    
    #Ratios between the values from ALARA and OpenMC:
    Ratio_Day.append(Density_OpenMC_Day / Density_ALARA_Day)
    Ratio_Month.append(Density_OpenMC_Month / Density_ALARA_Month)
    Ratio_Shutdown.append(Density_OpenMC_Shutdown / Density_ALARA_Shutdown)  

    #Absolute values of the differences between ALARA and OpenMC
    Common_Day_Diff.append(abs(Density_OpenMC_Day - Density_ALARA_Day))
    Common_Month_Diff.append(abs(Density_OpenMC_Month - Density_ALARA_Month))
    Common_Shutdown_Diff.append(abs(Density_OpenMC_Shutdown - Density_ALARA_Shutdown))   

for j, com_iso in enumerate(Common_Nuclides):
    common_diff = [Common_Day_Diff[j], Common_Month_Diff[j], Common_Shutdown_Diff[j]]
    common_ratios = [Ratio_Day[j], Ratio_Month[j], Ratio_Shutdown[j]]
    #plt.plot(time_steps, common_diff, marker='o', label=com_iso)
    plt.plot(time_steps, common_ratios, marker = 'o', label = com_iso) 
    
# plt.xlabel('Time [seconds]')
# plt.xlim(1, 1E+7)
# #plt.ylim(1E+12, 1E+15)
# plt.ylabel('Number Density Difference [atoms/cm^3]')
# plt.title('Number Density [atoms/cm^3] vs. Time After Shutdown [s]')
# plt.subplots_adjust(right=0.7)
# plt.xscale('log')
# plt.yscale('log')
# plt.grid(True)
#plt.savefig('Nuclide_density_diff')
# plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
# plt.show()      
    
plt.xlabel('Time [seconds]')
plt.ylabel('Ratio Between OpenMC and ALARA')
plt.title('Number Density Ratio vs. Time After Shutdown [s]')
plt.subplots_adjust(right=0.7)
plt.legend(loc='upper left', bbox_to_anchor=(1.15, 1.1))
plt.grid(True)
plt.savefig('Nuclide_density_ratio')
plt.show()    
plt.close()
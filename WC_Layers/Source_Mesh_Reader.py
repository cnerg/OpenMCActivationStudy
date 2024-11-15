# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:16:44 2024

@author: Anupama Rajendra
"""

import h5py
import numpy as np

# Open the HDF5 files
File_1 = h5py.File("source_mesh_1.h5m", 'r')
File_2 = h5py.File("source_mesh_2.h5m", 'r')
Files = [File_1, File_2]

# Number of photon groups
photon_groups = 24

def extract_source_data(Files):
    '''
    Identifies the location of the source density dataset within each mesh file.
    '''
    sd_list = []
    for file in Files:
        # Extract data from the current file
        tstt = file['tstt']
        elements = tstt['elements']
        tet_four = elements['Tet4']
        tags = tet_four['tags']
        sd = tags['source_density'][:]
        sd_list.append(sd)
    return sd_list
    
def arrange_source_strengths(sd_list, photon_groups):
    '''
    Writes individual mesh element strengths (separated by energy group) to density_list, 
    then writes density_list to all_strengths_list for each set of source densities. 
    '''
    all_strengths_list = []
    for density in sd_list:
        density_list = []
        # Calculate (individual mesh) strengths for each photon group
        for group in range(photon_groups):
            strengths = density[:,group]
            density_list.append(strengths)
        all_strengths_list.append(density_list)    
    return all_strengths_list    

def write_source_density(sd_list, sd_filename):
    '''
    Writes each of the specified source density datasets to a text file.
    
    sd_list: list of source density datasets from h5m files
    sd_filename: user-provided filename that appears as a prefix for the text files
    '''
    for sd_index, source_density in enumerate(sd_list):
        # Write source density to a separate file for each mesh
        with open(f'{sd_filename}_{sd_index + 1}.txt', 'w') as source:
            for tet_element in source_density:
                source.write(' '.join(map(str, tet_element)) + '\n')
                
def extract_arrange_write_sd(Files, photon_groups, sd_filename):
    data_extractor = extract_source_data(Files)
    strength_arranger = arrange_source_strengths(data_extractor, photon_groups)
    write_source_density(data_extractor, sd_filename)
    return data_extractor, strength_arranger

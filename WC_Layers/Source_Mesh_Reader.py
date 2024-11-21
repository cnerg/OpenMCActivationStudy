# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:16:44 2024

@author: Anupama Rajendra
"""

import h5py
import numpy as np

source_density_name = 'source_density'

# Open the HDF5 files
File_1 = h5py.File("source_mesh_1.h5m", 'r')
File_2 = h5py.File("source_mesh_2.h5m", 'r')
Files = [File_1, File_2]

def extract_source_data(file_list):
    '''
    Identifies the location of the source density dataset within each mesh file.
    
    input: list of .h5m mesh files (opened by h5py) containing photon source information
    output: numpy array of source density data with rows = # of mesh elements and 
    columns = # number of photon groups, with one array per source mesh
    '''
    sd_list = np.ndarray((len(file_list), num_elements, photon_groups))

    for file_num, file in enumerate(file_list):
        sd_list[file_num,:] = file['tstt']['elements']['Tet4']['tags']['source_density'][:]
    return sd_list  

def save_source_density(sd_list, sd_filename):
    '''
    Saves source density data as a text file.
    
    sd_list: list of source density datasets from h5m files
    sd_filename: user-provided filename that appears as a prefix for the text files
    '''
    for sd_index, source_density in enumerate(sd_list):
        with open(f'{sd_filename}_{sd_index + 1}.txt', 'w') as source:
            for tet_element in source_density:
                source.write(' '.join(map(str, tet_element)) + '\n')
                
def extract_save_sd(file_list, sd_filename):
    data_extractor = extract_source_data(file_list)
    save_source_density(data_extractor, sd_filename)
    return data_extractor

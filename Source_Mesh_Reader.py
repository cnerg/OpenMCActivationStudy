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

# Process each file
for file_index, file in enumerate(Files):
    # Extract data from the current file
    tstt = file['tstt']
    elements = tstt['elements']
    tet_four = elements['Tet4']
    tags = tet_four['tags']
    sd = tags['source_density'][:]
    
    # Write source density to a separate file for each mesh
    with open(f'source_density_{file_index + 1}.txt', 'w') as source:
        for tet_element in sd:
            source.write(' '.join(map(str, tet_element)) + '\n')

    # Calculate summed strengths for each photon group
    summed_strengths = []
    for group in range(photon_groups):
        strengths = np.sum(sd[:, group])  # Sum over all mesh elements
        summed_strengths.append(strengths)

    # Write summed strengths to a separate file for each mesh
    with open(f'summed_photon_strengths_{file_index + 1}.txt', 'w') as sums:
        sums.write(' '.join(map(str, summed_strengths)) + '\n')
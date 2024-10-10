# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 23:18:20 2024

@author: Anupama Rajendra
"""
import openmc
import matplotlib.pyplot as plt

with openmc.StatePoint("statepoint.10.h5") as sp:
    flux_spectrum = sp.get_tally(id=2) 
    mesh=sp.meshes[1]
    # Get the reshaped tally data
    tally_data_reshaped = flux_spectrum.get_reshaped_data(value='mean')

    # Print the shape of the tally data
    print("Tally data shape:", tally_data_reshaped.shape)

    flux_sum_en = tally_data_reshaped.sum(axis=(0,1,3,4,5))   
    #Vitamin-J energy filter:
    e_filter = flux_spectrum.filters[2]
    #Lower bounds of the energy bins
    e_filter_lower = e_filter.bins[:, 0]

    # Plot flux spectrum
    fix, ax = plt.subplots()
    ax.loglog(e_filter_lower, flux_sum_en, drawstyle='steps-post')
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Flux [photon-cm/source]')
    ax.grid(True, which='both')
    plt.savefig('Photon_flux_vs_energy.png')
    plt.show()
    
    mesh_data = tally_data_reshaped.sum(axis=(0,2,3,4,5))
    vtk = mesh.write_data_to_vtk(filename="Photon_Flux.vtk", datasets={"mean":mesh_data.flatten()})
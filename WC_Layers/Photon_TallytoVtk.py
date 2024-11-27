# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 23:18:20 2024

@author: Anupama Rajendra
"""
import openmc
import matplotlib.pyplot as plt

def read_statepoint(sp_filename, flux_spectrum_tally_id, mesh_number, summed_axes_energy_filter, energy_filter_index):
    with openmc.StatePoint(sp_filename) as sp:
        flux_spectrum = sp.get_tally(id=flux_spectrum_tally_id) 
        mesh=sp.meshes[mesh_number]
        # Return tally data condensed into 1 dimension per filter
        tally_data_reshaped = flux_spectrum.get_reshaped_data(value='mean')
        flux_sum_energy = tally_data_reshaped.sum(axis=summed_axes_energy_filter)   
        #Vitamin-J energy filter:
        e_filter = flux_spectrum.filters[energy_filter_index]
        #Lower bounds of the energy bins
        e_filter_lower = e_filter.bins[:, 0]
    return e_filter_lower, flux_sum_energy, tally_data_reshaped, mesh

def plot_flux_data(e_filter_lower, flux_sum_energy, figure_filename):
    # Plot flux spectrum
    fix, ax = plt.subplots()
    ax.loglog(e_filter_lower, flux_sum_energy, drawstyle='steps-post')
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Flux [photon-cm/source]')
    ax.grid(True, which='both')
    plt.savefig(figure_filename)
    plt.show()

def write_to_vtk(tally_data_reshaped, summed_axes_mesh_filter, vtk_filename, mesh):
    mesh_data = tally_data_reshaped.sum(axis=summed_axes_mesh_filter)
    mesh.write_data_to_vtk(filename=vtk_filename, datasets={"mean":mesh_data.flatten()})

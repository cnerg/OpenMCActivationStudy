import openmc
import matplotlib.pyplot as plt
import yaml
import argparse
import numpy as np

def read_statepoint(sp_filename, flux_spectrum_tally_id, summed_axes_energy_filter, energy_filter_index, mesh_number):
    '''
    Reads OpenMC Statepoint file and returns energy bins and flux from reshaped data.
    
    inputs: 
        sp_filename : path to OpenMC Statepoint file
        flux_spectrum_tally_id : id of flux tally with energy filter
        summed_axes_energy_filter: n-tuple of axes over which tally data is summed in order to isolate axis with energy bounds
        energy_filter_index: index of energy filter within flux tally
        mesh_number: mesh id of OpenMC unstructued mesh
        
    '''
    with openmc.StatePoint(sp_filename) as sp:
        flux_spectrum_tally = sp.get_tally(id=flux_spectrum_tally_id) 
        mesh=sp.meshes[mesh_number]
        # Return tally data condensed into 1 dimension per filter
        flux_spectrum_mean = flux_spectrum_tally.get_reshaped_data(value='mean')
        total_flux = flux_spectrum_mean.sum(axis=eval(summed_axes_energy_filter))   
        #Vitamin-J energy filter:
        e_filter = flux_spectrum_tally.filters[energy_filter_index]
        #Lower bounds of the energy bins
        e_filter_lower = e_filter.bins[:, 0]
    return e_filter_lower, total_flux, flux_spectrum_mean, mesh

def plot_flux_data(e_filter_lower, total_flux, figure_filename):
    '''
    Plots flux data as a function of energy.
    
    inputs:
        e_filter_lower: iterable of lower bounds of each energy bin
        total_flux: 
        figure_filename: file name of the plot of flux vs energy    
    '''
    fix, ax = plt.subplots()
    ax.loglog(e_filter_lower, total_flux, drawstyle='steps-post')
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Flux [photon-cm/source]')
    ax.grid(True, which='both')
    plt.savefig(figure_filename)
    plt.show()

def calculate_dose(flux_spectrum_mean, summed_axes_dose_filter):
    '''
    Calculates effective dose [pSv-cm^3/source] based on energy function filter
    inputs:
        flux_spectrum_mean: reshaped flux spectrum data from read_statepoint
        summed_axes_dose_filter: n-tuple of axes over which tally data is summed in order to isolate dose filter axis
    '''    
    total_dose = flux_spectrum_mean.sum(axis=eval(summed_axes_dose_filter))
    return total_dose

def write_to_vtk(flux_spectrum_mean, summed_axes_mesh_filter, vtk_filename, mesh):
    '''
    Saves mesh tally data in vtk format.
    inputs:
        flux_spectrum_mean: reshaped flux spectrum data from read_statepoint
        summed_axes_mesh_filter: n-tuple of axes over which tally data is summed in order to isolate mesh axis
        vtk_filename: name of file saved in vtk format
        mesh: OpenMC Unstructured Mesh object
    '''
    mesh_data = flux_spectrum_mean.sum(axis=eval(summed_axes_mesh_filter))
    mesh.write_data_to_vtk(filename=vtk_filename, datasets={"mean":mesh_data.flatten()})

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--Photon_Transport_YAML', default = "PhotonTransport_Inputs.yaml", help="Path (str) to YAML containing inputs for Photon_TallytoVtk")
    args = parser.parse_args()
    transport_yaml = args.Photon_Transport_YAML
    with open(transport_yaml, 'r') as file:
        inputs = yaml.safe_load(file)
    sp_filename = inputs['filename_dict']['sp_filename']
    flux_spectrum_tally_id = inputs['file_indices']['flux_spectrum_tally_id']
    summed_axes_energy_filter = inputs['summed_axes_energy_filter']
    energy_filter_index = inputs['file_indices']['energy_filter_index']
    mesh_number = inputs['file_indices']['mesh_number']
    rs = read_statepoint(sp_filename, flux_spectrum_tally_id, summed_axes_energy_filter, energy_filter_index, mesh_number)
    figure_filename = inputs['filename_dict']['figure_filename']
    pfd = plot_flux_data(rs[0], rs[1], figure_filename)
    
    summed_axes_dose_filter = inputs['summed_axes_dose_filter']
    cd = calculate_dose(rs[2], summed_axes_dose_filter)
    
    summed_axes_mesh_filter = inputs['summed_axes_mesh_filter']
    vtk_filename = inputs['filename_dict']['vtk_filename']
    wtv = write_to_vtk(rs[2], summed_axes_mesh_filter, vtk_filename, rs[3])

if __name__ == "__main__":
    main()

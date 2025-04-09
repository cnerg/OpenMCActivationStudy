import openmc
import matplotlib.pyplot as plt
import yaml
import argparse
import numpy as np

def read_statepoint(sp_filename, flux_spectrum_tally_id, axes_without_energy_bins, energy_filter_index, mesh_number):
    '''
    Reads OpenMC Statepoint file and returns energy bins and flux from reshaped data.
    
    inputs: 
        sp_filename : path to OpenMC Statepoint file
        flux_spectrum_tally_id : id of flux tally with energy filter
        axes_without_energy_bins: n-tuple of axes over which tally data is summed in order to isolate axis with energy bounds
        energy_filter_index: index of energy filter within flux tally
        mesh_number: mesh id of OpenMC unstructured mesh
        
    '''
    with openmc.StatePoint(sp_filename) as sp:
        photon_tally = sp.get_tally(id = photon_tally_id)
        phtn_tally_e_filter_lower = photon_tally.filters[0].bins[:,0]
        
        flux_spectrum_tally = sp.get_tally(id=flux_spectrum_tally_id) 
        mesh=sp.meshes[mesh_number]
        # Return tally data condensed into 1 dimension per filter
        flux_spectrum_mean = flux_spectrum_tally.get_reshaped_data(value='mean')
        energy_binned_flux = flux_spectrum_mean.sum(axis=eval(axes_without_energy_bins))   
        #Vitamin-J energy filter:
        e_filter = flux_spectrum_tally.filters[energy_filter_index]
        #Lower bounds of the energy bins
        e_filter_lower = e_filter.bins[:, 0]
    return photon_tally, phtn_tally_e_filter_lower, e_filter_lower, energy_binned_flux, flux_spectrum_mean, mesh

def plot_photon_tally(photon_tally, phtn_tally_e_filter_lower,  photon_tally_figname):
    '''
    Plots flux and absorption tallies as a function of energy bounds from ALARA
    inputs:
        photon_tally : OpenMC Tally object with applied energy filter
        phtn_e_filter_lower : iterable of lower bounds of each energy bin
        photon_tally_figname : file name of the plot of flux & absorption vs energy
    '''
    fix, ax = plt.subplots()
    for score in photon_tally.scores:
        score_slice = photon_tally.get_slice(scores=[score])
        score_slice_values = score_slice.get_values(value='mean').ravel()
        plt.loglog(phtn_tally_e_filter_lower, score_slice_values, label=score)
    plt.legend()
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Tally Value [photon-cm/source]')
    ax.grid(True, which='both')    
    plt.savefig(photon_tally_figname) 

def plot_flux_data(e_filter_lower, energy_binned_flux, figure_filename):
    '''
    Plots flux data as a function of energy.
    
    inputs:
        e_filter_lower: iterable of lower bounds of each energy bin
        energy_binned_flux: photon flux spectrum (summed over any other dimensions)
        figure_filename: file name of the plot of flux vs energy    
    '''
    fix, ax = plt.subplots()
    ax.loglog(e_filter_lower, energy_binned_flux, drawstyle='steps-post')
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Flux [photon-cm/source]')
    ax.grid(True, which='both')
    plt.savefig(figure_filename)
    plt.show()

def save_summed_data_to_vtk(flux_spectrum_mean, axes_without_mesh, vtk_filename, mesh):
    '''
    Saves mesh tally data in vtk format.
    inputs:
        flux_spectrum_mean: reshaped flux spectrum data from read_statepoint
        axes_without_mesh: n-tuple of axes over which tally data is summed in order to isolate axis containing mesh elements
        vtk_filename: name of file saved in vtk format
        mesh: OpenMC Unstructured Mesh object
    '''
    mesh_data = flux_spectrum_mean.sum(axis=eval(axes_without_mesh))
    mesh.write_data_to_vtk(filename=vtk_filename, datasets={"mean":mesh_data.flatten()})

def main():
    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument('--Photon_Transport_YAML', default = "PhotonTransport_Inputs.yaml", help="Path (str) to YAML containing inputs for Photon_TallytoVtk")
        args = parser.parse_args()
        return args

    def read_yaml(args):
        with open(args.Photon_Transport_YAML, 'r') as file:
            inputs = yaml.safe_load(file)
        return inputs
    
    def save_photon_tally_vtk(inputs):
       photon_tally, phtn_tally_e_filter_lower, e_filter_lower, energy_binned_flux, flux_spectrum_mean, mesh = read_statepoint(inputs['filename_dict']['sp_filename'], 
                                                                              inputs['file_indices']['photon_tally_id'],
                                                                               inputs['file_indices']['flux_spectrum_tally_id'],
                                                                               inputs['axes_without_energy_bins'], 
                                                                               inputs['file_indices']['energy_filter_index'],
                                                                               inputs['file_indices']['mesh_number'])
        photon_tally_plot = plot_photon_tally(photon_tally, phtn_tally_e_filter_lower,
                                              inputs['filename_dict']['photon_tally_figname'])
        flux_data_plot = plot_flux_data(e_filter_lower, energy_binned_flux,
                                        inputs['filename_dict']['figure_filename'])

        total_dose = flux_spectrum_mean.sum(axis=eval(inputs['axes_without_dose_filter_bin']))
        
        vtk_file = save_summed_data_to_vtk(flux_spectrum_mean, 
                          inputs['axes_without_mesh'], 
                          inputs['filename_dict']['vtk_filename'], 
                          mesh)
        
    args = parse_args()
    inputs = read_yaml(args)
    save_photon_tally_vtk(inputs)

if __name__ == "__main__":
    main()

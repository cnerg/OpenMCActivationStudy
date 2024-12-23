import openmc
import numpy as np
import matplotlib.pyplot as plt

def extract_tally_values(statepoint_file_path) :
    '''
    Stores values from each tally in a numpy array
    
    inputs:
        statepoint_file_path : path to OpenMC Statepoint .h5 file
    '''
    with openmc.StatePoint(statepoint_file_path) as sp:
        tallies = sp.tallies
        tally_array = []
        for tally_id, tally in tallies.items():
            tally_values = tally.get_values(value="mean").ravel()
            tally_array.append(np.array(tally_values))
        tally_array = np.array(tally_array)   
    return tallies, tally_array
    
def plot_flux_tally(tallies, flux_tally_id, energy_filter_index) :
    '''
    Plots flux tally as a function of energy
    
    inputs :
        tallies : dictionary where keys = OpenMC Tally IDs and values = OpenMC Tally objects
        flux_tally_id : id of OpenMC Tally object that scores flux with energy filter
        energy_filter_index : index of EnergyFilter within the list of applied filters
    '''
    flux_tally = tallies[flux_tally_id]
    flux_tally_values = flux_tally.get_values(value = 'mean').ravel()
    energy_filter = flux_tally.filters[energy_filter_index]
    energy_bins = energy_filter.bins[:, 0]
    fix, ax = plt.subplots()
    ax.loglog(energy_bins, flux_tally_values, drawstyle='steps-post')
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Flux [neutron-cm/source]')
    ax.grid(True, which='both')
    plt.savefig('Flux_spectrum_vs_energy.png')
    plt.close()
    return flux_tally_values, energy_bins

def save_data(flux_tally_values, mat_vol, tally_array) :
    '''
    Saves average values of all tallies and neutron flux in a format appropriate for ALARA
    
    inputs :
        flux_tally_values : numpy array of neutron flux values over energy filter
        mat_vol : volume of tallied material (float)
        tally_array : numpy array containing separate numpy arrays for each tally
    '''
    np.savetxt('tally_averages', tally_array, fmt='%s')
    flux = flux_tally_values/mat_vol
    #ALARA flux inputs go from high energy to low energy
    flux_alara = flux[::-1]
    np.savetxt('flux_for_alara', flux_alara)

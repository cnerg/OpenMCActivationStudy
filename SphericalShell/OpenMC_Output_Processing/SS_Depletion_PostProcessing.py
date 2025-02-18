import openmc
import openmc.deplete
import numpy as np
import matplotlib.pyplot as plt

def extract_nuclides(dep_file_path, time_units, depletable_mat_index) : 
    '''
    Identifies nuclides from each depletion timestep and stores them in an array
    
    inputs: 
        dep_file_path : path to .h5 file with OpenMC depletion simulation results
        time_units : units ('s', 'd', 'min', 'h', 'a') used when retrieving depletion timesteps
        depletable_mat_index : index of depletable material in OpenMC Materials object
    '''    
    dep_results = openmc.deplete.Results(filename=dep_file_path)
    time_steps = dep_results.get_times(time_units = time_units)
    nuclide_set = set()
    stable_init_nuc = set()
    for step in range(len(time_steps)):
        # Obtain a materials object with depletion information at each timestep
        materials_object = dep_results.export_to_materials(step)[depletable_mat_index]
        for nuclide in materials_object.get_nuclides() :
            if step == 0:
                half_life = openmc.data.half_life(nuclide)
                if half_life == None:
                    stable_init_nuc.add(nuclide)
            nuclide_set.add(nuclide)     
    nuclide_set = nuclide_set - stable_init_nuc
    return nuclide_set, materials_object, dep_results, time_steps

def extract_data(nuclide_set, materials_object, dep_results, time_steps, nuc_units):
    '''
    Extract and store nuclide density data to be accessed in other functions.
    
    inputs :
        nuclide_set : iterable of nuclide names (str)
        nuc_units : units of nuclide concentration ('atoms', 'atom/b-cm', 'atom/cm3')
        (All other inputs from the output of extract_nuclides())
    '''
    num_dens = {}
    for nuclide in nuclide_set:    
        times, num_dens[nuclide] = dep_results.get_atoms(materials_object, nuclide, nuc_units = nuc_units)
    return times, num_dens

def plot_data(times, num_dens, nuclide_set):
    '''
    Plots nuclide density vs. time.      
    '''
    for nuclide in nuclide_set:
        plt.plot(times, num_dens[nuclide], marker='.', linestyle='solid', label=nuclide)

def save_dep_data(times, num_dens, nuclide_set):
    '''
    Saves nuclide density vs. time data to text file.
    '''
    with open(r'number_density_vs_time.txt', 'w') as density_file:
        for nuclide in nuclide_set:
            density_file.write(f'{nuclide} : ' + '\n')
            for time_step, density in zip(times, num_dens[nuclide]):
                density_file.write(f'  {time_step} : {density}\n')
    plt.xlabel('Time after beginning of operation [s]')
    plt.ylabel('Nuclide density [atoms/cm^3]')
    plt.xscale("log")
    plt.yscale("log")
    plt.title('Plot of number density vs time')
    plt.legend()
    plt.savefig('Nuclide_density_OpenMC')
    plt.close()

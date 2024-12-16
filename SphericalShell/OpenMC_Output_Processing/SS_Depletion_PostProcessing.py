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
    nuclide_array = np.ndarray([], dtype = object)
    stable_init_nuc = np.ndarray([], dtype = object)
    
    for step in range(len(time_steps)):
        # Obtain a materials object with depletion information at each timestep
        materials_object = dep_results.export_to_materials(step)[depletable_mat_index]
        all_nuc = materials_object.get_nuclides()      
        for nuclide in all_nuc :   
            # Remove stable nuclides at beginning of operation
            if step == 0:
                half_life = openmc.data.half_life(nuclide)
                if half_life == None:
                    stable_init_nuc = np.append(stable_init_nuc, nuclide)
            if nuclide not in nuclide_array and nuclide not in stable_init_nuc:        
                nuclide_array = np.append(nuclide_array, nuclide) 
                       
    # Remove None element from array initialization    
    nuclide_array = [nuc for nuc in nuclide_array if nuc is not None]
    return nuclide_array, materials_object, dep_results, time_steps

def plot_save_data(nuclide_array, materials_object, dep_results, time_steps, nuc_units) :
    '''
    Plot nuclide density vs time data, and save data to a text file
    
    inputs :
        nuclide_array : iterable of nuclide names (str)
        nuc_units : units of nuclide concentration ('atoms', 'atom/b-cm', 'atom/cm3')
        (All other inputs from the output of extract_nuclides())    
    '''
    with open(r'number_density_vs_time.txt', 'w') as density_file:
        num_dens = {}
        for nuclide in nuclide_array : 
            times, num_dens[nuclide] = dep_results.get_atoms(materials_object, nuclide, nuc_units = nuc_units)
            plt.plot(times, num_dens[nuclide], marker='.', linestyle='solid', label=nuclide)
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

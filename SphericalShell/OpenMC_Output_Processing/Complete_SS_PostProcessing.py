import openmc
import openmc.deplete
import numpy as np
import matplotlib.pyplot as plt
import argparse
import yaml

#Read tally data from statepoint:

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
    
def plot_flux_spectrum(flux_tally, energy_filter_index) :
    '''
    Plots flux tally as a function of energy
    
    inputs :
        flux_tally : OpenMC Tally object that scores flux with energy filter
        energy_filter_index : index of EnergyFilter within the list of applied filters
    '''
    flux_tally_values = flux_tally.get_values(value = 'mean').ravel()
    energy_bins = flux_tally.filters[energy_filter_index].bins[:, 0]
    fix, ax = plt.subplots()
    ax.loglog(energy_bins, flux_tally_values, drawstyle='steps-post')
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Flux [neutron-cm/source]')
    ax.grid(True, which='both')
    plt.savefig('Flux_spectrum_vs_energy.png')
    plt.close()
    return flux_tally_values, energy_bins

def save_tally_data(flux_tally_values, tally_array, inner_radius, thickness) :
    '''
    Saves average values of all tallies and neutron flux in a format appropriate for ALARA
    
    inputs :
        flux_tally_values : numpy array of neutron flux values over energy filter
        tally_array : numpy array containing separate numpy arrays for each tally
        inner_radius : inner radius (float) of material 
        thickness: radial thickness (float) of material
    '''
    np.savetxt('tally_averages', tally_array, fmt='%s')
    mat_vol = 4.0/3.0 * np.pi * ((inner_radius+thickness)**3 - inner_radius**3)
    flux = flux_tally_values/mat_vol
    #ALARA flux inputs go from high energy to low energy
    flux_alara = flux[::-1]
    np.savetxt('flux_for_alara', flux_alara)
    
#Read depletion data from depletion results    
    
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

def extract_dep_data(nuclide_set, materials_object, dep_results, time_steps, nuc_units): #**
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

def plot_dep_data(times, num_dens, nuclide_set):
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
    
# Define inputs and execute all functions:    

def parse_args():
   parser = argparse.ArgumentParser()
   parser.add_argument('--yaml_postprocess_path', default = 'SS_Post_Processing_YAML.yaml', help="Path to YAML file containing inputs for post-processing (str)")
   parser.add_argument('--yaml_model_path', default = 'OpenMC_SS_YAML.yaml', help="Path to YAML file containing required inputs to build model (str)")
   args = parser.parse_args()
   return args

def read_yamls(args):
    with open(args.yaml_postprocess_path, 'r') as yaml_pp :
        pp_inputs = yaml.safe_load(yaml_pp)
    with open(args.yaml_model_path, 'r') as model_file:
        model_inputs = yaml.safe_load(model_file)
    return pp_inputs, model_inputs
    
def post_process_tallies(pp_inputs, model_inputs):
    filepaths = pp_inputs['filepaths']  
    indices = pp_inputs['indices']
    geom_info = model_inputs['geom_info']
        
    tallies, tally_array = extract_tally_values(pp_inputs['filepaths']['statepoint_file_path'])
        
    flux_tally = tallies[indices['flux_tally_id']]
    flux_tally_values, energy_bins = plot_flux_spectrum(flux_tally, 
                      indices['energy_filter_index']) 
        
    tally_averages =  save_tally_data(flux_tally_values, tally_array, 
                                    geom_info['inner_radius'], 
                                    geom_info['thickness'])

def post_process_dep(pp_inputs):
    nuclide_set, materials_object, dep_results, time_steps = extract_nuclides(pp_inputs['filepaths']['dep_file_path'], 
                      pp_inputs['units']['time_units'], 
                      pp_inputs['indices']['depletable_mat_index'])  
    times, num_dens = extract_dep_data(nuclide_set, materials_object, dep_results, time_steps, 
                      pp_inputs['units']['nuc_units'])
    dep_plot = plot_dep_data(times, num_dens, nuclide_set)
    dep_data = save_dep_data(times, num_dens, nuclide_set)

def main() : 
    args = parse_args()
    pp_inputs, model_inputs = read_yamls(args)
    post_process_tallies(pp_inputs, model_inputs)
    post_process_dep(pp_inputs)    

if __name__ == "__main__":
    main()
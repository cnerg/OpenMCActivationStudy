import openmc
import openmc.deplete
import yaml
import argparse
from SS_Tally_PostProcessing import extract_tally_values, plot_flux_spectrum, save_tally_data
from SS_Depletion_PostProcessing import extract_nuclides, extract_dep_data, plot_dep_data, save_dep_data

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
    units = pp_inputs['units']        
    indices = pp_inputs['indices']
    geom_info = model_inputs['geom_info']
        
    tallies, tally_array = extract_tally_values(pp_inputs['filepaths']['statepoint_file_path'])
        
    flux_tally = tallies[pp_inputs['indices']['flux_tally_id']]
    flux_tally_values, energy_bins = plot_flux_spectrum(flux_tally, 
                      pp_inputs['indices']['energy_filter_index']) 
        
    tally_averages =  save_tally_data(flux_tally_values, tally_array, 
                                    model_inputs['geom_info']['inner_radius'], 
                                    model_inputs['geom_info']['thickness'])

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

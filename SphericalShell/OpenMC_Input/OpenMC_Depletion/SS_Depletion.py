import openmc
import openmc.deplete
import numpy as np
import argparse
import yaml

def deplete_ss(chain_file_path, model, inner_radius, thickness, times_post_boc, particle_source_rates, norm_mode, timestep_units):
    chain_file = chain_file_path
    
    material = model.materials[0]
    material.depletable = True
    material.volume = 4.0/3.0 * np.pi * ((inner_radius+thickness)**3 - inner_radius**3)
    timesteps = times_post_boc
    source_rates = particle_source_rates
    operator = openmc.deplete.CoupledOperator(model, chain_file, normalization_mode = norm_mode)
    integrator = openmc.deplete.PredictorIntegrator(operator, timesteps, source_rates = source_rates, timestep_units = timestep_units)
    return integrator

def main():
    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument('--yaml_filepath', default = 'OpenMC_SS_YAML.yaml', help="Path to YAML file containing required inputs (str)")
        args = parser.parse_args()
        return args
    
    def read_yaml(args):   
        with open(args.yaml_filepath, 'r') as yaml_dep:
            dep_inputs = yaml.safe_load(yaml_dep)    
        return dep_inputs
    
    def run_depletion(dep_inputs):
        geom_info = dep_inputs['geom_info']
        dep_params = dep_inputs['depletion_params']    
        model_file = dep_params['model_file']
        model = openmc.model.Model.from_model_xml(path=model_file)
    
        integrator = deplete_ss(dep_params['chain_file'],
                     model, 
                     geom_info['inner_radius'], 
                     geom_info['thickness'], 
                     dep_params['times_post_boc'], 
                     dep_params['source_rates'],
                     dep_params['norm_mode'], 
                     dep_params['timestep_units'])     
    
        integrator.integrate()
        
    args = parse_args()
    dep_inputs = read_yaml(args)
    depletion = run_depletion(dep_inputs)
    
if __name__ == "__main__":
    main()

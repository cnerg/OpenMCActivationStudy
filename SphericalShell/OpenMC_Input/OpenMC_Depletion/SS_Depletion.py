import openmc
import openmc.deplete
import numpy as np
import argparse
import yaml

def deplete_ss(chain_file_path, model, inner_radius, thickness, times_post_boc, particle_source_rates, norm_mode, timestep_units):
    material = model.materials[0]
    chain_file = chain_file_path
    material.depletable = True
    timesteps = times_post_boc
    source_rates = particle_source_rates
    operator = openmc.deplete.CoupledOperator(model, chain_file, normalization_mode = norm_mode)
    integrator = openmc.deplete.PredictorIntegrator(operator, timesteps, source_rates = source_rates, timestep_units = timestep_units)
    integrator.integrate()
    return integrator

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--yaml_filepath', default = 'OpenMC_SS_YAML.yaml', help="Path to YAML file containing required inputs (str)")
    args = parser.parse_args()
    yaml_filepath = args.yaml_filepath
    with open(yaml_filepath, 'r') as yaml_dep:
        dep_inputs = yaml.safe_load(yaml_dep)    
    geom_info = dep_inputs['geom_info']
    dep_params = dep_inputs['depletion_params']    
        
    chain_file_path = dep_params['chain_file']
    model_file = dep_params['model_file']
    model = openmc.model.Model.from_model_xml(path=model_file)
    inner_radius = geom_info['inner_radius']
    thickness = geom_info['thickness']
    times_post_boc = dep_params['times_post_boc']
    particle_source_rates = dep_params['source_rates']
    norm_mode = dep_params['norm_mode']
    timestep_units = dep_params['timestep_units'] 
    
    dss = deplete_ss(chain_file_path, model, inner_radius, thickness, times_post_boc, particle_source_rates, norm_mode, timestep_units)     
    dss.integrate()
    
if __name__ == "__main__":
    main()

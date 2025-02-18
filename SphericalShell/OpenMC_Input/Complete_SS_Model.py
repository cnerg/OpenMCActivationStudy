import openmc
import openmc.deplete
import numpy as np
import yaml
import argparse

def alara_element_densities(elelib_fp):  
    '''
    Create a dictionary of element names and their corresponding densities using the ALARA element library.
    
    inputs:
        elelib_fp: path to file containing ALARA element library (str)
        
    '''
    with open(elelib_fp) as ALARA_Lib:
        libLines = ALARA_Lib.readlines()
    num_lines = len(libLines)
    density_dict = {}
    line_num = 0
    while line_num < num_lines:
        element_data = libLines[line_num].strip().split()
        element_name = element_data[0].lower()
        density_dict[element_name] = float(element_data[3])
        line_num += int(element_data[4]) + 1
    return density_dict

def make_materials(element, density_dict):
    '''
    inputs:
        element: elemental symbol of chosen element (str)
        density_dict: dictionary with key = elemental symbol & value = density [g/cm^3]
        
    outputs:
        mats : OpenMC Materials object
    '''
    mat = openmc.Material(material_id=1, name=element)
    mat.add_element(element, 1.00)
    mat.set_density('g/cm3', density_dict.get(element.lower()))
    mats = openmc.Materials([mat])
    return mats

def make_spherical_shell(material, thickness, inner_radius):  
    '''
    Creates a spherical shell with its own material and inner/outer radius.
    
    inputs:
        material: OpenMC Material object/iterable of OpenMC Material
        thickness: radial thickness (float) of material
        inner_radius : inner radius of material
    
    outputs:
        geometry: OpenMC Geometry object
    '''
    inner_sphere = openmc.Sphere(r = inner_radius)
    cells = [openmc.Cell(fill = None, region = -inner_sphere)]
    outer_radius = inner_radius + thickness
    outer_sphere = openmc.Sphere(r = outer_radius, boundary_type = 'vacuum')
    cells.append(openmc.Cell(fill = material, region = +inner_sphere & -outer_sphere))   
    cells.append(openmc.Cell(fill = None, region = +outer_sphere)) 
    geometry = openmc.Geometry(cells)    
    return geometry

def make_source(energy):
    point_source = openmc.stats.Point(xyz=(0.0, 0.0, 0.0))
    energy_dist = openmc.stats.Discrete(energy, 1.0)
    source = openmc.Source(space = point_source, energy = energy_dist, strength = 1.0, particle = 'neutron')
    return source

def make_tallies(tallied_cells):
    particle_filter = openmc.ParticleFilter('neutron')
    cell_filter = openmc.CellFilter(tallied_cells)
    energy_filter_flux = openmc.EnergyFilter.from_group_structure("VITAMIN-J-175")
    
    neutron_tally = openmc.Tally(name="Neutron tally")
    neutron_tally.scores = ['flux', 'elastic', 'absorption']
    neutron_tally.filters = [cell_filter, particle_filter]

    spectrum_tally = openmc.Tally(name="Neutron flux spectrum")
    spectrum_tally.filters = [energy_filter_flux, cell_filter, particle_filter]
    spectrum_tally.scores = ['flux']
    
    return openmc.Tallies([neutron_tally, spectrum_tally])

def make_settings(source, total_batches, inactive_batches, num_particles, run_mode):
    sets = openmc.Settings()
    sets.batches = total_batches
    sets.inactive = inactive_batches
    sets.particles = num_particles
    sets.source = source
    sets.run_mode = run_mode
    return sets

def deplete_ss(chain_file_path, model, inner_radius, thickness, timesteps, source_rates, norm_mode, timestep_units):
    chain_file = chain_file_path 
    material = model.materials[0]
    material.depletable = True
    material.volume = 4.0/3.0 * np.pi * ((inner_radius+thickness)**3 - inner_radius**3)
    operator = openmc.deplete.CoupledOperator(model, chain_file, normalization_mode = norm_mode)
    integrator = openmc.deplete.PredictorIntegrator(operator, timesteps, source_rates = source_rates, timestep_units = timestep_units)
    return integrator

# Specify inputs and execute all functions:

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--yaml_filepath', default = 'OpenMC_SS_YAML.yaml', help="Path to YAML file containing required inputs (str)")
    parser.add_argument(
        '--run_depletion',
        default = 'True',
        choices=['True', 'False'],
        help='Specify whether to run depletion simulation (true/false)',
        )
    args = parser.parse_args()
    return args
   
def read_yaml(args):
    with open(args.yaml_filepath, 'r') as file:
        inputs = yaml.safe_load(file)
    return inputs

def create_model(inputs):
    geometry = inputs['geom_info']
    settings_info = inputs['settings_info']
    densities = alara_element_densities(inputs['elelib_fp'])
    materials = make_materials(inputs['element'], 
                        densities)
    element = materials[0]
    spherical_shell_geom = make_spherical_shell(element, 
                    geometry['thickness'], 
                    geometry['inner_radius'])
    source_dists = make_source(inputs['particle_energy'])
    sets = make_settings(source_dists,
                    settings_info['total_batches'], 
                    settings_info['inactive_batches'], 
                    settings_info['num_particles'], 
                    settings_info['run_mode'])
    # tallied cells = all cells with non-void material
    tallied_cells = list(spherical_shell_geom.get_all_material_cells().values())
    talls = make_tallies(tallied_cells)
    model = openmc.model.Model(geometry = spherical_shell_geom, materials = materials, settings = sets, tallies = talls)
    return model

def run_depletion(inputs):
    geom_info = inputs['geom_info']
    dep_params = inputs['depletion_params']
    model_file = dep_params['model_file']
    model = openmc.model.Model.from_model_xml(path=model_file)
    
    integrator = deplete_ss(dep_params['chain_file'],
              model, 
              geom_info['inner_radius'], 
              geom_info['thickness'], 
              dep_params['timesteps'], 
              dep_params['source_rates'],
              dep_params['norm_mode'], 
              dep_params['timestep_units'])     
    
    integrator.integrate() 

def main():        
    args = parse_args()
    inputs = read_yaml(args)
    model = create_model(inputs)
    model.export_to_model_xml()
    if args.run_depletion.lower() == 'true':
        run_depletion(inputs)

if __name__ == "__main__":
    main()

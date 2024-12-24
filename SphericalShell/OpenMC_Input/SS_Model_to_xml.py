import openmc
import yaml
import argparse
from OpenMC_SS_Material import alara_element_densities, make_materials
from OpenMC_SS_Geometry import make_spherical_shell
from OpenMC_Source_Tallies_Model import make_source, settings, tallies, create_openmc_model

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--yaml_filepath', required=True, help="Path to YAML file containing required inputs (str)")
    args = parser.parse_args()
    yaml_filepath = args.yaml_filepath
    with open(yaml_filepath, 'r') as file:
        inputs = yaml.safe_load(file)
    geometry = inputs['geom_info']    
    settings_info = inputs['settings_info']
    # for OpenMC_SS_Material
    aed = alara_element_densities(inputs['elelib_fp'])
    materials = make_materials(inputs['element'], 
                           aed, 
                           geometry['inner_radius'], 
                           geometry['thickness'])
    # for OpenMC_SS_Geometry
    mss = make_spherical_shell(materials[0][0], 
                               geometry['thickness'], 
                               geometry['inner_radius'])
    # for OpenMC_Source_Tallies_Model
    source_dists = make_source(inputs['particle_energy'])
    sets = make_settings(ms[0], ms[1], 
                    settings_info['total_batches'], 
                    settings_info['inactive_batches'], 
                    settings_info['num_particles'], 
                    settings_info['run_mode'])
    # tallied cells = all cells with non-void material
    tallied_cells = list(mss.get_all_material_cells().values())
    talls = make_tallies(tallied_cells)
    com = create_openmc_model(materials[0], mss, talls, sets)
    com.export_to_model_xml()

if __name__ == "__main__":
    main()

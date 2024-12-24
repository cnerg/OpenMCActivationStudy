import openmc
import yaml
import argparse
from OpenMC_SS_Material import alara_element_densities, make_materials
from OpenMC_SS_Geometry import make_spherical_shell
from OpenMC_Source_Tallies_Model import make_source, settings, tallies

def main():
    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument('--yaml_filepath', default = 'OpenMC_SS_YAML.yaml', help="Path to YAML file containing required inputs (str)")
        args = parser.parse_args()
        return args
   
    def read_yaml(args):
        with open(args.yaml_filepath, 'r') as file:
            inputs = yaml.safe_load(file)
        return inputs
    
    def create_model(inputs):
        geometry = inputs['geom_info']    
        settings_info = inputs['settings_info']
        # for OpenMC_SS_Material
        densities = alara_element_densities(inputs['elelib_fp'])
        materials = make_materials(inputs['element'], 
                            densities)
        # for OpenMC_SS_Geometry
        element = materials[0]
        spherical_shell_geom = make_spherical_shell(element, 
                                geometry['thickness'], 
                                geometry['inner_radius'])
        # for OpenMC_Source_Tallies_Model
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
        
    arg_parser = parse_args()
    yaml_reader = read_yaml(arg_parser)
    model_creator = create_model(yaml_reader)
    model_creator.export_to_model_xml()

if __name__ == "__main__":
    main()

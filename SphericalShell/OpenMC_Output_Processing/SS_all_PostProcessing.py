import openmc
import openmc.deplete
import yaml
import argparse
from SS_Tally_PostProcessing import extract_tally_values, plot_flux_tally, save_data
from SS_Depletion_PostProcessing import extract_nuclides, plot_save_data
from OpenMC_SS_Material import alara_element_densities, make_element
#from SS_Model_to_xml_gitedit import main

def main() : 
    parser = argparse.ArgumentParser()
    parser.add_argument('--yaml_postprocess_path', default = 'SS_Post_Processing_YAML.yaml', help="Path to YAML file containing inputs for post-processing (str)")
    parser.add_argument('--yaml_model_path', default = 'OpenMC_SS_YAML.yaml', help="Path to YAML file containing required inputs to build model (str)")
    args = parser.parse_args()
    yaml_postprocess_path = args.yaml_postprocess_path
    yaml_model_path = args.yaml_model_path
    
    with open(yaml_postprocess_path, 'r') as yaml_pp :
        pp_inputs = yaml.safe_load(yaml_pp)
    indices = pp_inputs['indices']
    filepaths = pp_inputs['filepaths']
    units = pp_inputs['units']
    
    with open(yaml_model_path, 'r') as model_file:
        model_inputs = yaml.safe_load(model_file)
    elelib_fp = model_inputs['elelib_fp']    
    element = model_inputs['element']
    inner_radius = model_inputs['geom_info']['inner_radius']
    thickness = model_inputs['geom_info']['thickness'] 
    
    # Tally post-processing:
    etv = extract_tally_values(filepaths['statepoint_file_path'])
    pft = plot_flux_tally(etv[0], 
                          indices['flux_tally_id'], 
                          indices['energy_filter_index']) 
    aed = alara_element_densities(elelib_fp)  
    me = make_element(element, aed, inner_radius, thickness)
    mat_vol = me[1]
    sd = save_data(pft[0], mat_vol, etv[1])
    
    # Depletion post-processing:
    en = extract_nuclides(filepaths['dep_file_path'], 
                          units['time_units'], 
                          indices['depletable_mat_index'])    
    psd = plot_save_data(en[0], en[1], en[2], en[3], units['nuc_units'])

if __name__ == "__main__":
    main()

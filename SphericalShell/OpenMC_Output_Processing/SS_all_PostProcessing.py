import openmc
import openmc.deplete
import yaml
import argparse
from SS_Tally_PostProcessing import extract_tally_values, plot_flux_tally, save_data
from SS_Depletion_PostProcessing import extract_nuclides, plot_save_data
from OpenMC_SS_Material import alara_element_densities, make_element

def main() : 
    parser = argparse.ArgumentParser()
    parser.add_argument('--yaml_postprocess_path', required=True, help="Path to YAML file containing inputs for post-processing (str)")
    parser.add_argument('--yaml_model_path', required = True, help = "Path to YAML file containing inputs for openmc model (str)")
    args = parser.parse_args()
    yaml_postprocess_path = args.yaml_postprocess_path
    yaml_model_path = args.yaml_model_path
    
    with open(yaml_postprocess_path, 'r') as yaml_pp :
        pp_inputs = yaml.safe_load(yaml_pp)
    indices = pp_inputs['indices']
    filepaths = pp_inputs['filepaths']
    units = pp_inputs['units']
    # tally post-processing:
    etv = extract_tally_values(filepaths['statepoint_file_path'])
    pft = plot_flux_tally(etv[0], 
                          indices['flux_tally_id'], 
                          indices['energy_filter_index'])
    with open(yaml_model_path, 'r') as yaml_model :
        inputs = yaml.safe_load(yaml_model)
    geom_info = inputs['geom_info']    
    aed = alara_element_densities(inputs['elelib_fp'])   
    me = make_element(inputs['element'], aed, geom_info['inner_radius'], geom_info['thickness']) #needed to find material volume
    mat_vol = me[1]
    sd = save_data(pft[0], mat_vol, etv[1])
    # depletion post-processing:
    en = extract_nuclides(filepaths['dep_file_path'], 
                          units['time_units'], 
                          indices['depletable_mat_index'])    
    psd = plot_save_data(en[0], en[1], en[2], en[3], 
                         units['nuc_units'])

if __name__ == "__main__":
    main()    

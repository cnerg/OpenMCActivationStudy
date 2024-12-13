import openmc
import openmc.deplete
import yaml
from SS_Tally_PostProcessing import extract_tally_values, plot_flux_tally, save_data
from SS_Depletion_PostProcessing import extract_nuclides, plot_save_data
from SS_Model_to_xml import rs 

file = open("SS_Post_Processing_YAML.yaml")
inputs = yaml.safe_load(file)

# SS_Tally_Post_Processing :
statepoint_file_path = inputs['filepaths']['statepoint_file_path']    
flux_tally_id = inputs['indices']['flux_tally_id']
energy_filter_index = inputs['indices']['energy_filter_index']

# SS_Depletion_Post_Processing :
dep_file_path = inputs['filepaths']['dep_file_path']
time_units = inputs['units']['time_units']   
nuc_units = inputs['units']['nuc_units']
depletable_mat_index = inputs['indices']['depletable_mat_index']

def post_process(statepoint_file_path, flux_tally_id, energy_filter_index, dep_file_path, time_units, depletable_mat_index) :
    # Tally Post-Processing :
    etv = extract_tally_values(statepoint_file_path)
    pft = plot_flux_tally(etv[0], flux_tally_id, energy_filter_index)
    mat_vol = rs[1][1]
    sd = save_data(pft[0], mat_vol, etv[1])
    # Depletion Post-Processing :
    en = extract_nuclides(dep_file_path, time_units, depletable_mat_index)
    psd = plot_save_data(en[0], en[1], en[2], en[3], nuc_units)
    
post_processing = post_process(statepoint_file_path, flux_tally_id, energy_filter_index, dep_file_path, time_units, depletable_mat_index)   

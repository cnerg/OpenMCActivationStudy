import openmc
import openmc.deplete
import yaml
from SS_Depletion import deplete_ss
from SS_Model_to_xml import rs

file = open("OpenMC_SS_YAML.yaml")
inputs = yaml.safe_load(file)

chain_file_path = inputs['depletion_params']['chain_file']
model = rs[0]
inner_radius = inputs['geom_info']['inner_radius']
thickness = inputs['geom_info']['thickness']
times_post_boc = inputs['depletion_params']['times_post_boc']
particle_source_rates = inputs['depletion_params']['source_rates']
norm_mode = inputs['depletion_params']['norm_mode']
timestep_units = inputs['depletion_params']['timestep_units']

dss = deplete_ss(chain_file_path, model, inner_radius, thickness, times_post_boc, particle_source_rates, norm_mode, timestep_units)
integration = dss.integrate()
import openmc
import yaml
from OpenMC_SS_Material import alara_element_densities, make_element
from OpenMC_SS_Geometry import make_spherical_shell
from OpenMC_Source_Tallies_Model import make_source, settings, tallies, create_openmc_model

file = open("OpenMC_SS_YAML.yaml")
inputs = yaml.safe_load(file)

# for OpenMC_SS_Material :    
elelib_fp = inputs['elelib_fp']
element = inputs['element']

# for OpenMC_SS_Geometry :
thickness = inputs['geom_info']['thickness']
inner_radius = inputs['geom_info']['inner_radius']

# for OpenMC_Source_Tallies_Model :
energy = inputs['particle_energy']
total_batches = inputs['settings_info']['total_batches']
inactive_batches = inputs['settings_info']['inactive_batches']
num_particles = inputs['settings_info']['num_particles']
run_mode = inputs['settings_info']['run_mode']

def run_ss_scripts(elelib_fp, element, thickness, inner_radius, energy, total_batches, inactive_batches, num_particles, run_mode) :
    #OpenMC_SS_Material :
    aed = alara_element_densities(elelib_fp)
    element = make_element(element, aed, inner_radius, thickness)
    #OpenMC_SS_Geometry :
    # The 1st argument is an OpenMC Material object, derived from a Materials object
    mss = make_spherical_shell(element[0][0], thickness, inner_radius)
    #OpenMC_Source_Tallies_Model :
    ms = make_source(energy)
    sets = settings(ms[0], ms[1], total_batches, inactive_batches, num_particles, run_mode)
    # tallied cells = all cells with non-void material
    tallied_cells = list(mss.get_all_material_cells().values())
    talls = tallies(tallied_cells)
    com = create_openmc_model(element[0], mss, talls, sets)
    com.export_to_model_xml()
    return com, element
    
rs = run_ss_scripts(elelib_fp, element, thickness, inner_radius, energy, total_batches, inactive_batches, num_particles, run_mode)
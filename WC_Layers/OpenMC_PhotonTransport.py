import openmc
import numpy as np
import yaml
import argparse
from Source_Mesh_Reader import extract_source_data
from TwoLayers_Materials import alara_element_densities, make_materials
from TwoLayers_Geometry import make_spherical_shells

def make_source(bounds, cells, mesh_file, source_mesh_index, esd):
    '''
    Creates a list of OpenMC sources, complete with the relevant space and energy distributions
    
    inputs:
        bounds : iterable of photon energy bounds (float)
        cells: list of OpenMC Cell objects
        mesh_file: .h5/.h5m mesh onto which photon source will be distributed
        source_mesh_index: index specifying the photon source from which data is extracted
        
    output:
        source_list: list of OpenMC independent sources
        total_mesh: OpenMC Unstructured Mesh object
    '''

    source_list = []
    total_mesh = openmc.UnstructuredMesh(mesh_file, library='moab')
    for index, (lower_bound, upper_bound) in enumerate(zip(bounds[:-1],bounds[1:])):
        mesh_dist = openmc.stats.MeshSpatial(total_mesh, strengths=esd[source_mesh_index][:,index], volume_normalized=False)
        energy_dist = openmc.stats.Uniform(a=lower_bound, b=upper_bound)
        source_list.append(openmc.IndependentSource(space=mesh_dist, energy=energy_dist, strength=np.sum(esd[source_mesh_index][:, index]), particle='photon', domains=cells))
    return source_list, total_mesh

def make_tallies(total_mesh, tallied_cells, coeff_geom):
    '''
    Creates tallies and assigns energy, spatial, and particle filters.
    
    inputs: 
        total_mesh: OpenMC unstructured mesh object
        tallied_cells: OpenMC Cell/iterable of OpenMC Cell objects/iterable of Cell ID #
        
    outputs:
        OpenMC Tallies object
    '''
    particle_filter = openmc.ParticleFilter('photon')
    total_filter = openmc.MeshFilter(total_mesh)
    
    photon_tally = openmc.Tally(tally_id=1, name="Photon tally")
    photon_tally.scores = ['flux', 'elastic', 'absorption']
    
    cell_filter = openmc.CellFilter(tallied_cells)
    photon_tally.filters = [cell_filter, total_filter, particle_filter]
    
    # Obtain coefficients to calculate effective dose
    dose_energy, dose = openmc.data.dose_coefficients('photon', geometry=coeff_geom)
    dose_filter = openmc.EnergyFunctionFilter(dose_energy, dose)

    # An energy filter is created to assign to the flux tally.
    energy_filter_flux = openmc.EnergyFilter.from_group_structure("VITAMIN-J-42")

    spectrum_tally = openmc.Tally(tally_id=2, name="Flux spectrum")
    # Implementing energy and cell filters for flux spectrum tally
    spectrum_tally.filters = [cell_filter, total_filter, energy_filter_flux, particle_filter, dose_filter]
    spectrum_tally.scores = ['flux']
    
    talls = openmc.Tallies([photon_tally, spectrum_tally])
    return talls

def make_settings(source_list, tot_batches, inactive_batches, num_particles, run_mode):
    '''
    Creates an OpenMC Settings object
    
    inputs:
        source_list: iterable of OpenMC SourceBase objects
    outputs:
        sets: OpenMC Settings object
    '''
    sets = openmc.Settings()
    sets.batches = tot_batches
    sets.inactive = inactive_batches
    sets.particles = num_particles
    sets.source = source_list
    sets.run_mode = run_mode
    return sets

def export_to_xml(geom_object, mats_object, sets_object, talls_object):
    model = openmc.model.Model(geometry = geom_object, materials = mats_object, settings = sets_object, tallies = talls_object)
    return model    
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--Photon_Transport_YAML', default = "PhotonTransport_Inputs.yaml", help="Path (str) to YAML containing inputs for OpenMC_PhotonTransport")
    parser.add_argument('--Mesh_Reader_YAML', default = 'Source_Mesh_Reader_Inputs.yaml', help="Path (str) to YAML containing inputs for Source_Mesh_Reader")
    args = parser.parse_args()
    transport_yaml = args.Photon_Transport_YAML
    smr_yaml = args.Mesh_Reader_YAML
    
    with open(transport_yaml, 'r') as transport_file:
        transport_data = yaml.safe_load(transport_file)
    with open(smr_yaml, 'r') as smr_file:
        smr_data = yaml.safe_load(smr_file)    
        
    # TwoLayers_Materials:    
    alara_fp = transport_data['filename_dict']['alara_el_lib']
    elements = transport_data['mat_info']['element_list']
    
    aed = alara_element_densities(alara_fp)
    mm = make_materials(elements, aed)
    
    # TwoLayers_Geometry:
    materials = mm
    geom_info = transport_data['geom_info']  
    mss = make_spherical_shells(materials, geom_info['thicknesses'], geom_info['inner_radius'], geom_info['outer_boundary_type'])
    
    mesh_file = transport_data['filename_dict']['mesh_file_openmc']
    source_mesh_index = transport_data['file_indices']['source_mesh_index']
    bounds = transport_data['source_info']['phtn_e_bounds']
    cells = list(mss.get_all_cells().values())
    tallied_cells = list(mss.get_all_material_cells().values())
    coeff_geom = transport_data['coeff_geom']
    
    source_mesh_list = smr_data['source_meshes']
    num_elements = smr_data['num_elements']
    photon_groups = smr_data['photon_groups']
    esd = extract_source_data(source_mesh_list, num_elements, photon_groups)
   
    ms = make_source(bounds,cells, mesh_file, source_mesh_index, esd)
    mt = make_tallies(ms[1], tallied_cells, coeff_geom)
    
    s = transport_data['settings_info']
    
    ms = make_settings(ms[0], s['batches'], s['inactive_batches'], s['particles'], s['run_mode'])
    
    etx = export_to_xml(mss, mm, ms, mt)
    etx.export_to_model_xml()
    
if __name__ == "__main__":
    main()

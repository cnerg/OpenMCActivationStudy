import openmc
import numpy as np
import yaml
import argparse
from Source_Mesh_Reader import extract_source_data
from TwoLayers_Materials import alara_element_densities, make_materials
from TwoLayers_Geometry import make_spherical_shells

def make_source(bounds, cells, mesh_file, source_mesh_index, source_data):
    '''
    Creates a list of OpenMC sources, complete with the relevant space and energy distributions
    
    inputs:
        bounds : iterable of photon energy bounds (float)
        cells: list of OpenMC Cell objects
        mesh_file: .h5/.h5m mesh onto which photon source will be distributed
        source_mesh_index: index specifying the photon source from which data is extracted
        
    output:
        source_list: list of OpenMC independent sources
        unstructured_mesh: OpenMC Unstructured Mesh object
    '''

    source_list = []
    unstructured_mesh = openmc.UnstructuredMesh(mesh_file, library='moab')
    for index, (lower_bound, upper_bound) in enumerate(zip(bounds[:-1],bounds[1:])):
        mesh_dist = openmc.stats.MeshSpatial(unstructured_mesh, strengths=source_data[source_mesh_index][:,index], volume_normalized=False)
        energy_dist = openmc.stats.Uniform(a=lower_bound, b=upper_bound)
        source_list.append(openmc.IndependentSource(space=mesh_dist, energy=energy_dist, strength=np.sum(source_data[source_mesh_index][:, index]), particle='photon', domains=cells))
    return source_list, unstructured_mesh

def make_tallies(unstructured_mesh, tallied_cells, coeff_geom, bounds):
    '''
    Creates tallies and assigns energy, spatial, and particle filters.
    
    inputs: 
        unstructured_mesh: OpenMC unstructured mesh object
        tallied_cells: OpenMC Cell/iterable of OpenMC Cell objects/iterable of Cell ID #
        coeff_geom : irradation geometry associated with dose coefficient calculation
        bounds : energy bounds associated with photon source from ALARA
    outputs:
        OpenMC Tallies object
    '''
    particle_filter = openmc.ParticleFilter('photon')
    mesh_filter = openmc.MeshFilter(unstructured_mesh)
    
    photon_tally = openmc.Tally(tally_id=1, name="Photon tally")
    photon_tally.scores = ['flux', 'absorption']
    energy_filter_photon_tally = openmc.EnergyFilter(bounds)  
    
    cell_filter = openmc.CellFilter(tallied_cells)
    photon_tally.filters = [energy_filter_photon_tally, particle_filter]
    
    # Obtain coefficients to calculate effective dose
    dose_energy, dose = openmc.data.dose_coefficients('photon', geometry=coeff_geom)
    dose_filter = openmc.EnergyFunctionFilter(dose_energy, dose)

    # An energy filter is created to assign to the flux tally.
    energy_filter_flux = openmc.EnergyFilter.from_group_structure("VITAMIN-J-42")

    spectrum_tally = openmc.Tally(tally_id=2, name="Flux spectrum")
    # Implementing energy and cell filters for flux spectrum tally
    spectrum_tally.filters = [cell_filter, mesh_filter, energy_filter_flux, particle_filter, dose_filter]
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
    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument('--Photon_Transport_YAML', default = "PhotonTransport_Inputs.yaml", help="Path (str) to YAML containing inputs for OpenMC_PhotonTransport")
        parser.add_argument('--Mesh_Reader_YAML', default = 'Source_Mesh_Reader_Inputs.yaml', help="Path (str) to YAML containing inputs for Source_Mesh_Reader")
        args = parser.parse_args()
        return args

    def read_yamls(args):
        with open(args.Photon_Transport_YAML, 'r') as transport_file:
            transport_data = yaml.safe_load(transport_file)
        with open(args.Mesh_Reader_YAML, 'r') as smr_file:
            mesh_data = yaml.safe_load(smr_file)    
        return transport_data, mesh_data

    def model_xml(transport_data, mesh_data):
    
        # TwoLayers_Materials:   
        filenames = transport_data['filename_dict']
        alara_fp = filenames['alara_el_lib']    
        density_dict = alara_element_densities(filenames['alara_el_lib'])
        materials = make_materials(transport_data['mat_info']['element_list'], 
                               density_dict)
    
        # TwoLayers_Geometry:
        geom_info = transport_data['geom_info']  
        coeff_geom = transport_data['coeff_geom']
        settings_info = transport_data['settings_info']
        source_mesh_list = mesh_data['source_meshes']
        num_elements = mesh_data['num_elements']
        photon_groups = mesh_data['photon_groups']  

        layers = zip(materials, geom_info['thicknesses'])
        
        spherical_shell_geom = make_spherical_shells(geom_info['inner_radius'], 
                                                     layers,
                                                     geom_info['outer_boundary_type'])
        cells = list(spherical_shell_geom.get_all_cells().values())
        tallied_cells = list(spherical_shell_geom.get_all_material_cells().values())
        
        source_data = extract_source_data(source_mesh_list, num_elements, photon_groups)
        source_list, unstructured_mesh = make_source(transport_data['source_info']['phtn_e_bounds'],
                     cells, transport_data['filename_dict']['mesh_file_openmc'], 
                     transport_data['file_indices']['source_mesh_index'], 
                     source_data)
        tallies = make_tallies(unstructured_mesh, tallied_cells, coeff_geom, transport_data['source_info']['phtn_e_bounds'])
        settings = make_settings(source_list, settings_info['batches'], settings_info['inactive_batches'], settings_info['particles'], settings_info['run_mode'])
    
        model_obj = export_to_xml(spherical_shell_geom, materials, settings, tallies)
        model_obj.export_to_model_xml()
    
    args = parse_args()
    transport_data, mesh_data = read_yamls(args)
    model_xml(transport_data, mesh_data)
    
if __name__ == "__main__":
    main()

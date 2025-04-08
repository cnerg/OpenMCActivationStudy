import openmc
import numpy as np
import yaml
import argparse
import h5py

#Set up materials for model:

def alara_element_densities(alara_fp):
    '''
    Creates a dictionary where keys = element names (str) and values = element density (float)
    inputs:
        alara_filepath : path to ALARA element library
    '''
    with open(alara_fp) as ALARA_Lib:
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
    
def make_materials(elements, density_dict):
    '''
    Creates an OpenMC Materials object using user-specified elements
    inputs:
        elements: iterable of element names (str)
        density_dict: dictionary with keys = element names (str) and values = element density (float)
    '''    
    mats = openmc.Materials([])
    for element_id, element in enumerate(elements):
        mat = openmc.Material(material_id=element_id+1, name=element)
        mat.add_element(element, 1.00)
        mat.set_density('g/cm3', density_dict.get(element.lower()))
        mats.append(mat)
    return mats

#Set up geometry for model:

def make_spherical_shells(inner_radius, layers, outer_boundary_type):    
    '''
    Creates a set of concentric spherical shells, each with its own material & inner/outer radius.
    inputs:
        inner_radius: the radius of the innermost spherical shell
        layers: iterable of tuples of OpenMC Material object and its respective thickness (float)
    '''
    inner_sphere = openmc.Sphere(r = inner_radius)
    cells = [openmc.Cell(fill = None, region = -inner_sphere)]
    for (material, thickness) in layers:
        outer_radius = inner_radius + thickness
        outer_sphere = openmc.Sphere(r = outer_radius)
        cells.append(openmc.Cell(fill = material, region = +inner_sphere & -outer_sphere))
        inner_radius = outer_radius
        inner_sphere = outer_sphere
    outer_sphere.boundary_type = outer_boundary_type     
    cells.append(openmc.Cell(fill = None, region = +outer_sphere)) 
    geometry = openmc.Geometry(cells)    
    return geometry

def make_neutron_source(energy):

    point_source = openmc.stats.Point(xyz=(0.0, 0.0, 0.0))
    energy_dist = openmc.stats.Discrete(energy, 1.0)
    source = openmc.Source(space = point_source, energy = energy_dist, strength = 1.0, particle = 'neutron')
    return source

def make_neutron_tallies(mesh_file):
    mesh_filter = openmc.MeshFilter(openmc.UnstructuredMesh(mesh_file, library='moab'))
    particle_filter = openmc.ParticleFilter('neutron')
    energy_filter_flux = openmc.EnergyFilter.from_group_structure("VITAMIN-J-175")
    
    neutron_tally = openmc.Tally(name="Neutron tally")
    neutron_tally.scores = ['flux', 'absorption']
    neutron_tally.filters = [particle_filter]

    spectrum_tally = openmc.Tally(name="Neutron flux spectrum")
    spectrum_tally.filters = [mesh_filter, energy_filter_flux, particle_filter]
    spectrum_tally.scores = ['flux']
    
    return openmc.Tallies([neutron_tally, spectrum_tally])

def extract_source_data(source_mesh_list, num_elements, num_photon_groups):
    '''
    Identifies the location of the source density dataset within each mesh file.
    
    input: 
        source_mesh_list: iterable of .h5m filenames (str), whose files contain photon source information
    output: 
        numpy array of source density data with rows = # of mesh elements and columns = # number of photon groups, with one array per source mesh
    ''' 
    sd_list = np.ndarray((len(source_mesh_list), num_elements, num_photon_groups))
    for source_index, source_name in enumerate(source_mesh_list):
         file = h5py.File(source_name, 'r')
         sd_list[source_index,:] = file['tstt']['elements']['Tet4']['tags']['source_density'][:]
    return sd_list   
    
def make_photon_sources(bounds, cells, mesh_file, source_mesh_index, sd_list):
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
        mesh_dist = openmc.stats.MeshSpatial(unstructured_mesh, strengths=sd_list[source_mesh_index][:,index], volume_normalized=False)
        energy_dist = openmc.stats.Uniform(a=lower_bound, b=upper_bound)
        source_list.append(openmc.IndependentSource(space=mesh_dist, energy=energy_dist, strength=np.sum(sd_list[source_mesh_index][:, index]), particle='photon', domains=cells))
    return source_list, unstructured_mesh

def make_photon_tallies(unstructured_mesh, tallied_cells, coeff_geom, bounds):
    '''
    Creates tallies and assigns energy, spatial, and particle filters.
    
    inputs: 
        unstructured_mesh: OpenMC unstructured mesh object
        tallied_cells: OpenMC Cell/iterable of OpenMC Cell objects/iterable of Cell ID #
        coeff_geom : irradation geometry associated with dose coefficient calculation
        bounds : energy bounds associated with photon source from ALARA
        
    outputs:
        talls: OpenMC Tallies object
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
    spectrum_tally.filters = [cell_filter, mesh_filter, energy_filter_flux, particle_filter, dose_filter] ##
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

#--------------
#Execute all functions:

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--yaml_filepath', default = 'OpenMC_ALARA_WC.yaml', help="Path to YAML file containing required inputs for OpenMC-ALARA R2S workflow (str)")
    parser.add_argument("--neutron_transport", default=True, help="Create neutron transport model")
    parser.add_argument("--photon_transport", default=True, help="Create photon transport model")
    args = parser.parse_args()
    return args
   
def read_yaml(args):
    with open(args.yaml_filepath, 'r') as file:
        inputs = yaml.safe_load(file)
    return inputs

def create_materials_obj(inputs):
    densities = alara_element_densities(inputs['filename_dict']['elelib_fp'])
    materials = make_materials(inputs['mat_info']['element_list'], 
                    densities)
    return materials

def create_geometry_obj(materials, inputs):
    geom_info = inputs['geom_info']
    layers = zip(materials, geom_info['thicknesses'])
    geometry = make_spherical_shells(geom_info['inner_radius'], 
                    layers,
                    geom_info['outer_boundary_type'])
    return geometry

#Build neutron transport model:

def create_neutron_model(inputs, materials, geometry):
    settings_info = inputs['settings_info']
    
    source = make_source(inputs['particle_energy'])
    neutron_settings = make_settings(source,
                    settings_info['total_batches'], 
                    settings_info['inactive_batches'], 
                    settings_info['num_particles'], 
                    settings_info['run_mode'])
    neutron_tallies = make_neutron_tallies(inputs['filename_dict']['mesh_file'])
    neutron_model = openmc.model.Model(geometry = geometry, materials = materials, settings = neutron_settings, tallies = neutron_tallies)
    return neutron_model

#Convert the output of R2S Step2 to a format suitable for OpenMC photon transport:

def read_source_mesh(inputs):
    #Find the size of the first source density dataset (assumed to be the same for all other datasets):
    sd_data = h5py.File(inputs['source_meshes'][0], 'r')['tstt']['elements']['Tet4']['tags']['source_density'][:]
    sd_list = extract_photon_source_data(inputs['source_meshes'],

                                      sd_data.shape[0],
                                      sd_data.shape[1])
    return sd_list

#Build photon transport model:

def create_photon_model(inputs, materials, geometry, sd_list):
    settings_info = inputs['settings_info']                                            
    cells = list(geometry.get_all_cells().values())
    tallied_cells = list(geometry.get_all_material_cells().values())
    source_list, unstructured_mesh = make_photon_sources(inputs['source_info']['phtn_e_bounds'],
                cells, 
                inputs['filename_dict']['mesh_file'], 
                inputs['file_indices']['source_mesh_index'], 
                sd_list)
    photon_settings = make_settings(source_list, 
                settings_info['total_batches'], 
                settings_info['inactive_batches'], 
                settings_info['num_particles'], 
                settings_info['run_mode'])
    photon_tallies = make_photon_tallies(unstructured_mesh, tallied_cells, inputs['coeff_geom'], inputs['source_info']['phtn_e_bounds'])
    photon_model = openmc.model.Model(geometry = geometry, materials = materials, settings = photon_settings, tallies = photon_tallies) 
    return photon_model                             

def main():
    args = parse_args()
    inputs = read_yaml(args)
    materials = create_materials_obj(inputs)
    geometry = create_geometry_obj(materials, inputs)

    if args.neutron_transport == True:
        neutron_model = create_neutron_model(inputs, materials, geometry)
        neutron_model.export_to_model_xml(path="neutron_model.xml")
    
    if args.photon_transport == True:
        sd_list = read_source_mesh(inputs)
        photon_model = create_photon_model(inputs, materials, geometry, sd_list)
        photon_model.export_to_model_xml(path="photon_model.xml")

if __name__ == "__main__":
    main()     

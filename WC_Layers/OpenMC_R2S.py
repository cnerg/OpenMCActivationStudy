import openmc
import openmc.deplete
from pathlib import Path
import argparse
import yaml
import numpy as np
import h5py

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
        material.volume = 4.0/3.0 * np.pi * ((outer_radius)**3 - (inner_radius)**3)
        inner_radius = outer_radius
        inner_sphere = outer_sphere
    outer_sphere.boundary_type = outer_boundary_type    
    geometry = openmc.Geometry(cells)    
    return geometry

def make_neutron_source(energy):
    point_source = openmc.stats.Point(xyz=(0.0, 0.0, 0.0))
    energy_dist = openmc.stats.Discrete(energy, 1.0)
    neutron_source = openmc.Source(space = point_source, energy = energy_dist, strength = 1.0)
    return neutron_source

def make_neutron_tallies(mesh_file):
    mesh_filter = openmc.MeshFilter(openmc.UnstructuredMesh(mesh_file, library='moab'))
    particle_filter = openmc.ParticleFilter('neutron')
    energy_filter_flux = openmc.EnergyFilter.from_group_structure("VITAMIN-J-175")

    neutron_flux_tally = openmc.Tally(name="Neutron flux spectrum")
    neutron_flux_tally.filters = [mesh_filter, energy_filter_flux, particle_filter]
    neutron_flux_tally.scores = ['flux']
    
    return openmc.Tallies([neutron_flux_tally])

def make_settings(source, total_batches, inactive_batches, num_particles, run_mode):
    sets = openmc.Settings()
    sets.batches = total_batches
    sets.inactive = inactive_batches
    sets.particles = num_particles
    sets.source = source
    sets.run_mode = run_mode
    return sets

#Only executed if external geometry and materials are imported
def make_depletion_volumes(neutron_model, mesh_file):
    materials = neutron_model.materials
    mesh_file = Path(mesh_file).resolve()   
    unstructured_mesh = openmc.UnstructuredMesh(mesh_file, library='moab') 
    mat_vols = unstructured_mesh.material_volumes(neutron_model, n_samples=25000000)
    #returns dict-like object that maps mat ids to array of volumes equal to # of mesh elements

    total_volumes = {}
    for mat_id, volumes in mat_vols.items():
        total_volumes[mat_id] = np.sum(volumes)     
    for material in materials:
        material.volume = total_volumes[material.id]
    return neutron_model 

def deplete_model(neutron_model, mesh_file, chain_file, timesteps, source_rates, norm_mode, timestep_units):
    materials = neutron_model.materials
    mesh_file = Path(mesh_file).resolve()   
    unstructured_mesh = openmc.UnstructuredMesh(mesh_file, library='moab') 
    model_nuclide_names = []
    for material in materials:
       material.depletable = True
       material_nuclides = material.nuclides
       for material_nuclide in material_nuclides:
            model_nuclide_names.append(material_nuclide.name)

    activation_mats = unstructured_mesh.get_homogenized_materials(neutron_model, n_samples=7000000)
    activation_mats_object = openmc.Materials(activation_mats)
    activation_mats_object.export_to_xml("Activation_Materials.xml")   

    depletion_chain = openmc.deplete.Chain.from_xml(Path(chain_file).resolve())

    nuclide_names = []
    for nuclide in depletion_chain.nuclides:
        if nuclide.name in model_nuclide_names:
            nuclide_names.append(nuclide.name)

    fluxes, micros = openmc.deplete.get_microxs_and_flux(neutron_model, unstructured_mesh, nuclide_names,
                                                          depletion_chain.reactions,
                                                          openmc.mgxs.GROUP_STRUCTURES["VITAMIN-J-175"])
    operator = openmc.deplete.IndependentOperator(openmc.Materials(activation_mats), fluxes, micros, chain_file=chain_file, normalization_mode = norm_mode)    
    copy_timesteps = list(timesteps)
    copy_source_rates = list(source_rates)
    integrator_list = []
    for timestep in timesteps:
        while copy_source_rates[-1] == 0:
            integrator = openmc.deplete.PredictorIntegrator(operator, copy_timesteps, source_rates = copy_source_rates, timestep_units = timestep_units)
            integrator_list.append(integrator)
            copy_timesteps.pop()
            copy_source_rates.pop()

    #Letting integrator_list go from the fewest number of decay intervals to the most        
    integrator_list.reverse() 
    for int_index, integrator in enumerate(integrator_list):
        integrator.integrate(path=f"depletion_results_decay_set_{int_index}.h5")        
    return activation_mats, unstructured_mesh, neutron_model

def extract_source_data(inputs):
    '''
    Identifies the location of the source density dataset within each mesh file.
    output: 
        numpy array of source density data with rows = # of mesh elements and columns = # number of photon groups, with one array per source mesh
    '''
    source_mesh_list = inputs['source_meshes']
    sd_data = h5py.File(source_mesh_list[0], 'r')['tstt']['elements']['Tet4']['tags']['source_density'][:]
    num_elements = sd_data.shape[0]
    num_photon_groups = sd_data.shape[1]
    sd_list = np.ndarray((len(source_mesh_list), num_elements, num_photon_groups))
    for source_index, source_name in enumerate(source_mesh_list):
         file = h5py.File(source_name, 'r')
         sd_list[source_index,:] = file['tstt']['elements']['Tet4']['tags']['source_density'][:]
    return sd_list  
    
def make_alara_photon_sources(bounds, cells, mesh_file, source_mesh_indices, sd_list):
    '''
    Creates a list of OpenMC sources, complete with the relevant space and energy distributions
    
    inputs:
        bounds : iterable of photon energy bounds (float)
        cells: list of OpenMC Cell objects
        mesh_file: .h5/.h5m mesh onto which photon source will be distributed
        source_mesh_indices: iterable of indices specifying the photon source from which data is extracted
        
    output:
        all_sources: dictionary of OpenMC independent sources
    '''
    all_sources = {}
    unstructured_mesh = openmc.UnstructuredMesh(mesh_file, library='moab')
    for source_mesh_index in source_mesh_indices:
        source_list = []
        for index, (lower_bound, upper_bound) in enumerate(zip(bounds[:-1],bounds[1:])):
            mesh_dist = openmc.stats.MeshSpatial(unstructured_mesh, strengths=sd_list[source_mesh_index][:,index], volume_normalized=False)
            energy_dist = openmc.stats.Uniform(a=lower_bound, b=upper_bound)
            source_list.append(openmc.IndependentSource(space=mesh_dist, energy=energy_dist, strength=np.sum(sd_list[source_mesh_index][:, index]), particle='photon', domains=cells))
            all_sources[source_mesh_index] = source_list
    return all_sources

def make_openmc_photon_sources(num_cooling_steps, activation_mats, unstructured_mesh, neutron_model, inputs):
    sd_data = h5py.File(inputs['source_meshes'][0], 'r')['tstt']['elements']['Tet4']['tags']['source_density'][:]
    for int_index in range(num_cooling_steps):
        results = openmc.deplete.Results(f"depletion_results_decay_set_{int_index}.h5")
        photon_sources = np.empty(sd_data.shape[0], dtype=object)      
        activated_mats_list = set(results[-1].index_mat.keys())

        for mat_index, mat in enumerate(activation_mats) :   
            if str(mat.id) in activated_mats_list :
                mat = results[-1].get_material(str(mat.id))
                energy = mat.get_decay_photon_energy()
                if energy == None:
                    photon_source = openmc.IndependentSource(
                        energy = energy,
                        particle = 'photon',
                        strength = 0.0)
                else:    
                    photon_source = openmc.IndependentSource(
                        energy = energy,
                        particle = 'photon',
                        strength = energy.integral())
                
            else:
                photon_source = openmc.IndependentSource(
                    energy = None,
                    particle = 'photon',
                    strength = 0.0)
            
            photon_sources[mat_index] = photon_source
        photon_model = neutron_model
        photon_model.settings.source = openmc.MeshSource(unstructured_mesh, photon_sources)
        photon_model.settings.export_to_xml(f'settings_{int_index}.xml')

    return photon_model
    
def make_photon_tallies(coeff_geom, photon_model, num_cooling_steps):
    dose_energy, dose = openmc.data.dose_coefficients('photon', geometry=coeff_geom)
    dose_filter = openmc.EnergyFunctionFilter(dose_energy, dose)

    #Could also use an unstructured mesh filter here - spherical mesh filter for visualization purposes
    spherical_mesh = openmc.SphericalMesh(np.arange(0, 3800, 10), origin = (0.0, 0.0, 0.0), mesh_id=2, name="spherical_mesh")
    spherical_mesh_filter = openmc.MeshFilter(spherical_mesh)
    
    particle_filter = openmc.ParticleFilter('photon')
    energy_filter_flux = openmc.EnergyFilter.from_group_structure("VITAMIN-J-42")

    flux_tally = openmc.Tally(tally_id=1, name="photon flux tally")
    flux_tally.filters = [spherical_mesh_filter, energy_filter_flux, particle_filter]
    flux_tally.scores = ['flux']

    dose_tally = openmc.Tally(tally_id=2, name="dose tally")
    dose_tally.filters = [spherical_mesh_filter, energy_filter_flux, particle_filter, dose_filter]
    dose_tally.scores = ['flux']

    photon_model.tallies = [flux_tally, dose_tally]
    #Change boundary condition:
    list(photon_model.geometry.get_all_surfaces().keys())[-1].boundary_type="white"    

    for int_index in range(num_cooling_steps):
        #Reassign settings (which contains source) for each decay step
        photon_model.settings = openmc.Settings.from_xml(f'settings_{int_index}.xml')
        photon_model.export_to_model_xml(path=f'photon_model_{int_index}.xml')
        sp_path = photon_model.run(f'photon_model_{int_index}.xml')
        sp_path.rename(f'statepoint_openmc_{int_index}.h5')

#----------------------------------------------------------------------------------
#Define variables and execute all functions:

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--OpenMC_YAML', default = "R2S.yaml", help="Path (str) to YAML containing inputs")
    parser.add_argument('--ext_mat_geom', default = True, help="Specify whether materials and geometry come from external model")
    parser.add_argument('--pyne_r2s', default = False, help="Choose to run pyne r2s steps")
    parser.add_argument('--openmc_r2s', default = False, help="Choose to run openmc r2s steps")
    if args.pyne_r2s and args.openmc_r2s:
        parser.error("Cannot run PyNE and OpenMC workflows at the same time")
    parser.add_argument("--pyne_neutron_transport", default=True, help="If True, run only neutron transport on PyNE model. If False, run only photon transport on PyNE model")
    args = parser.parse_args()
    return args

def read_yaml(yaml_arg):
    '''
    input:
        yaml_arg : output of parse_args() corresponding to args.OpenMC_YAML
    '''
    with open(yaml_arg, 'r') as transport_file:
        inputs = yaml.safe_load(transport_file)
    return inputs
    
def create_neutron_model(inputs, materials, geometry):
    neutron_settings_info = inputs['neutron_settings_info']
    neutron_source = make_neutron_source(inputs['particle_energy'])
    neutron_settings = make_settings(neutron_source,
                    neutron_settings_info['total_batches'],
                    neutron_settings_info['inactive_batches'],
                    neutron_settings_info['num_particles'],
                    neutron_settings_info['run_mode'])   
    neutron_model = openmc.model.Model(geometry = geometry, materials = materials, settings = neutron_settings)
    return neutron_model

#Convert the output of R2S Step2 to a format suitable for OpenMC photon transport:

def create_alara_photon_model(inputs, neutron_model, sd_list):
    '''
    photon_model: neutron_model to which this function adds photon-relevant settings
    '''
    cells = list(photon_model.geometry.get_all_cells().values())
    all_sources = make_alara_photon_sources(
        inputs['source_info']['phtn_e_bounds'],
        cells,
        inputs['filename_dict']['mesh_file'],
        inputs['file_indices']['source_mesh_indices'],
        sd_list,         
    )
    photon_settings_info = inputs['photon_settings_info']
    for source_mesh_index, photon_source in enumerate(list(all_sources.values())):
        photon_settings = make_settings(photon_source,
                        photon_settings_info['total_batches'],
                        photon_settings_info['inactive_batches'],
                        photon_settings_info['num_particles'],
                        photon_settings_info['run_mode'])
        photon_model.settings = photon_settings
        photon_model.settings.export_to_xml(f'settings_{source_mesh_index}.xml')
    return photon_model
    
def import_ext_model(inputs):
    '''
    Import openmc materials and geometry objects from some external openmc model
    '''
    ext_model = openmc.model.Model.from_model_xml(inputs['filename_dict']['ext_model'])
    materials = ext_model.materials
    geometry = ext_model.geometry
    
    #Settings also assigned here
    neutron_model = create_neutron_model(inputs, materials, geometry)

    return neutron_model

def make_native_model(inputs):
    '''
    Run make_materials() and make_spherical_shells()
    '''
    densities = alara_element_densities(inputs['filename_dict']['elelib_fp'])
    materials = make_materials(inputs['mat_info']['element_list'],
        densities)  
    geom_info = inputs['geom_info']
    layers = zip(materials, geom_info['thicknesses'])
    geometry = make_spherical_shells(geom_info['inner_radius'],
        layers,
        geom_info['outer_boundary_type'])
    
    #Settings also assigned here
    neutron_model = create_neutron_model(inputs, materials, geometry)

    return neutron_model

def run_pyne_r2s(inputs, args):
    if args.ext_geom == True: #Import materials and geometry from external model
        neutron_model = import_ext_model(inputs) 
    if args.ext_geom == False: #Run make_materials() and make_spherical_shells() 
        neutron_model = make_native_model(inputs)

    if args.pyne_neutron_transport == True:
        neutron_model.tallies = make_neutron_tallies(inputs['filename_dict']['mesh_file'])
        neutron_model.export_to_model_xml('neutron_model_pyne.xml')
        neutron_model_sp = neutron_model.run('neutron_model_pyne.xml')
        neutron_model_sp.rename('neutron_model_pyne.statepoint.h5')

    if args.pyne_neutron_transport == False:
        sd_list = extract_source_data(inputs)
        photon_model = create_alara_photon_model(inputs, neutron_model, sd_list)
        num_cooling_steps = len(inputs['file_indices']['source_mesh_indices'])
        photon_tallies = make_photon_tallies(inputs['coeff_geom'], photon_model, num_cooling_steps)            

def run_openmc_r2s(inputs, args): 
    dep_params = inputs['dep_params'] 
    if args.ext_geom == True: #Import materials and geometry from external model
        neutron_model = import_ext_model(inputs) 
        neutron_model = make_depletion_volumes(neutron_model, inputs['filename_dict']['mesh_file'])
    if args.ext_geom == False: #Run make_materials() and make_spherical_shells() 
        neutron_model = make_native_model(inputs)

    activation_mats, unstructured_mesh, neutron_model = deplete_model(neutron_model,
        inputs['filename_dict']['mesh_file'],
        dep_params['chain_file'],
        dep_params['timesteps'],
        dep_params['source_rates'],
        dep_params['norm_mode'],
        dep_params['timestep_units'])   
    
    num_cooling_steps = (dep_params['source_rates']).count(0)
    photon_model = make_openmc_photon_sources(num_cooling_steps, activation_mats, unstructured_mesh, neutron_model, inputs)
    photon_tallies = make_photon_tallies(inputs['coeff_geom'], photon_model, num_cooling_steps)

def main():        
    args = parse_args()
    inputs = read_yaml(args.OpenMC_YAML)

    openmc.config['chain_file'] = inputs['dep_params']['chain_file']

    if pyne_r2s == True:
        run_pyne_r2s(inputs, args)
    if openmc_r2s == True:
        run_openmc_r2s(inputs, args)
   
if __name__ == "__main__":
    main()

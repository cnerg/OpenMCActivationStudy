import openmc
import numpy as np
from Source_Mesh_Reader import extract_source_data, Files

# Photon source densities from ALARA output 
esd = extract_source_data(Files)

def make_source(cells, mesh_file):
      '''
    Creates a list of OpenMC sources, complete with the relevant space and energy distributions
    
    inputs:
        cells: list of OpenMC Cell objects
        mesh_file: .h5/.h5m mesh onto which photon source will be distributed
        
    output:
        source_list: list of OpenMC independent sources
        total_mesh: OpenMC Unstructured Mesh object
    '''
    source_list = []
    total_mesh = openmc.UnstructuredMesh(mesh_file, library='moab')
    for index, (lower_bound, upper_bound) in enumerate(zip(bounds[:-1],bounds[1:])):
        mesh_dist = openmc.stats.MeshSpatial(total_mesh, strengths=esd[mesh_index][:,index], volume_normalized=False)
        energy_dist = openmc.stats.Uniform(a=lower_bound, b=upper_bound)
        source_list.append(openmc.IndependentSource(space=mesh_dist, energy=energy_dist, strength=np.sum(esd[mesh_index][:, index]), particle='photon', domains=cells))
    return source_list, total_mesh

def tallies(total_mesh, tallied_cells):
      '''
    Creates tallies and assigns energy, spatial, and particle filters.
    
    inputs: 
        total_mesh: OpenMC unstructured mesh object
        tallied_cells: OpenMC Cell/iterable of OpenMC Cell objects/iterable of Cell ID #
        
    outputs:
        talls: OpenMC Tallies object
    '''
    particle_filter = openmc.ParticleFilter('photon')
    total_filter = openmc.MeshFilter(total_mesh)
    
    neutron_tally = openmc.Tally(tally_id=1, name="Neutron tally")
    neutron_tally.scores = ['flux', 'elastic', 'absorption']
    
    # Implementing filter for neutron tally through shells with material
    cell_filter = openmc.CellFilter(tallied_cells)
    neutron_tally.filters = [cell_filter, total_filter]

    # Vitamin-J energy filter created to assign to the flux tally.
    energy_filter_flux = openmc.EnergyFilter.from_group_structure("VITAMIN-J-42")

    spectrum_tally = openmc.Tally(tally_id=2, name="Flux spectrum")
    spectrum_tally.filters = [cell_filter, total_filter, energy_filter_flux, particle_filter]
    spectrum_tally.scores = ['flux']
    
    talls = openmc.Tallies([neutron_tally, spectrum_tally])
    return talls
  
def settings(source_list, total_batches, inactive_batches, num_particles, run_mode):
      '''
    Creates an OpenMC Settings object
    
    inputs:
        source_list: iterable of OpenMC SourceBase objects
    outputs:
        sets: OpenMC Settings object
    '''
    sets = openmc.Settings()
    sets.batches = total_batches
    sets.inactive = inactive_batches
    sets.particles = num_particles
    sets.source = source_list
    sets.run_mode = run_mode
    return sets

def create_openmc_model(geom_object, mats_object, sets_object, talls_object):
    model = openmc.model.Model(geometry = geom_object, materials = mats_object, settings = sets_object, tallies = talls_object)
    return model

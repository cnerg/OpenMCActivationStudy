# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 10:07:33 2024

@author: Anupama Rajendra
"""

import openmc
import argparse
import numpy as np
from Source_Mesh_Reader import extract_source_data, Files
from TwoLayers_Materials import *
from TwoLayers_Geometry import *

# These will be moved to a yaml file
openmc_mesh_file = 'OpenMC_Mesh.h5m'
mesh_index = 0
esd = extract_source_data(Files)

parser = argparse.ArgumentParser(description="Specify required inputs: file path to ALARA Element Library, element name, inner radius [cm], outer_radius [cm]")

#Required user inputs:
parser.add_argument('--filepath', type=str, required=True)    
parser.add_argument('--element_1', type=str, required=True)
parser.add_argument('--element_2', type=str, required=True)
#Shell radii are left as user inputs, but the mesh is specific to W_inner_radius = 1000, W_outer_radius = 1005, C_inner_radius = 995
parser.add_argument('--W_inner_radius', type=float, required=True)
parser.add_argument('--W_outer_radius', type=float, required=True)
parser.add_argument('--C_inner_radius', type=float, required=True)

args = parser.parse_args()

fp = args.filepath
E_1 = args.element_1
E_2 = args.element_2
R_W_1 = args.W_inner_radius
R_W_2 = args.W_outer_radius
R_C_1 = args.C_inner_radius

# fp = '../../ALARA/data/elelib.std'
# E_1 = 'W'
# E_2 = 'C'
# R_W_1 = 1000
# R_W_2 = 1005
# R_C_1 = 995
bounds = np.array([0, 
                 1.00e+4, 
                 2.00e+4, 
                 5.00e+4, 
                 1.00e+5, 
                 2.00e+5, 
                 3.00e+5, 
                 4.00e+5, 
                 6.00e+5, 
                 8.00e+5, 
                 1.00e+6, 
                 1.22e+6, 
                 1.44e+6, 
                 1.66e+6, 
                 2.00e+6, 
                 2.50e+6, 
                 3.00e+6, 
                 4.00e+6, 
                 5.00e+6, 
                 6.50e+6, 
                 8.00e+6, 
                 1.00e+7, 
                 1.20e+7, 
                 1.40e+7, 
                 2.00e+7])

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
    for index, (lower_bound, upper_bound) in enumerate(zip(bounds[:-1],bounds[1:[)):
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
  
def settings(source_list):
      '''
    Creates an OpenMC Settings object
    
    inputs:
        source_list: iterable of OpenMC SourceBase objects
    outputs:
        sets: OpenMC Settings object
    '''
    sets = openmc.Settings()
    sets.batches = 10
    sets.inactive = 1
    sets.particles = 100000
    sets.source = source_list
    sets.run_mode = 'fixed source'
    return sets

def create_openmc_model(geom_object, mats_object, sets_object, talls_object):
    model = openmc.model.Model(geometry = geom_object, materials = mats_object, settings = sets_object, tallies = talls_object)
    return model

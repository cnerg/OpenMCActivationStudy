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
Mesh_File = 'OpenMC_Mesh.h5m'
mesh_index = 0
esd = extract_source_data(Files)
inner_radius = 995

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
bounds = [0, 
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
                 2.00e+7]

#Define source:
def make_source(cells):
    source_list = []
    total_mesh = openmc.UnstructuredMesh(Mesh_File, library='moab')
    for index, bound in enumerate(bounds[:-1]):
        mesh_dist = openmc.stats.MeshSpatial(total_mesh, strengths=esd[mesh_index][:,index], volume_normalized=False)
        energy_dist = openmc.stats.Uniform(a=bounds[index], b=bounds[index + 1])
        source_list.append(openmc.IndependentSource(space=mesh_dist, energy=energy_dist, strength=np.sum(esd[mesh_index][:, index]), particle='photon', domains=cells))
    return source_list, total_mesh

# Define tallies
def tallies(total_mesh):
    particle_filter = openmc.ParticleFilter('photon')
    total_filter = openmc.MeshFilter(total_mesh)
    
    neutron_tally = openmc.Tally(tally_id=1, name="Neutron tally")
    neutron_tally.scores = ['flux', 'elastic', 'absorption']
    
    # Implementing filter for neutron tally through shells with material
    cell_filter = openmc.CellFilter(cells_w_mats)
    neutron_tally.filters = [cell_filter, Total_Filter]

    # Vitamin-J energy filter created to assign to the flux tally.
    energy_filter_flux = openmc.EnergyFilter.from_group_structure("VITAMIN-J-42")

    spectrum_tally = openmc.Tally(tally_id=2, name="Flux spectrum")
    # Implementing energy and cell filters for flux spectrum tally
    spectrum_tally.filters = [cell_filter, total_filter, energy_filter_flux, particle_filter]
    spectrum_tally.scores = ['flux']
    
    tall = openmc.Tallies([neutron_tally, spectrum_tally])
    tall.export_to_xml()
    return tall, neutron_tally, spectrum_tally, particle_filter, total_filter, cell_filter

def settings(source_list):
    sets = openmc.Settings()
    sets.batches = 10
    sets.inactive = 1
    sets.particles = 100000
    sets.source = Source_List
    sets.run_mode = 'fixed source'
    sets.export_to_xml()
    return sets

def plot_universe(Void, W_Shell, C_Shell, Cells):
    universe_ss = openmc.Universe(cells=Cells)
    Plot = universe_ss.plot(width=(2500.0, 2500.0), basis='xz',
              colors={Void: 'blue', W_Shell: 'red', C_Shell: 'green'}, legend=True)
    Plot.figure.savefig('Universe')
    return universe_ss

# Exporting materials, geometry, and tallies to .xml
def export_to_xml(filepath, element_1, element_2, inner_radius_W, outer_radius_W, inner_radius_C):
    OpenMC_SF = mat_lib(filepath)
    materials = []
    for material_id, element in enumerate(elements):
         materials.append(make_element(element, material_id+1, OpenMC_SF))
    OpenMC_Mat = all_mat(OpenMC_W, OpenMC_C)
    OpenMC_Geometry = make_spherical_shell(inner_radius_W, outer_radius_W, inner_radius_C, OpenMC_W, OpenMC_C)
    OpenMC_Source = make_source(OpenMC_Geometry[2], OpenMC_Geometry[3], OpenMC_Geometry[4])
    OpenMC_Settings = settings(OpenMC_Source[0])
    OpenMC_Tallies = tallies(OpenMC_Geometry[2], OpenMC_Geometry[3], OpenMC_Source[1], OpenMC_Source[2])
    OpenMC_Universe = plot_universe(OpenMC_Geometry[1], OpenMC_Geometry[2], OpenMC_Geometry[3], OpenMC_Geometry[4])
    #return OpenMC_Materials, OpenMC_Geometry, OpenMC_Source, OpenMC_Settings, OpenMC_Tallies, OpenMC_Universe
    return OpenMC_SF, OpenMC_W, OpenMC_C, *OpenMC_Geometry, *OpenMC_Source, OpenMC_Settings, *OpenMC_Tallies
    
Lib_Lines, M_1, M_2, geometry, Void, W_Shell, C_Shell, Cells, Source_List, Particle_Filter, Total_Mesh, tall, neutron_tally, spectrum_tally, Total_Filter, cell_filter, sets = export_to_xml(fp, E_1, E_2, R_W_1, R_W_2, R_C_1)

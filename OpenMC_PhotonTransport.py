# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 10:07:33 2024

@author: Anupama Rajendra
"""

import openmc
import argparse
import numpy as np

Mesh_File = 'OpenMC_Mesh.h5m'
Strengths = []
#Performing simulation for the first decay time
with open('Sum_1.txt', 'r') as Sum_1:
    Lines = Sum_1.readlines()
    strengths = [float(strength) for strength in Lines[0].split()]

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

#ALARA Element Library to read density for specified element:
# with open(""+fp+"") as ALARA_Lib:
with open("elelib.std") as ALARA_Lib:
    Lib_Lines = ALARA_Lib.readlines()

# Create materials & export to XML:
#Simulating tungsten shell:
  
def make_W(element_1):
    M_1 = openmc.Material(material_id=1, name=element_1)
    for line in Lib_Lines:
        if line.split()[0].lower() == element_1.lower():
            Density_M_1 = float(line.strip().split()[3])
            M_1.set_density('g/cm3', Density_M_1)
    M_1.add_element(element_1, 1.00)
    return M_1

def make_C(element_2):
    M_2 = openmc.Material(material_id=2, name=element_2)
    for line in Lib_Lines:
        if line.split()[0].lower() == element_2.lower():
        #if line.startswith(element_2.lower()):
            Density_M_2 = float(line.strip().split()[3])
            M_2.set_density('g/cm3', Density_M_2)
    M_2.add_element(element_2, 1.00)
    return M_2

def all_mat(M_1, M_2):
    all_materials = openmc.Materials([M_1, M_2])
    all_materials.cross_sections = '../fendl-3.2-hdf5/cross_sections.xml'
    all_materials.export_to_xml()
    
# Create geometry
#Spherical shell:
def make_spherical_shell(inner_radius_W, outer_radius_W, inner_radius_C, M_1, M_2):    
    S_W_1= openmc.Sphere(r=inner_radius_W) #sphere of radius 1000cm
    inside_W_sphere_1 = -S_W_1
    outside_W_sphere_1 = +S_W_1
    S_W_2 = openmc.Sphere(r=outer_radius_W, boundary_type='vacuum')
    inside_W_sphere_2 = -S_W_2
    outside_W_sphere_2 = +S_W_2
    S_W_3 = outside_W_sphere_1 & inside_W_sphere_2 #filled with specified material
      
    S_C_1= openmc.Sphere(r=inner_radius_C) #sphere of radius 995cm
    inside_C_sphere_1 = -S_C_1
    outside_C_sphere_1 = +S_C_1
    S_C_3 = outside_C_sphere_1 & inside_W_sphere_1 #filled with specified material    

    # Mapping materials to geometry:
    Void = openmc.Cell(fill=None, region = inside_C_sphere_1)
    W_Shell = openmc.Cell(fill=M_1, region=S_W_3)
    C_Shell = openmc.Cell(fill=M_2, region=S_C_3)
    Cells = [Void, W_Shell, C_Shell]
    geometry = openmc.Geometry(Cells)
    geometry.export_to_xml()
    return geometry, Void, W_Shell, C_Shell, Cells

#Define source:
def make_source(W_Shell, C_Shell, Cells):
    Source_List = []
    Particle_Filter = openmc.ParticleFilter('photon')
    Total_Mesh = openmc.UnstructuredMesh(Mesh_File, library='moab')
    Mesh_Dist = openmc.stats.MeshSpatial(Total_Mesh, volume_normalized=False)  
    for index, bound in enumerate(bounds[:-1]):
        Energy_Dist = openmc.stats.Uniform(a=bounds[index], b=bounds[index + 1])
        #Source strengths given by Strengths list at the top
        Source = Source_List.append(openmc.IndependentSource(space=Mesh_Dist, energy=Energy_Dist, strength=Strengths[index], particle='photon', domains=Cells))
    np.savetxt('Energy_List.txt', Energy_List, fmt='%s')
    return Source_List, Particle_Filter, Total_Mesh, Mesh_Dist

# Define tallies
def tallies(W_Shell, C_Shell, Particle_Filter, Total_Mesh):
    # global Total_Mesh
    # global Total_Filter
    Total_Filter = openmc.MeshFilter(Total_Mesh)
    
    neutron_tally = openmc.Tally(tally_id=1, name="Neutron tally")
    neutron_tally.scores = ['flux', 'elastic', 'absorption']
    
    # Implementing filter for neutron tally through W shell
    global cell_filter
    cell_filter = openmc.CellFilter([W_Shell, C_Shell])
    neutron_tally.filters = [cell_filter, Total_Filter]

    # Creating a tally to get the flux energy spectrum.
    # An energy filter is created to assign to the flux tally.
    energy_filter_flux = openmc.EnergyFilter.from_group_structure("VITAMIN-J-42")

    spectrum_tally = openmc.Tally(tally_id=2, name="Flux spectrum")
    # Implementing energy and cell filters for flux spectrum tally
    spectrum_tally.filters = [cell_filter, Total_Filter, energy_filter_flux, Particle_Filter]
    spectrum_tally.scores = ['flux']
    
    tall = openmc.Tallies([neutron_tally, spectrum_tally])
    tall.export_to_xml()
    return tall, neutron_tally, spectrum_tally, Total_Filter, cell_filter

# Assign simulation settings
def settings(Source_List):
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
def export_to_xml(element_1, element_2, inner_radius_W, outer_radius_W, inner_radius_C):
    OpenMC_W = make_W(element_1)
    OpenMC_C = make_C(element_2)
    OpenMC_Mat = all_mat(OpenMC_W, OpenMC_C)
    OpenMC_Geometry = make_spherical_shell(inner_radius_W, outer_radius_W, inner_radius_C, OpenMC_W, OpenMC_C)
    OpenMC_Source = make_source(OpenMC_Geometry[2], OpenMC_Geometry[3], OpenMC_Geometry[4])
    OpenMC_Settings = settings(OpenMC_Source[0])
    OpenMC_Tallies = tallies(OpenMC_Geometry[2], OpenMC_Geometry[3], OpenMC_Source[1], OpenMC_Source[2])
    OpenMC_Universe = plot_universe(OpenMC_Geometry[1], OpenMC_Geometry[2], OpenMC_Geometry[3], OpenMC_Geometry[4])
    #return OpenMC_Materials, OpenMC_Geometry, OpenMC_Source, OpenMC_Settings, OpenMC_Tallies, OpenMC_Universe
    return OpenMC_W, OpenMC_C, *OpenMC_Geometry, *OpenMC_Source, OpenMC_Settings, *OpenMC_Tallies

M_1, M_2, geometry, Void, W_Shell, C_Shell, Cells, Source_List, Particle_Filter, Total_Mesh, Mesh_Dist, tall, neutron_tally, spectrum_tally, Total_Filter, cell_filter, sets = export_to_xml(E_1, E_2, R_W_1, R_W_2, R_C_1)

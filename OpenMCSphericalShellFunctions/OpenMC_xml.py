# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 10:07:33 2024

@author: Anupama Rajendra
"""

import openmc
import argparse

parser = argparse.ArgumentParser(description="Specify required inputs: file path to ALARA Element Library, element name, inner radius [cm], outer_radius [cm], particle energy [eV]")

#Required user inputs:
parser.add_argument('--filepath', type=str, required=True)    
parser.add_argument('--element', type=str, required=True)
parser.add_argument('--inner_radius', type=float, required=True)
parser.add_argument('--outer_radius', type=float, required=True)
parser.add_argument('--energy', type=float, required=True)

args = parser.parse_args()

fp = args.filepath
E_1 = args.element
R_1 = args.inner_radius
R_2 = args.outer_radius
eg = args.energy

#ALARA Element Library to read density for specified element:
with open(""+fp+"") as ALARA_Lib:
    Lib_Lines = ALARA_Lib.readlines()

# Create materials & export to XML:
#Simulating tungsten shell:
print(E_1)   
def make_material(element):
    M_1 = openmc.Material()
    for line in Lib_Lines:
        if line.startswith(element.lower()):
            Density = float(line.strip().split()[3])
            M_1.set_density('g/cm3', Density)
    M_1.add_element(element, 1.00)
    all_materials = openmc.Materials([M_1])
    all_materials.export_to_xml()
    return M_1

# Create geometry
#Spherical shell:
def make_spherical_shell(inner_radius, outer_radius, M_1):    
    S_1= openmc.Sphere(r=R_1) #sphere of radius 1000cm
    inside_sphere_1 = -S_1
    outside_sphere_1 = +S_1
    S_2 = openmc.Sphere(r=R_2, boundary_type='vacuum')
    inside_sphere_2 = -S_2
    outside_sphere_2 = +S_2
    S_3 = outside_sphere_1 & inside_sphere_2 #filled with specified material

   # Mapping materials to geometry:
    Void = openmc.Cell(fill=None, region = inside_sphere_1)
    Shell = openmc.Cell(fill=M_1, region=S_3)
    Cells = [Void, Shell]
    geometry = openmc.Geometry(Cells)
    geometry.export_to_xml()
    return geometry, Void, Shell, Cells
    
# Source distribution:
#Point-source is hard-coded, particle energy is a required input    
def make_source(energy):
    PointSource = openmc.stats.Point(xyz=(0.0, 0.0, 0.0))
    Prob = openmc.stats.Discrete(energy, 1.0)
    return PointSource, Prob

# Assign simulation settings
def settings(PointSource, Prob):
    sets = openmc.Settings()
    sets.batches = 10
    sets.inactive = 1
    sets.particles = 10000
    sets.source = openmc.Source(space=PointSource, energy=Prob, strength = 1.0, particle = 'neutron')
    sets.run_mode = 'fixed source'
    sets.export_to_xml()
    return sets

# Define tallies
def tallies(Shell):
    neutron_tally = openmc.Tally(name="Neutron tally")
    neutron_tally.scores = ['flux', 'elastic', 'absorption']
    
    # Implementing filter for neutron tally through W shell
    cell_filter = openmc.CellFilter([Shell])
    neutron_tally.filters = [cell_filter]

    # Creating a tally to get the flux energy spectrum.
    # An energy filter is created to assign to the flux tally.
    energy_filter_flux = openmc.EnergyFilter.from_group_structure("VITAMIN-J-175")

    spectrum_tally = openmc.Tally(name="Flux spectrum")
    # Implementing energy and cell filters for flux spectrum tally
    spectrum_tally.filters = [energy_filter_flux, cell_filter]
    spectrum_tally.scores = ['flux']
    
    tall = openmc.Tallies([neutron_tally, spectrum_tally])
    tall.export_to_xml()
    return tall, neutron_tally, spectrum_tally

def plot_universe(Void, Shell, Cells):
    universe_ss = openmc.Universe(cells=Cells)
    Plot = universe_ss.plot(width=(2500.0, 2500.0), basis='xz',
              colors={Void: 'blue', Shell: 'red'}, legend=True)
    Plot.figure.savefig('Universe')
    return universe_ss

# Exporting materials, geometry, and tallies to .xml
def export_to_xml(element, inner_radius, outer_radius, energy):
    OpenMC_Materials = make_material(element)
    OpenMC_Geometry = make_spherical_shell(inner_radius, outer_radius, OpenMC_Materials)
    OpenMC_Energy = make_source(energy)
    OpenMC_Settings = settings(*OpenMC_Energy)
    OpenMC_Tallies = tallies(OpenMC_Geometry[2])
    OpenMC_Universe = plot_universe(OpenMC_Geometry[1], OpenMC_Geometry[2], OpenMC_Geometry[3])
    #return OpenMC_Materials, OpenMC_Geometry, OpenMC_Energy, OpenMC_Settings, OpenMC_Tallies, OpenMC_Universe
    return OpenMC_Materials, *OpenMC_Geometry, *OpenMC_Energy, OpenMC_Settings, *OpenMC_Tallies

M_1, geometry, Void, Shell, Cells, PointSource, Prob, sets, tall, neutron_tally, spectrum_tally = export_to_xml(E_1, R_1, R_2, eg)

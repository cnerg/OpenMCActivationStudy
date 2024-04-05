# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 08:16:12 2024

@author: Anupama Rajendra
"""
import openmc
import numpy as np
#import tkinter as tk

# Create materials & export to XML:
#Simulating tungsten shell:
W = openmc.Material(name='W_Shell')
W.set_density('g/cm3', 19.28)
W.add_element('W', 1.0)
materials = openmc.Materials([W])
materials.export_to_xml()

# Create geometry
#Spherical shell:
R_1= openmc.Sphere(r=1000) #sphere of radius 1000cm
inside_sphere_1 = -R_1
outside_sphere_1 = +R_1
R_2 = openmc.Sphere(r=1005, boundary_type='vacuum')
inside_sphere_2 = -R_2
outside_sphere_2 = +R_2
R_3 = outside_sphere_1 & inside_sphere_2 #filled with tungsten

# Mapping materials to geometry:
Void = openmc.Cell(fill=None, region = inside_sphere_1)
Shell = openmc.Cell(fill=W, region=R_3)
geometry = openmc.Geometry([Void, Shell])
geometry.export_to_xml()


# # Source distribution:
PointSource = openmc.stats.Point(xyz=(0.0, 0.0, 0.0))
Prob = openmc.stats.Discrete(14E+06, 1.0)

# Assign simulation settings
settings = openmc.Settings()
settings.batches = 10
settings.inactive = 1
settings.particles = 100000
settings.source = openmc.Source(space=PointSource, energy=Prob, strength = 10.0, particle = 'neutron')
settings.run_mode = 'fixed source'
settings.export_to_xml()

# Define tallies
neutron_tally = openmc.Tally(name="Neutron tally")
neutron_tally.scores = ['flux', 'elastic', 'absorption']
# Implementing filter for neutron tally through W shell
cell_filter = openmc.CellFilter([Shell])
neutron_tally.filters = [cell_filter]

# Creating a tally to get the flux energy spectrum.
# An energy filter is created to assign to the flux tally.
e_min, e_max = 5e2, 14.001e6
groups = 500
energies = np.logspace(np.log10(e_min), np.log10(e_max), groups + 1)
energy_filter = openmc.EnergyFilter(energies)

spectrum_tally = openmc.Tally(name="Flux spectrum")
# Implementing energy and cell filters for flux spectrum tally
spectrum_tally.filters = [energy_filter, cell_filter]
spectrum_tally.scores = ['flux']

# Collecting and exporting tallies to .xml
tallies = openmc.Tallies([neutron_tally, spectrum_tally])
tallies.export_to_xml()

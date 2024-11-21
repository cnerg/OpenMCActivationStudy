# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 15:56:57 2024

@author: Anupama Rajendra
"""

import openmc

# Create geometry
#Spherical shell:
def make_spherical_shells(layers, inner_radius):    
    '''
    Creates a set of concentric spherical shells, each with its own material & inner/outer radius.
    
    layers: list of tuples with OpenMC Material name and thickness: (material, thickness)
    inner_radius: the radius of the innermost spherical shell
    '''
    inner_sphere = openmc.Sphere(r = inner_radius)
    cells = [openmc.Cell(fill = None, region = -inner_sphere)]
    for (material, thickness) in layers:
        outer_radius = inner_radius + thickness
        outer_sphere = openmc.Sphere(r = outer_radius)
        cells.append(openmc.Cell(fill = material, region = +inner_sphere & -outer_sphere))
        outer_radius = inner_radius
        outer_sphere = inner_sphere
    outer_sphere.boundary_type = 'vacuum'     
    cells.append(openmc.Cell(fill = None, region = +outer_sphere)) 
    geometry = openmc.Geometry(cells)    
    return geometry

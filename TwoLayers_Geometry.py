# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 15:56:57 2024

@author: Anupama Rajendra
"""

import openmc

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
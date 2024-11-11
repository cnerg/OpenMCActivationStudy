# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 13:49:23 2024

@author: Anupama Rajendra
"""
import openmc

#ALARA Element Library to read density for specified element:
def alara_element_densities(filepath):  
    with open(""+filepath+"") as ALARA_Lib:
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
    
# Create materials & export to XML:
#Simulating tungsten shell:
def make_element(elements, density_dict):
    mats = openmc.Materials([])
    for element_id, element in enumerate(elements):
        mat = openmc.Material(material_id=element_id+1, name=element)
        mat.add_element(element, 1.00)
        mat.set_density('g/cm3', density_dict.get(element.lower()))
        mats.append(mat)
    return mats


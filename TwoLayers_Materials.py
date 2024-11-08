# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 13:49:23 2024

@author: Anupama Rajendra
"""
import openmc

#ALARA Element Library to read density for specified element:
def mat_lib(filepath):  
    with open(""+filepath+"") as ALARA_Lib:
        Lib_Lines = ALARA_Lib.readlines()
    return Lib_Lines
    
def create_densities(elements, Lib_Lines):
    density_dict = {}
    for element in elements:
        for line in Lib_Lines:
            if line.split()[0].lower() == element.lower():
                density_dict[element] = float(line.strip().split()[3])
    return density_dict  
    
# Create materials & export to XML:
#Simulating tungsten shell:
def make_element(elements, density_dict):
    mats = openmc.Materials([])
    for element_id, element in enumerate(elements):
        mat = openmc.Material(material_id=element_id+1, name=element)
        mat.add_element(element, 1.00)
        mat.set_density('g/cm3', density_dict.get(f'{element}'))
        mats.append(mat)
    return mats

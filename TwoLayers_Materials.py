# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 13:49:23 2024

@author: Anupama Rajendra
"""
import openmc

#ALARA Element Library to read density for specified element:
def mat_lib(filepath):  
    with open(""+filepath+"") as ALARA_Lib:
    #with open("elelib.std") as ALARA_Lib:
        Lib_Lines = ALARA_Lib.readlines()
    return Lib_Lines

# Create materials & export to XML:
#Simulating tungsten shell:
  
def make_W(element_1, Lib_Lines):
    M_1 = openmc.Material(material_id=1, name=element_1)
    for line in Lib_Lines:
        if line.split()[0].lower() == element_1.lower():
            Density_M_1 = float(line.strip().split()[3])
            M_1.set_density('g/cm3', Density_M_1)
    M_1.add_element(element_1, 1.00)
    return M_1

def make_C(element_2, Lib_Lines):
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
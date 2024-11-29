import openmc

def alara_element_densities(filepath):  
    '''
    Create a dictionary of element names and their corresponding densities using the ALARA element library.
    
    inputs:
        filepath: path to file containing ALARA element library
        
    '''
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

def make_element(element, density_dict):
    '''
    inputs:
        element: elemental symbol of chosen element (str)
        density_dict: dictionary with key = elemental symbol & value = density [g/cm^3]
        
    outputs:
        mats : OpenMC Materials object
    '''
    mat = openmc.Material(material_id=1, name=element)
    mat.add_element(element, 1.00)
    mat.set_density('g/cm3', density_dict.get(element.lower()))
    mats = openmc.Materials([mat])
    return mats
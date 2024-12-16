import h5py
import numpy as np
import yaml
import argparse

def extract_source_data(source_mesh_list, num_elements, photon_groups):
    '''
    Identifies the location of the source density dataset within each mesh file.
    
    input: iterable of .h5m filenames (str), whose files contain photon source information
    output: numpy array of source density data with rows = # of mesh elements and 
    columns = # number of photon groups, with one array per source mesh
    '''
    
    sd_list = np.ndarray((len(source_mesh_list), num_elements, photon_groups))
    for source_index, source_name in enumerate(source_mesh_list):
         file = h5py.File(source_name, 'r')
         sd_list[source_index,:] = file['tstt']['elements']['Tet4']['tags']['source_density'][:]
    return sd_list

def save_source_density(sd_list, sd_filename):
    '''
    Saves source density data as a text file.
    
    sd_list: list of source density datasets from h5m files
    sd_filename: user-provided filename that appears as a prefix for the text files (str)
    '''
    for sd_index, source_density in enumerate(sd_list):
        with open(f'{sd_filename}_{sd_index + 1}.txt', 'w') as source:
            for tet_element in source_density:
                source.write(' '.join(map(str, tet_element)) + '\n')
                
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--Mesh_Reader_YAML', required = True, help="Path (str) to YAML file containing required inputs for Source_Mesh_Reader")
    args = parser.parse_args()
    smr_yaml = args.Mesh_Reader_YAML
    
    with open(smr_yaml, 'r') as yaml_file:
        inputs = yaml.safe_load(yaml_file)
        
    source_mesh_list = inputs['source_meshes']
    num_elements = inputs['num_elements']
    photon_groups = inputs['photon_groups']
    sd_filename = inputs['sd_filename']
    
    esd = extract_source_data(source_mesh_list, num_elements, photon_groups)
    ssd = save_source_density(esd, sd_filename)

if __name__ == "__main__":
    main()

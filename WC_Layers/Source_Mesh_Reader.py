import h5py
import numpy as np
import yaml
import argparse
    
def extract_source_data(source_mesh_list, num_elements, num_photon_groups):
    '''
    Identifies the location of the source density dataset within each mesh file.
    
    input: 
        source_mesh_list: iterable of .h5m filenames (str), whose files contain photon source information
    output: numpy array of source density data with rows = # of mesh elements and columns = # number of photon groups, with one array per source mesh
    ''' 
    sd_list = np.ndarray((len(source_mesh_list), num_elements, num_photon_groups))
    for source_index, source_name in enumerate(source_mesh_list):
         file = h5py.File(source_name, 'r')
         # directly accessing the data stored in the HDF5 as specified by the MOAB mesh library
         # MOAB data is stored in the `tstt` group, and this data is in a `tag` with name `source_density`
         # applied to each of the `Tet4` type of `elements`
         sd_list[source_index,:] = file['tstt']['elements']['Tet4']['tags']['source_density'][:]
    return sd_list  

def save_source_density(sd_list, sd_filename):
    '''
    Saves source density data as a separate text file for each source mesh. Each file is of size num_photon_groups * num_elements.
    inputs:
        sd_list: list of source density datasets from h5m files
        sd_filename: user-provided filename that appears as a prefix for the text files (str)
    '''
    for sd_index, source_density in enumerate(sd_list):
        with open(f'{sd_filename}_{sd_index + 1}.txt', 'w') as source:
            for tet_element in source_density:
                source.write(' '.join(map(str, tet_element)) + '\n')

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--Mesh_Reader_YAML', default = 'Source_Mesh_Reader_Inputs.yaml', help="Path (str) to YAML file containing required inputs for Source_Mesh_Reader")
    args = parser.parse_args()
    return args
    
def read_yaml(args): 
    with open(args.Mesh_Reader_YAML, 'r') as yaml_file:
        mesh_reader_inputs = yaml.safe_load(yaml_file)
    return mesh_reader_inputs

def read_source_mesh(mesh_reader_inputs):
    #Find the size of the first source density dataset (assumed to be the same for all other datasets):
    sd_data = h5py.File(mesh_reader_inputs['source_meshes'][0], 'r')['tstt']['elements']['Tet4']['tags']['source_density'][:]
    sd_list = extract_source_data(mesh_reader_inputs['source_meshes'],
                                      sd_data.shape[0],
                                      sd_data.shape[1])
    source_density = save_source_density(sd_list, 
                              mesh_reader_inputs['sd_filename'])

def main():
    args = parse_args()
    mesh_reader_inputs = read_yaml(args)
    read_source_mesh(mesh_reader_inputs)

if __name__ == "__main__":
    main()

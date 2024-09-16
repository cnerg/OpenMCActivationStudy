import openmc
import numpy as np
import pymoab
from pymoab import core, types, rng
import matplotlib.pyplot as plt

# Load the statepoint file
sp = openmc.StatePoint("statepoint.10.h5")

mb = core.Core()
mb.load_file('OpenMC_Mesh.h5m')

# Retrieve the tally
tally = sp.get_tally(id=2)

# Load the mesh
mesh = sp.meshes[1]

tally_data = tally.get_reshaped_data(value='mean')

#Summing over cell filter:
flux_sum_en = tally_data.sum(axis=0)

ebounds = tally.find_filter(openmc.EnergyFilter).bins

plt.xlabel('Energy [eV]')
plt.ylabel('Flux [n/cm^2-s]')
plt.loglog(ebounds[:,0],flux_sum_en[0,:,0,0])
plt.savefig("Flux_Graph.png")

# get all tets from the MOAB mesh
all_tets = mb.get_entities_by_type(0, types.MBTET)
n_flux_tag = mb.tag_get_handle('FLUX_MESH', 175, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, create_if_missing=True)
mb.tag_set_data(n_flux_tag, all_tets, flux_sum_en[:,:,0,0]) # the shape of the data will need to be checked here

mb.write_file('Mesh_with_Tally.h5m')

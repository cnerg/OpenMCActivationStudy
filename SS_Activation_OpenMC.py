# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 08:16:12 2024

@author: Anupama Rajendra
"""
import openmc
import openmc.deplete
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import random

# Importing Vitamin-J energy group structure:
# This excel file contains the energy bounds of the Vitamin J structure
# Vit_J = pd.read_excel('VitaminJEnergyGroupStructure.xlsx')
# ebounds = Vit_J.iloc[:, 1]

openmc.config['chain_file'] = 'chain_endfb71_sfr.xml'

# Create materials & export to XML:
#Simulating tungsten shell:
W = openmc.Material(material_id=1, name='W_Shell')
W.set_density('g/cm3', 19.28)
W.add_element('W', 1.00)
materials = openmc.Materials([W])
materials.export_to_xml()

# Create geometry
#Spherical shell:
    
R_1 = 1000
S_1= openmc.Sphere(r=R_1) #sphere of radius 1000cm
inside_sphere_1 = -S_1
outside_sphere_1 = +S_1
R_2 = 1005
S_2 = openmc.Sphere(r=R_2, boundary_type='vacuum')
inside_sphere_2 = -S_2
outside_sphere_2 = +S_2
S_3 = outside_sphere_1 & inside_sphere_2 #filled with tungsten

# Mapping materials to geometry:
Void = openmc.Cell(fill=None, region = inside_sphere_1)
Shell = openmc.Cell(fill=W, region=S_3)
geometry = openmc.Geometry([Void, Shell])
geometry.export_to_xml()

# Source distribution:
PointSource = openmc.stats.Point(xyz=(0.0, 0.0, 0.0))
Prob = openmc.stats.Discrete(14E+06, 1.0)

# Assign simulation settings
settings = openmc.Settings()
settings.batches = 10
settings.inactive = 1
settings.particles = 10000
settings.source = openmc.Source(space=PointSource, energy=Prob, strength = 1.0, particle = 'neutron')
settings.run_mode = 'fixed source'
settings.export_to_xml()

# Define tallies
neutron_tally = openmc.Tally(name="Neutron tally")
neutron_tally.scores = ['flux', 'elastic', 'absorption']
# Implementing filter for neutron tally through W shell
cell_filter = openmc.CellFilter([Shell])
neutron_tally.filters = [cell_filter]

# Creating a tally to get the flux energy spectrum.
# An energy filter is created to assign to the flux tally using the Vitamin J structure.
energy_filter_flux = openmc.EnergyFilter.from_group_structure("VITAMIN-J-175")

spectrum_tally = openmc.Tally(name="Flux spectrum")
# Implementing energy and cell filters for flux spectrum tally
spectrum_tally.filters = [energy_filter_flux, cell_filter]
spectrum_tally.scores = ['flux']

# Collecting and exporting tallies to .xml
tallies = openmc.Tallies([neutron_tally, spectrum_tally])
tallies.export_to_xml()

model = openmc.model.Model(geometry=geometry,settings=settings)
#Depletion calculation
W.depletable = True
W.volume = 4.0/3.0 * np.pi * (R_2**3 - R_1**3) #volume of W wall material
operator = openmc.deplete.CoupledOperator(model, normalization_mode='source-rate')
time_steps = [3e8, 86400, 2.6e6]
source_rates = [1E+18, 0, 0]
integrator = openmc.deplete.PredictorIntegrator(operator=operator, timesteps=time_steps, source_rates=source_rates, timestep_units='s')
integrator.integrate()

#Opening statepoint file to read tallies:
with openmc.StatePoint('statepoint.10.h5') as sp:
    fl = sp.get_tally(name="Flux spectrum")
    nt = sp.get_tally(name="Neutron tally")
    
# Get the neutron energies from the energy filter
energy_filter_fl = fl.filters[0]
energies_fl = energy_filter_fl.bins[:, 0]

# Get the neutron flux values
flux = fl.get_values(value='mean').ravel()

#Neutron flux/elastic/absorption tallies:
tal = nt.get_values(value='mean').ravel()
print(tal)

Flux_Data = np.c_[energies_fl, flux]
#Creating an excel file that stores flux data for each energy bin (used as input for ALARA)
FD_Excel = pd.DataFrame(Flux_Data, columns=['Energy [eV]', 'Flux [n-cm/sp]'])
FD_Excel.to_excel('Neutron_Flux.xlsx', index=False)

Tallies_Excel = pd.DataFrame(tal)
#Creating an excel file that stores total tally value data
Tallies_Excel.to_excel('Tally_Values.xlsx', index=False)

# Depletion results file
results = openmc.deplete.Results(filename='depletion_results.h5')

# Stable W nuclides present at beginning of operation (will not be plotted)
stable_nuc = ['W180', 'W182', 'W183', 'W184', 'W186']  

#Store list of nuclides from last timestep as a Materials object
materials_last = results.export_to_materials(-1)
# Storing depletion data from 1st material
mat_dep = materials_last[0]
# Obtaining the list of nuclides from the results file
nuc_last = mat_dep.get_nuclides()

# Removing stable W nuclides from list so that they do not appear in the plot
for j in stable_nuc :
    nuc_last.remove(j)
print(nuc_last)

colors = list(mcolors.CSS4_COLORS.keys())
num_dens= {}
pair_list = {}

for nuclide in nuc_last:
    plot_color = random.choice(colors)
    time, num_dens[nuclide] = results.get_atoms('1', nuclide, nuc_units = 'atom/cm3')
    print(time, num_dens[nuclide])
    plt.plot(time, num_dens[nuclide], marker='.', linestyle='solid', color=plot_color, label=nuclide)

# Adding labels and title
plt.xlabel('Time after shutdown [s]')
plt.xlim(time_steps[1], time_steps[-1])
#plt.gca().set_ylim(bottom=0)
plt.ylabel('Nuclide density [atoms/cm^3]')
plt.xscale("log")
plt.yscale("log")
plt.title('Plot of number density vs time after shutdown')

# Adding a legend
plt.legend()

plt.savefig('Nuclide_density_OpenMC')
# Display the plot
plt.show()

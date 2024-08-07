# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 12:34:26 2024

@author: Anupama Rajendra
"""

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
from matplotlib.pyplot import figure

openmc.config['chain_file'] = 'chain_endfb71_sfr.xml'

#Importing Vitamin-J energy group structure:
Vit_J = pd.read_excel('VitaminJEnergyGroupStructure.xlsx')
ebounds = Vit_J.iloc[:, 1]

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
Cells = [Void, Shell]
geometry = openmc.Geometry(Cells)
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
# An energy filter is created to assign to the flux tally.
energy_filter_flux = openmc.EnergyFilter.from_group_structure("VITAMIN-J-175")

spectrum_tally = openmc.Tally(name="Flux spectrum")
# Implementing energy and cell filters for flux spectrum tally
spectrum_tally.filters = [energy_filter_flux, cell_filter]
spectrum_tally.scores = ['flux']

universe_ss = openmc.Universe(cells=Cells)
Plot = universe_ss.plot(width=(2500.0, 2500.0), basis='xz',
              colors={Void: 'blue', Shell: 'red'}, legend=True)
Plot.figure.savefig('Universe.png')

# Collecting and exporting tallies to .xml
tallies = openmc.Tallies([neutron_tally, spectrum_tally])
tallies.export_to_xml()

#----------------------------------------------------------------------------
#Beginning of the depletion model:

model = openmc.model.Model(geometry=geometry,settings=settings)
#Depletion calculation
W.depletable = True
W.volume = 4.0/3.0 * np.pi * (R_2**3 - R_1**3) #volume of W wall material
operator = openmc.deplete.CoupledOperator(model, normalization_mode='source-rate')
time_steps = [3E+8, 86400, 2.6E+6]
source_rates = [1E+18, 0,0]
integrator = openmc.deplete.PredictorIntegrator(operator, time_steps, source_rates=source_rates, timestep_units='s')
integrator.integrate()

#Statepoint file to read tallies:
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

#Saving neutron flux values to csv file
#Dividing by volume to obtain proper units of flux (#/cm^2-s)
Flux_Data = np.c_[energies_fl, flux/W.volume]
#print(Flux_Data[1])
#ALARA flux inputs go from high energy to low energy
Flux_Data_ALARA = Flux_Data[::-1]
FD_CSV = pd.DataFrame(Flux_Data_ALARA, columns=['Energy [eV]', 'Flux [n-cm/sp]'])
FD_CSV.to_csv('Neutron_Flux.csv', index=False)

Tallies_CSV = pd.DataFrame(tal)
Tallies_CSV.to_csv('Tally_Values.csv')
                                           
# Depletion results file
results = openmc.deplete.Results(filename='depletion_results.h5')

# Stable W nuclides present at beginning of operation:
stable_nuc = ['W180', 'W182', 'W183', 'W184', 'W186']    

materials_list=[]
#Store list of nuclides from each timestep as a Materials object
for step in range(len(time_steps)):
    materials_object = results.export_to_materials(step)
    # Storing depletion data from 1st material
    mat_dep = materials_object[0]
    # Obtaining the list of nuclides in the results file
    rad_nuc = mat_dep.get_nuclides()
    materials_list.append(rad_nuc)

#Removing stable W nuclides from list so that they do not appear in the plot
for j in stable_nuc :
    rad_nuc.remove(j)

print(rad_nuc)
print(len(rad_nuc))    

colors = list(mcolors.CSS4_COLORS.keys())
num_dens= {}
pair_list = {}
g_list = {}

plt.figure(figsize=(9,10))

with open(r'Densities_CSV.csv', 'w') as density_file:

    for nuclide in rad_nuc:
        plot_color = random.choice(colors)
        time, num_dens[nuclide] = results.get_atoms('1', nuclide, nuc_units = 'atom/cm3')
        print(time, num_dens[nuclide])
        density_file.write(f'{nuclide},')
        density_file.write(','.join(map(str, num_dens[nuclide])) + '\n')
        plt.plot(time, num_dens[nuclide], marker='.', linestyle='solid', color=plot_color, label=nuclide)


# Adding labels and title
plt.xlabel('Time after start of operation [s]')
plt.xlim(1, sum(time_steps))
#plt.ylim(1e-09, 1e+20)
#plt.gca().set_ylim(bottom=0)
plt.ylabel('Nuclide density [atoms/cm^3]')
plt.xscale("log")
plt.yscale("log")
plt.title('Plot of number density vs time after operation')
# Adding a legend
plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.15), fontsize='8')

plt.savefig('Nuclide_density_OpenMC')
plt.close()
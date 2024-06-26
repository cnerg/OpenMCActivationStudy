# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 08:16:12 2024

@author: Anupama Rajendra
"""
import openmc
import openmc.deplete
import numpy as np
import matplotlib.pyplot as plt

# Importing Vitamin-J energy group structure:
# The text file contains the energy bounds of the Vitamin J structure
with open('VitJ.txt', 'r') as ebounds:
   VitJ  = ebounds.readlines()

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
energy_filter_flux = openmc.EnergyFilter(VitJ)

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

# Depletion results file
results = openmc.deplete.Results(filename='depletion_results.h5')
# Obtain U235 concentration as a function of time

time, n_He4 = results.get_atoms('1', 'He4', nuc_units = 'atom/cm3')
time, n_W179 = results.get_atoms('1', 'W179', nuc_units = 'atom/cm3')
time, n_W180 = results.get_atoms('1', 'W180', nuc_units = 'atom/cm3')
time, n_W181 = results.get_atoms('1', 'W181', nuc_units = 'atom/cm3')
time, n_W182 = results.get_atoms('1', 'W182', nuc_units = 'atom/cm3')
time, n_W183 = results.get_atoms('1', 'W183', nuc_units = 'atom/cm3')
time, n_W184 = results.get_atoms('1', 'W184', nuc_units = 'atom/cm3')
time, n_W185 = results.get_atoms('1', 'W185', nuc_units = 'atom/cm3')
time, n_W186 = results.get_atoms('1', 'W186', nuc_units = 'atom/cm3')
time, n_W187 = results.get_atoms('1', 'W187', nuc_units = 'atom/cm3')
time, n_Re185 = results.get_atoms('1', 'Re185', nuc_units = 'atom/cm3')
time, n_Re187 = results.get_atoms('1', 'Re187', nuc_units = 'atom/cm3')
time, n_Os187 = results.get_atoms('1', 'Os187', nuc_units = 'atom/cm3')

# time, n_He4 = results.get_atoms('1', 'He4')
# time, n_W180 = results.get_atoms('1', 'W180')
# time, n_W181 = results.get_atoms('1', 'W181')
# time, n_W182 = results.get_atoms('1', 'W182')
# time, n_W183 = results.get_atoms('1', 'W183')
# time, n_W184 = results.get_atoms('1', 'W184')
# time, n_W185 = results.get_atoms('1', 'W185')
# time, n_W186 = results.get_atoms('1', 'W186')
# time, n_W187 = results.get_atoms('1', 'W187')
# time, n_Re185 = results.get_atoms('1', 'Re185')
# time, n_Re187 = results.get_atoms('1', 'Re187')
# time, n_Os187 = results.get_atoms('1', 'Os187')

print(time, n_He4)
print(time, n_W179)
print(time, n_W180)
print(time, n_W181)
print(time, n_W182)
print(time, n_W183)
print(time, n_W184)
print(time, n_W185)
print(time, n_W186)
print(time, n_W187)
print(time, n_Re185)
print(time, n_Re187)
print(time, n_Os187)

# Plotting the arrays with customizations
plt.plot(time, n_He4, marker='.', linestyle='solid', color='red', label='He4')
plt.plot(time, n_W179, marker='.', linestyle='solid', color='blue', label="W179")
plt.plot(time, n_W180, marker='.', linestyle='solid', color='green', label="W180")
plt.plot(time, n_W181, marker='.', linestyle='solid', color='deepskyblue', label="W181")
plt.plot(time, n_W182, marker='.', linestyle='solid', color='peru', label="W182")
plt.plot(time, n_W183, marker='.', linestyle='solid', color='plum', label="W183")
plt.plot(time, n_W184, marker='.', linestyle='solid', color='crimson', label="W184")
plt.plot(time, n_W185, marker='.', linestyle='solid', color='rosybrown', label="W185")
plt.plot(time, n_W186, marker='.', linestyle='solid', color='gray', label="W186")
plt.plot(time, n_W187, marker='.', linestyle='solid', color='lime', label="W187")
plt.plot(time, n_Re185, marker='.', linestyle='solid', color='gold', label="Re185")
plt.plot(time, n_Re187, marker='.', linestyle='solid', color='tomato', label="Re187")
plt.plot(time, n_Os187, marker='.', linestyle='solid', color='cyan', label="Os187")

# Adding labels and title
plt.xlabel('Time after shutdown [s]')
plt.ylabel('Nuclide density [atoms/cm^3]')
plt.title('Plot of number density vs time after shutdown')

# Adding a legend
plt.legend()

plt.savefig('Nuclide_density_OpenMC')
# Display the plot
plt.show()

#DH = results.get_decay_heat('1', units='W/cm3', by_nuclide=False)

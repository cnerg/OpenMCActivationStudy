import matplotlib.pyplot as plt
import openmc


# Get results from statepoint
with openmc.StatePoint('statepoint.100.h5') as sp:
    t = sp.get_tally(name="Flux spectrum")
    k = sp.get_tally(name="Neutron tally")

    # Get the energies from the energy filter
    energy_filter = t.filters[0]
    energies = energy_filter.bins[:, 0]

    # Get the flux values
    mean = t.get_values(value='mean').ravel()
    
    #Flux/elastic/absorption tallies:
    tal = k.get_values(value='mean').ravel()
    print(tal)
    
# Plot flux spectrum
fix, ax = plt.subplots()
ax.loglog(energies, mean, drawstyle='steps-post')
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('Flux [neutron-cm/source]')
ax.grid(True, which='both')
plt.savefig('Neutron_flux_vs_energy.png')
plt.show()

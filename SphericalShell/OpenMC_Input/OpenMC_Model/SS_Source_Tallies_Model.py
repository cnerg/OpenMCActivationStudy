import openmc
    
def make_source(energy):
    point_source = openmc.stats.Point(xyz=(0.0, 0.0, 0.0))
    energy_dist = openmc.stats.Discrete(energy, 1.0)
    source = openmc.Source(space = point_source, energy = energy_dist, strength = 1.0, particle = 'neutron')
    return source

def make_settings(source, total_batches, inactive_batches, num_particles, run_mode):
    sets = openmc.Settings()
    sets.batches = total_batches
    sets.inactive = inactive_batches
    sets.particles = num_particles
    sets.source = source
    sets.run_mode = run_mode
    return sets

def make_tallies(tallied_cells):
    particle_filter = openmc.ParticleFilter('neutron')
    cell_filter = openmc.CellFilter(tallied_cells)
    energy_filter_flux = openmc.EnergyFilter.from_group_structure("VITAMIN-J-175")
    
    neutron_tally = openmc.Tally(name="Neutron tally")
    neutron_tally.scores = ['flux', 'elastic', 'absorption']
    neutron_tally.filters = [cell_filter, particle_filter]

    spectrum_tally = openmc.Tally(name="Neutron flux spectrum")
    spectrum_tally.filters = [energy_filter_flux, cell_filter, particle_filter]
    spectrum_tally.scores = ['flux']
    
    talls = openmc.Tallies([neutron_tally, spectrum_tally])
    return talls

def create_openmc_model(mats_object, geom_object, tallies_object, settings_object):
    model = openmc.model.Model(geometry = geom_object, materials = mats_object, settings = settings_object, tallies = tallies_object)
    return model

import openmc
import openmc.deplete
import numpy as np

def deplete_ss(chain_file_path, model, inner_radius, thickness, times_post_boc, particle_source_rates, norm_mode, timestep_units):
    material = model.materials[0]
    chain_file = chain_file_path
    material.depletable = True
    timesteps = times_post_boc
    source_rates = particle_source_rates
    operator = openmc.deplete.CoupledOperator(model, chain_file, normalization_mode = norm_mode)
    integrator = openmc.deplete.PredictorIntegrator(operator, timesteps, source_rates = source_rates, timestep_units = timestep_units)
    return integrator

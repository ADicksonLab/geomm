import numpy as np

def free_energy(weights, max_energy=100, zero_point_energy=1.0e-12):

    # assume the weights are normalized
    # -log(weights)
    free_energies = np.negative(np.log(weights))

    # set the min as the 0 energy
    free_energies = free_energies - free_energies.min()

    # then increase everything by the zero-point energy so there are
    # no zero energy states
    free_energies += zero_point_energy

    # energies greater than the max are set to the maximum
    free_energies[np.where(free_energies > max_energy)] = max_energy

    return free_energies

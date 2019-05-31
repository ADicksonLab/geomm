"""Methods for calculating temperatures from the velocities of atoms."""

import numpy as np

BOLTZMANN_CONSTANT = 1.3806504e-23

def kinetic_energy(velocities, masses):

    return 0.5 * sum(masses * np.array([np.dot(vel, vel) for vel in velocities]))

def temperature(velocities, masses):

    return (2. / (3. * BOLTZMANN_CONSTANT)) * kinetic_energy(velocities, masses)

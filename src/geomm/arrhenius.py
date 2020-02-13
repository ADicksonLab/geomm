import numpy as np

import pint

unit = pint.UnitRegistry()

BOLTZMANN_CONSTANT = 1.3806504e-23 * (unit.joule / unit.kelvin)
PLANCK_CONSTANT = 6.626070e-34 * (unit.joule * unit.second)
UNIVERSAL_GAS_CONSTANT = 8.314 * (unit.joule / (unit.kelvin * unit.mole))
TRANSMISSION_COEFFICIENT = 1.0

def arrhenius_preexponent(temperature):

    preexponent = (TRANSMISSION_COEFFICIENT * BOLTZMANN_CONSTANT * temperature) /\
                  PLANCK_CONSTANT

    return preexponent

def arrhenius_rate(activation_free_energy, temperature):

    preexponent = arrhenius_preexponent(temperature)

    exponent = -1 * (activation_free_energy / (UNIVERSAL_GAS_CONSTANT * temperature))

    rate = preexponent * np.exp(exponent)

    return rate

def arrhenius_activation_free_energy(rate, temperature):

    preexponent = arrhenius_preexponent(temperature)

    free_energy = -(UNIVERSAL_GAS_CONSTANT * temperature) * np.log(rate/preexponent)

    return free_energy

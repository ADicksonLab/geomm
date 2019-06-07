import numpy as np

def free_energy(weights, max_energy=None, zero_point_energy=1.0e-12,
                supermax_value=np.nan):
    """Perform a transformation of probability weights to free energies.

    The free energy of probabilities is a transformation of -ln(p)
    where p is the properly normalized probability.

    In addition to this simple transformation we set a maximum free
    energy (max_energy) and a mask value (supermax_value) for values
    greater than the max energy. This is to account for missing
    observations or extreme outliers. This can be disabled by setting
    max_energy to None.

    Free energy values are relative to a zero-point energy, which in
    the absence of a theoretically determined value can be simply
    determined by the minimum value in the dataset. If the
    zero_point_energy value is given this simply fixes this minimum
    free energy to the zero-point energy so that the rest of the
    values have some meaningful value. Typically, we set this to
    something to close to but greater than 0. Because a 0 free energy
    is equivalent to a probability of 1.0, which would make the
    dataset non-normalized.

    Parameters
    ----------

    weights : arraylike of float
        The normalized weights to transform

    max_energy : float or None
        Sets a maximum value for free energies. Free energies larger
        than this will be replaced by the 'supermax_value'. If None
        the masking is not performed.
       (Default = None)

    supermax_value : numpy scalar
        The value that replaces free energies that are larger than the
        max energy. Typically is a numpy NaN so as to omit the data.
       (Default = np.nan)

    zero_point_energy : float or None
        Free energies are relative and the minimum value will be set
        to this zero point energy. This is to avoid having 0 free
        energy values. If None this is not done.
       (Default = 1.0e-12)

    Returns
    -------

    free_energies : arraylike
       The free energy for each weight.

    """

    # assume the weights are normalized
    # -log(weights)
    free_energies = np.negative(np.log(weights))

    if zero_point_energy is not None:
        # set the min as the 0 energy
        free_energies = free_energies - free_energies.min()

        # then increase everything by the zero-point energy so there are
        # no zero energy states
        free_energies += zero_point_energy

    # energies greater than the max are set to the supermax_value
    if max_energy is not None:
        free_energies[np.where(free_energies > max_energy)] = supermax_value

    return free_energies

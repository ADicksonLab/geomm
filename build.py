import os
import sys

from Cython.Build import cythonize

import numpy as np

# use cythonize to build the extensions
modules = ["geomm/pyqcprot.pyx",]

extensions = cythonize(modules)

def build(setup_kwargs):
    """Needed for the poetry building interface."""

    setup_kwargs.update({
        'ext_modules' : extensions,
        'include_dirs' : [np.get_include()],
    })

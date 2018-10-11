import os
import sys

from distutils.core import Extension

from Cython.Build import cythonize

# use cythonize to build the extensions
modules = ["geomm/pyqcprot.pyx",]

extensions = cythonize(modules)

def build(setup_kwargs):
    """Needed for the poetry building interface."""

    setup_kwargs.update({
        'ext_modules' : extensions,
    })

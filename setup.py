#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import absolute_import
from __future__ import print_function

import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import setup, find_packages

import itertools as it

import versioneer

# package specific imports
from Cython.Build import cythonize
import numpy as np

# the basic needed requirements for a package
base_requirements = [
    'numpy',
    'scipy',
    'pint',
]

# extras requirements list

# SNIPPET: example extra requirement
# example_extra_requirements = ['requests']
# extras = [example_extra_requirements,]

# Add your extra requirements lists here:
extras = [
]

# combination of all the extras requirements
_all_requirements = [[base_requirements]] + extras
all_requirements = it.chain.from_iterable(_all_requirements)

setup(
    name='geomm',
    version=versioneer.get_version(),
    author="Samuel D. Lotz",
    author_email="samuel.lotz@salotz.info",
    description="A simple no-nonsense library for computing common geometry on macromolecular systems.",
    #long_description=open('README.org').read(),
    license="MIT",
    url="https://github.com/ADicksonLab/geomm",
    classifiers=[
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
        'Programming Language :: Python :: 3'
    ],

    keywords='geometry chemistry macromolecules protein structural-informatics',
    # building/dev
    setup_requires=[
        'pytest-runner',
        'numpy',
        'cython',
    ],
    tests_require=['pytest', 'tox'],

    cmdclass=versioneer.get_cmdclass(),

    include_dirs=[np.get_include()],
    ext_modules = cythonize("src/geomm/pyqcprot.pyx"),


    # package
    packages=find_packages(where='src'),

    package_dir={'' : 'src'},

    # if this is true then the package_data won't be included in the
    # dist. Use MANIFEST.in for this
    include_package_data=True,

    # pymodules is for single file standalone modules not part of the
    # package
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],

    install_requires=base_requirements,

    # SNIPPET: example of using extra requirements
    # extras_require={
    #     'extras' : example_extra_requirements
    #     'all' : all_requirements,
    # }

    # include your extra requirement sets here
    extras_require={
        'all' : all_requirements,
    }
)

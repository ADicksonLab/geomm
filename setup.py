#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import absolute_import
from __future__ import print_function

from setuptools import setup, find_packages
from Cython.Build import cythonize

import numpy as np

import itertools as it

# setuptools only specifies abstract requirements. For the concrete
# requirements i.e. index or repo URL see requirements.txt
abstract_requirements = [
    'numpy',
    'scipy',
]

proj_urls = {
    'Source' : 'https://github.com/ADicksonLab/geomm',
    'Tracker' : 'https://github.com/ADicksonLab/geomm/issues'
}

setup(
    name='geomm',
    version='0.1.7',
    author="Samuel D. Lotz",
    author_email="samuel.lotz@salotz.info",
    description="A simple no-nonsense library for computing common geometry on macromolecular systems.",
    #long_description=open('README.org').read(),
    license="MIT",
    url="https://github.com/ADicksonLab/geomm",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
        'Programming Language :: Python :: 3'
    ],
    keywords='geometry chemistry macromolecules protein structural-informatics',
    project_urls=proj_urls,
    # building/dev
    setup_requires=['pytest-runner', 'numpy', 'cython'],
    tests_require=['pytest', 'tox'],
    install_requires=abstract_requirements,
    include_dirs=[np.get_include()],
    ext_modules = cythonize("src/geomm/pyqcprot.pyx"),

    # package
    packages=find_packages('src'),

    package_dir={'' : 'src'},

    # if this is true then the package_data won't be included in the
    # dist, and I prefer this to MANIFEST
    include_package_data=False,


)

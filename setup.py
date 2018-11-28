from setuptools import setup, find_packages
from Cython.Build import cythonize

import numpy as np

setup(
    name='geomm',
    version='0.1.6',
    py_modules=['geomm'],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'scipy'
    ],
    include_dirs=[np.get_include()],
    ext_modules = cythonize("geomm/pyqcprot.pyx")
)

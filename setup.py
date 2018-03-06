from setuptools import setup, find_packages
from Cython.Build import cythonize

setup(
    name='geomm',
    version='0.1',
    py_modules=['geomm'],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
    ],
    ext_modules = cythonize("geomm/pyqcprot.pyx")
)

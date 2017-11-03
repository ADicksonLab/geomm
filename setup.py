from setuptools import setup, find_packages

setup(
    name='geomm',
    version='0.1',
    py_modules=['geomm'],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
    ],
)

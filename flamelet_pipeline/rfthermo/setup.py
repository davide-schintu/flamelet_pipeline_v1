from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import os

# Questa classe serve a trovare automaticamente gli header di pybind11
class get_pybind_include(object):
    def __str__(self):
        import pybind11
        return pybind11.get_include()

ext_modules = [
    Extension(
        'RFThermo',
        # Elenco dei file sorgente (aggiungi Mixture.cpp e utilities.cpp)
        ['src/bindings.cpp', 'src/Mixture.cpp', 'src/utilities.cpp', 'src/functions.cpp', 'src/mixingRules.cpp'], 
        include_dirs=[
            get_pybind_include(),
            'include'  # La tua cartella con gli header .h
        ],
        language='c++'
    ),
]

setup(
    name='RFThermo',
    version='0.1.0',
    author='Il Tuo Nome',
    description='Wrapper Python per RFThermoFunctions',
    ext_modules=ext_modules,
    install_requires=['pybind11>=2.5.0'],
    zip_safe=False,
)
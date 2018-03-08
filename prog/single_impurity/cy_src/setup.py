#!/opt/intel/intelpython2/bin/python

from distutils.core import setup
from Cython.Build import cythonize

setup(
        ext_modules = cythonize("algo.pyx")        
)

'''
Setup file for the spherical coordinate voxel traversal algorithm.
Requires Python 3.
Dependencies:
  1. Cython (https://cython.org/)
  2. Numpy (https://numpy.org/)
  3. distutils (https://docs.python.org/3/library/distutils.html)

Code must be compiled before use:
  > python3 cython_SVR_setup.py build_ext --inplace
'''

import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension(
    name="cython_SVR",
    sources=["cython_SVR.pyx", "../spherical_volume_rendering_util.cpp"],
    language="c++",
    extra_compile_args=["-std=c++11", "-O3", "-march=native", "-flto", "-fno-signed-zeros", "-funroll-loops"],
    define_macros = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')], # Hides deprecated Numpy warning.
    include_dirs = [numpy.get_include()],
)]

# Python-3
for e in ext_modules:
    e.cython_directives = {'language_level': "3"}

setup(
    name = 'cython_SVR',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
)
'''
Setup file for the spherical coordinate voxel traversal algorithm.

Code must be compiled before use:
> python Cython_SVR_setup.py build_ext --inplace
'''

import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension(
    name="CythonSVR",
    sources=["CythonSVR.pyx", "spherical_volume_rendering_util.cpp"],
    language="c++",
    extra_compile_args=["-std=c++11"],
    include_dirs = [numpy.get_include()],
)]

# Python-3
for e in ext_modules:
    e.cython_directives = {'language_level': "3"}

setup(
    name = 'CythonSVR',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
)
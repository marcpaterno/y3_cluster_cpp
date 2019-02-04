#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, Extension
from Cython.Build import cythonize
import os

# TODO: get these from CMake
include_dirs = [os.path.join(os.path.dirname(__file__), '..'),
                '/cosmosis/cubacpp/',
                '/cosmosis/cubacpp/cubacpp/']

setup(
    name='cluster_abundance',
    ext_modules=cythonize(
        Extension('cluster_abundance',
                  sources=['cluster_abundance.pyx'],
                  language='c++',
                  extra_objects=['../liblc_lt.a', '../libinterp.a'],
                  # TODO: get these from CMake
                  extra_compile_args=['-O3', '-Wall', '-std=c++17'],
                  include_dirs=include_dirs,
                  libraries=['cosmosis', 'gsl', 'cuba', 'gslcblas'],
                  library_dirs=['/cosmosis/cosmosis/datablock/'],
                  ),
        annotate=True
        )
    )

#!/usr/bin/env python

"""
setup.py file for aum
"""

from distutils.core import setup, Extension

cosmology_module = Extension('_cosmology',
                           sources=['src/cosmology.i', 'src/cosmology.cpp', 'src/haloes.cpp', 'src/powerspectrum.cpp', 'src/gauleg.cpp'],
                           swig_opts=["-c++"],
                           libraries=['m','gsl','gslcblas'],
                           )

hod_module = Extension('_hod',
                           sources=['src/cosmology.i', 'src/cosmology.cpp',
                               'src/haloes.cpp', 'src/powerspectrum.cpp',
                               'src/gauleg.cpp', 'src/hod.i', 'src/hod.cpp'],
                           swig_opts=["-c++"],
                           libraries=['m','gsl','gslcblas'],
                           )

setup (name        = 'aum',
       version     = '1.0rc',
       author      = "Surhud More",
       url         = "http://member.ipmu.jp/surhud.more/research",
       author_email= "surhud.more@ipmu.jp",
       description = """A Unified Modelling scheme for galaxy abundance, galaxy clustering and galaxy-galaxy lensing""",
       ext_modules = [cosmology_module, hod_module],
       license     = ['GPL'],
       py_modules  = ["cosmology", "hod"],
       package_dir = { '':'src'},
       )


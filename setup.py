#!/usr/bin/env python

"""
setup.py file for aum
"""
from setuptools import setup, Extension
from setuptools.command.build_py import build_py
from distutils.command.build import build

cosmology_module = Extension(name='_cosmology',
                           sources=['src/cosmology.i', 'src/cosmology.cpp', 'src/haloes.cpp', 'src/powerspectrum.cpp', 'src/gauleg.cpp'],
                           swig_opts=["-c++"],
                           libraries=['m','gsl','gslcblas'],
                           )

hod_module = Extension(name='_hod',
                           sources=['src/cosmology.i', 'src/cosmology.cpp',
                               'src/haloes.cpp', 'src/powerspectrum.cpp',
                               'src/gauleg.cpp', 'src/hod.i', 'src/hod.cpp'],
                           swig_opts=["-c++"],
                           libraries=['m','gsl','gslcblas'],
                           )

# Build extensions before python modules,
# or the generated SWIG python files will be missing.
class BuildPy(build_py):
    def run(self):
        self.run_command('build_ext')
        super(build_py, self).run()

setup (name        = 'aum',
       version     = '1.0rc',
       author      = "Surhud More",
       url         = "http://member.ipmu.jp/surhud.more/research",
       author_email= "surhud.more@ipmu.jp",
       description = """A Unified Modelling scheme for galaxy abundance, galaxy clustering and galaxy-galaxy lensing""",
       ext_modules = [cosmology_module, hod_module],
       py_modules  = ["cosmology", "hod"],
       package_dir = { '':'src'},
       cmdclass={'build_py': BuildPy,},
       )


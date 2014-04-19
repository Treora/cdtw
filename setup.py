from distutils.core import setup, Extension

module1 = Extension('cdtw', sources = ['cdtwmodule.c'])

setup (name = 'cDTW',
    version = '1.0',
    description = 'Constrained Dynamic Time Warping (cDTW) implementation',
    ext_modules = [module1])

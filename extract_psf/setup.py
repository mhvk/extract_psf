# File setup.py

# python setup.py build
# works but leaves in build directory

# Fortran source code
files = ['fortran/extract_psf.f95',]

## scypy_distutils Script
from numpy.distutils.core import setup, Extension

## setup the python module
setup(name="extrpsf",  # name of the package to import later
      ## Build fortran wrappers, uses f2py
      ## directories to search for libraries defined in setup.cfg
      ext_modules=[Extension('extrpsf',
                             files,
                             # libraries=[],
                             # library_dirs=[],
                             # include_dirs=[],
                             f2py_options=['only:'] + ['extract_psf',
                                                       'psfmany'] + [':'])],
      ## Install these to their own directory
      # package_dir={'sf':'Lib'},
      # packages=["sf"]
      )

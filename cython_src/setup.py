#/Users/chplatt/Library/Enthought/Canopy_64bit/User/bin/python setup.py build_ext --inplace
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
 
setup(
  name = 'Demos',
  ext_modules=[ 
    Extension("ham", 
              sources=["ham.pyx", "hamiltonian.cpp", "symmetry.cpp"], # Note, you can link against a c++ library instead of including the source
              language="c++",
              extra_link_args=["-lgsl", "-lgslcblas", "-lCGAL"]),
    ],
  cmdclass = {'build_ext': build_ext},
 
)
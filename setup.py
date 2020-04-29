from distutils.core import setup, Extension
import numpy

# define the extension module for the TDMA algorithm 
tdma_module = Extension('tdma',
                    include_dirs=[numpy.get_include()],
                    sources = ['msynchro/tdma.c'])


# run the setup
setup(name = 'msynchro',
	  version = '1.0',
	  packages = ["msynchro"],
	  description = 'Functions for synchrotron radiation and particle evolution',
	  author_email = 'matthews@ast.cam.ac.uk',
	  author = 'James Matthews',
	  ext_modules=[tdma_module], 
	  py_modules=["msynchro"])  

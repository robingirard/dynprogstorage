#https://packaging.python.org/guides/distributing-packages-using-setuptools/

from setuptools import setup, find_packages
#from distutils.core import setup
#from distutils.core import Extension
#from distutils.extension import Extension
#from setuptools import setup
from Cython.Build import cythonize

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


dynprogstorage_ext = Extension('dynprogstorage.Wrapper_dynprogstorage',
    sources=['dynprogstorage/Wrapper_dynprogstorage.pyx',
             'dynprogstorage/cplfunction.hpp',
             'dynprogstorage/cplfunction.cpp'],
    #libraries=['boost_python38'],
    extra_compile_args=['-std=c++11']
)


#libboost_python38.dylib
setup(
    cmdclass={'build_ext': build_ext},
    name='dynprogstorage',
    version='0.1.4',
    ext_modules=cythonize([dynprogstorage_ext]),
    #long_description=readme,
    #long_description_content_type='text/markdown',
    author='Robin Girard',
    author_email='robin.girard@mines-paristech.fr',license=license,
    packages=find_packages()
)



# from setuptools import setup, find_packages
# # from setuptools.command.
# import sys
# import os
# from Cython.Build import cythonize
# from Cython.Distutils import build_ext
# from setuptools import Extension
# import numpy as np


# extensions = [Extension("progdynstorageModule.src.mod_cplfunction",
#                         sources=['progdynstorageModule/src/dynprogstorage.cpp'],
#                         include_dirs=[np.get_include(),"../inst/include/"],
#                         libraries=['boost_python38'])]

# compiler_check_directives = ['boundscheck', 'initializedcheck', 'nonecheck']
# compiler_directives = {'language_level': 3}

# setup(
#     name='progdynstorageModule',
#     version='0.1',
#     author='',
#     author_email='',
#     packages=find_packages(),
#     url='',
#     license='See LICENSE.txt',
#     description="A dynamic optimisation program for storage optimisation",
#     extras_require={
#     },
#     classifiers=[
#         "Development Status ::  Alpha",
#     ],
#     cmdclass={'build_ext': build_ext
#               },
#     ext_modules=cythonize(extensions, annotate=True, compiler_directives=compiler_directives),
#     zip_safe=False
# )

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
    sources=['dynprogstorage/Wrapper_dynprogstorage.pyx'],
    #libraries=['boost_python38'],
    #extra_compile_args=['-std=c++11']
)

#set INCLUDE=C:\Program Files (x86)\Windows Kits\10\Include\10.0.17763.0\ucrt
#set INCLUDE=C:\Program Files (x86)\Windows Kits\10\Include\10.0.17763.0\ucrt;C:\Program Files (x86)\Windows Kits\10\Include\10.0.17763.0\shared;C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\lib\amd64;C:\Program Files (x86)\Windows Kits\10\bin\10.0.17763.0\x64
#set LIB=C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\lib\amd64;C:\Program Files (x86)\Windows Kits\10\Lib\10.0.17763.0\um\x64;C:\Program Files (x86)\Windows Kits\10\Lib\10.0.17763.0\ucrt\x64;
#libboost_python38.dylib
#C:\Users\robin.girard\AppData\Local\Programs\Python\Python38-32\Scripts\;C:\Users\robin.girard\AppData\Local\Programs\Python\Python38-32\;%PyCharm%;C:\Program Files (x86)\Windows Kits\10\Include
#set INCLUDE=C:\Program Files (x86)\Windows Kits\10\Include\10.0.18362.0\ucrt;C:\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.27.29110\include;C:\Program Files (x86)\Windows Kits\10\Include\10.0.18362.0\shared;C:\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.27.29110\lib\onecore\x64;C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\lib\amd64;C:\Program Files (x86)\Windows Kits\10\bin\10.0.18362.0\x64

#set LIB=C:\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.27.29110\lib\onecore\x64;C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\lib\amd64;C:\Program Files (x86)\Windows Kits\10\Lib\10.0.18362.0\um\x64;C:\Program Files (x86)\Windows Kits\10\Lib\10.0.18362.0\um\x64;C:\Program Files (x86)\Windows Kits\10\Lib\10.0.18362.0\ucrt\x64
setup(
    cmdclass={'build_ext': build_ext},
    name='dynprogstorage',
    version='0.1.1',
    ext_modules=cythonize([dynprogstorage_ext]),
    #long_description=readme,
    #long_description_content_type='text/markdown',
    author='Robin Girard',
    author_email='robin.girard@mines-paristech.fr',license=license,
    packages=find_packages(exclude=('tests', 'docs'))
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

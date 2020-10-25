#https://packaging.python.org/guides/distributing-packages-using-setuptools/

import os
from setuptools import setup, find_packages
#from distutils.core import setup
#from distutils.core import Extension
#from distutils.extension import Extension
#from setuptools import setup
try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    cythonize = None



from distutils.core import setup
from distutils.extension import Extension



# https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#distributing-cython-modules
def no_cythonize(extensions, **_ignore):
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in (".pyx", ".py"):
                if extension.language == "c++":
                    ext = ".cpp"
                else:
                    ext = ".c"
                sfile = path + ext
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions


indlude_dirs=['']

#indlude_dirs=['C:\\Program Files (x86)\\Windows Kits\\10\\Include\\10.0.18362.0\\ucrt',
              #'C:\\Microsoft Visual Studio\\2019\\Community\\VC\\Tools\\MSVC\\14.27.29110\\include',
#              'C:\\Program Files (x86)\\Windows Kits\\10\\Include\\10.0.18362.0\\shared',
              # 'C:\\Microsoft Visual Studio\\2019\\Community\\VC\\Tools\\MSVC\\14.27.29110\\lib\\onecore\\x64',
#              'C:\\Program Files (x86)\\Microsoft Visual Studio 10.0\\VC\\lib\\amd64']
              # 'C:\\Program Files (x86)\\Windows Kits\\10\\bin\\10.0.18362.0\\x64']

dynprogstorage_ext = [Extension('dynprogstorage.Wrapper_dynprogstorage',
    sources=['dynprogstorage/Wrapper_dynprogstorage.pyx',
            #'dynprogstorage/Wrapper_dynprogstorage.cpp',
            #'dynprogstorage/cplfunction.cpp',
            #'dynprogstorage/cplfunction.hpp',
             ],

#    include_dirs=indlude_dirs,  # put include paths here
#    library_dirs=['C:\\Program Files (x86)\\Microsoft Visual Studio 10.0\\VC\\lib\\amd64',
#                  'C:\\Microsoft Visual Studio\\2019\\Community\\VC\\Tools\\MSVC\\14.27.29110\\lib\\onecore\\x64',
                  #,
#                  'C:\\Program Files (x86)\\Windows Kits\\10\\Lib\\10.0.18362.0\\um\\x64',
#                  'C:\\Program Files (x86)\\Windows Kits\\10\\Lib\\10.0.18362.0\\ucrt\\x64'],  # usually need your Windows SDK stuff here
    language='c++',
    #libraries=['boost_python38'],
    #extra_compile_args=['-std=c++11']
)]

#extensions = [
#    Extension("cypack.utils", ["src/cypack/utils.pyx"]),
#    Extension("cypack.answer", ["src/cypack/answer.pyx"]),
#    Extension("cypack.fibonacci", ["src/cypack/fibonacci.pyx"]),
#    Extension(
#        "cypack.sub.wrong",
#        ["src/cypack/sub/wrong.pyx", "src/cypack/sub/helper.c"]
#    ),
#]

CYTHONIZE = bool(int(os.getenv("CYTHONIZE", 0))) and cythonize is not None

if CYTHONIZE:
    compiler_directives = {"language_level": 3, "embedsignature": True}
    extensions = cythonize(dynprogstorage_ext, compiler_directives=compiler_directives)
else:
    extensions = no_cythonize(dynprogstorage_ext)

#set INCLUDE=C:\Program Files (x86)\Windows Kits\10\Include\10.0.17763.0\ucrt
#set INCLUDE=C:\Program Files (x86)\Windows Kits\10\Include\10.0.17763.0\ucrt;C:\Program Files (x86)\Windows Kits\10\Include\10.0.17763.0\shared;C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\lib\amd64;C:\Program Files (x86)\Windows Kits\10\bin\10.0.17763.0\x64
#set LIB=C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\lib\amd64;C:\Program Files (x86)\Windows Kits\10\Lib\10.0.17763.0\um\x64;C:\Program Files (x86)\Windows Kits\10\Lib\10.0.17763.0\ucrt\x64;
#libboost_python38.dylib
#C:\Users\robin.girard\AppData\Local\Programs\Python\Python38-32\Scripts\;C:\Users\robin.girard\AppData\Local\Programs\Python\Python38-32\;%PyCharm%;C:\Program Files (x86)\Windows Kits\10\Include
#set INCLUDE=C:\Program Files (x86)\Windows Kits\10\Include\10.0.18362.0\ucrt;C:\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.27.29110\include;C:\Program Files (x86)\Windows Kits\10\Include\10.0.18362.0\shared;C:\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.27.29110\lib\onecore\x64;C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\lib\amd64;C:\Program Files (x86)\Windows Kits\10\bin\10.0.18362.0\x64

#set LIB=C:\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.27.29110\lib\onecore\x64;C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\lib\amd64;C:\Program Files (x86)\Windows Kits\10\Lib\10.0.18362.0\um\x64;C:\Program Files (x86)\Windows Kits\10\Lib\10.0.18362.0\um\x64;C:\Program Files (x86)\Windows Kits\10\Lib\10.0.18362.0\ucrt\x64
setup(
    cmdclass={'build_ext': build_ext},
 #   setup_args={'script_args': ["--compiler=msvc"]},
    name='dynprogstorage',
    version='0.1.10',
    ext_modules=extensions,
    #long_description=readme,
    #long_description_content_type='text/markdown',
    author='Robin Girard',
    author_email='robin.girard@mines-paristech.fr',license=license,
    include_dirs = indlude_dirs,
    packages=find_packages(),

)


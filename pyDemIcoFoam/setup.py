long_description = """
Python wrapper for demIcoFoam.

demIcoFoam is an OpenFOAM CFD solver which can interact with a
discrete element model.
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os

assert not os.getenv("FOAM_SRC") is None, "Cannot find OpenFOAM install. Did you source ~/OpenFOAM/OpenFOAM-v3.0+/etc/bashrc ?"

assert os.path.isfile(os.getenv("FOAM_APPBIN")+"/icoFoam"), "Cannot find OpenFOAM binaries. Did you build OpenFOAM in the default location? "


ext = [
    Extension("_pyDemIcoFoam",
              sources=["pyDemIcoFoam.pyx",
                       "demIcoFoam.C"],
              include_dirs = [
                  os.getenv("FOAM_SRC")+"/finiteVolume/lnInclude",
                  os.getenv("FOAM_SRC")+"/meshTools/lnInclude",
                  os.getenv("FOAM_SRC")+"/OpenFOAM/lnInclude",
                  os.getenv("FOAM_SRC")+"/OSspecific/POSIX/lnInclude"],
              extra_compile_args= ["-m64", "-Dlinux64", "-DWM_ARCH_OPTION=64",
                                   "-DWM_DP", "-DWM_LABEL_SIZE=32", "-Wall",
                                   "-Wextra", "-Wold-style-cast",
                                   "-Wnon-virtual-dtor",
                                   "-Wno-unused-parameter",
                                   "-Wno-invalid-offsetof", "-O3",
                                   "-DNoRepository", "-ftemplate-depth-100"],
              extra_link_args=["-Xlinker", "--add-needed",
                                  "-Xlinker", "--no-as-needed"],
              libraries = ["finiteVolume", "meshTools", "OpenFOAM", "dl", "m"],
              library_dirs = [os.getenv("FOAM_LIBBIN")],
              language="c++",             # generate C++ code
)]

setup(
    name = 'pyDemIcoFoam',
    packages = ["pyDemIcoFoam"], # this must be the same as the name above
    version = __import__('itasca').__version__,
    description = "Python wrapper for demIcoFoam.",
    long_description = long_description,
    author = 'Jason Furtney',
    requires = ['numpy'],
    author_email = 'jkfurtney@gmail.com',
    url = "https://github.com/jkfurtney/PFC3D_OpenFOAM",
    keywords = 'OpenFOAM,CFD,icoFoam,PFC3D,PFC,DEM'.split(","),
    license          = "BSD",
    classifiers = [
        'Programming Language :: Python :: 2',
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: BSD License",
        'Topic :: Scientific/Engineering :: Interface Engine/Protocol Translator',
        "Intended Audience :: Science/Research"
    ],
    ext_modules = cythonize(ext, language="c++"))

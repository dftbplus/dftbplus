#!/usr/bin/env python
from distutils.core import setup

setup(
    name="dptools",
    version="0.1p3",
    description="Tools to process DFTB+ related data",
    author="B. Aradi, B. Hourahine, A. Pecchia",
    url="http://www.dftb-plus.info",
    platforms="platform independent",
    package_dir={ "": "src" },
    packages=[ "dptools", ],
    scripts=[
        "bin/dp_bands",
        "bin/dp_dos",
        "bin/gen2xyz",
        "bin/gen2cif",
        "bin/xyz2gen",
    ],
    classifiers=[
        "Programming Language :: Python",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
    long_description="""
Processing and converting data related to the DFTB+ package
-----------------------------------------------------------
A few scripts which should make the life of DFTB+ users easier, by providing
functions to process and convert various DFTB+ data formats.
""",
    requires=[ "numpy" ]
)

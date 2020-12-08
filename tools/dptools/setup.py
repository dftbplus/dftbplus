#!/usr/bin/env python3
from distutils.core import setup

setup(
    name="dptools",
    version='20.2.1',
    description="Tools to process DFTB+ related data",
    author="DFTB+ developers",
    url="http://www.dftbplus.org",
    platforms="platform independent",
    package_dir={"": "src"},
    packages=["dptools", "dptools.scripts"],
    scripts=[
        "bin/dp_bands",
        "bin/dp_dos",
        "bin/gen2cif",
        "bin/gen2xyz",
        "bin/makecube",
        "bin/repeatgen",
        "bin/xyz2gen",
        "bin/straingen",
    ],
    classifiers=[
        "Programming Language :: Python",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: LGPL License",
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

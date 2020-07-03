#!/usr/bin/env python3
from distutils.core import setup

setup(
    name='pythonapi',
    version='0.1',
    description='Python interface to DFTB+',
    author='DFTB+ developers',
    url='http://www.dftbplus.org',
    platforms='platform independent',
    package_dir={'': 'src'},
    packages=[''],
    classifiers=[
        'Programming Language :: Python',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: LGPL',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering',
    ],
    long_description='''
        A ctypes based Python interface for DFTB+
        -----------------------------------------
        Interface module for the communication between DFTB+ and
        Python via the foreign function C-library ctypes. Provides
        methods for initializing and configuring a DFTB+ calculator.
        ''',
    requires=['numpy']
)

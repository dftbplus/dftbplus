#!/usr/bin/env python3

'''
Legacy setup.py file that gathers its
configuration from setup.cfg and pyproject.toml
'''

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if __name__ == '__main__':
    setup()

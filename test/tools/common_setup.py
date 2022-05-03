#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Common items for the dptools and ptools tests frameworks.'''

import unittest
import shutil
import os.path
import tempfile


class TestWithWorkDir(unittest.TestCase):
    '''Base class for tests which need work directory for file output.

    You should create test cases which need a work directory (as they write
    during tests) by deriving this class. In the respective setUp() method you
    have to set self.workroot and self.inputdir to the appropriate values and
    invoke the setUp() of this class after that.

    Each test will create a temporary directory within self.workroot. The name
    of the directory is prefixed with the name of the test class and the invoked
    test method.
    '''

    def setUp(self):
        '''Note: self.inputdir and self.workroot must have been already set!'''
        prefix = '.'.join(self.id().split('.')[1:]) + '_'
        self.workdir = tempfile.mkdtemp(prefix=prefix, dir=self.workroot)
        self.keepworkdir = False

    def tearDown(self):
        if not self.keepworkdir:
            shutil.rmtree(self.workdir)

    def get_input(self, fname):
        '''Returns file name prefixed with the input directory.

        Args:
            fname (str): File name to prefix with input directory.
        '''
        return os.path.join(self.inputdir, fname)

    def get_output(self, fname):
        '''Returns file name prefixed with the output directory.

        Args:
            fname (str): File name to prefix with output directory.
        '''
        return os.path.join(self.workdir, fname)

#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Common items for the dptools tests framework.'''

import sys
import unittest
import tempfile
import shutil
import os.path
if sys.version_info[0] >= 3:
    from io import StringIO
else:
    from cStringIO import StringIO
from dptools.gen import Gen


def gen_file_equals(current, reference):
    '''Checks whether a gen file equals to an other one.

    Args:
        current (str): Name of gen file to check.
        reference (str): Name of reference gen file.
    '''
    curgen = Gen.fromfile(current)
    refgen = Gen.fromfile(reference)
    return curgen.equals(refgen)


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


class OutputCatcher:
    '''Captures the stdandard output of a routine'''

    def __enter__(self):
        self._output = None
        self._stdout = sys.stdout
        self._stringio = StringIO()
        sys.stdout = self._stringio
        return self

    def __exit__(self, *args):
        sys.stdout = self._stdout
        self._output = self._stringio.getvalue()
        self._stringio.close()

    def get(self):
        '''Returns the catched output as string.'''
        return self._output

    def get_as_stringio(self):
        '''Returns the catched output as StringIO file object.'''
        return StringIO(self._output)

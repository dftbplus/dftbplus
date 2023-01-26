#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Tests for resultstag.'''

import unittest
import os.path
import sys
from dftbplus_ptools.resultstag import Output
import common
from resultstag.resultstag_ref import real_ref
from resultstag.resultstag_ref import integer_ref
from resultstag.resultstag_ref import complex_ref
from resultstag.resultstag_ref import logical_ref
from resultstag.resultstag_ref import scalar_ref


SCRIPTDIR = os.path.dirname(sys.argv[0])


class ResulstagTest(common.TestWithWorkDir):
    '''General tests for reading results.tag'''


    def setUp(self):
        self.inputdir = os.path.join(SCRIPTDIR, "resultstag")
        self.workroot = './'
        common.TestWithWorkDir.setUp(self)


    def test_real(self):
        """tests reading of real values"""
        self.assertTrue(self.common(real_ref, "real.tag"))


    def test_integer(self):
        """tests reading of integer values"""
        self.assertTrue(self.common(integer_ref, "integer.tag"))


    def test_complex(self):
        """tests reading of complex values"""
        self.assertTrue(self.common(complex_ref, "complex.tag"))


    def test_logical(self):
        """tests reading of logical values"""
        self.assertTrue(self.common(logical_ref, "logical.tag"))


    def test_scalar(self):
        """tests reading of scalar values"""
        self.assertTrue(self.common(scalar_ref, "scalar.tag"))


    def test_none(self):
        """tests the case if tag not present"""
        path = os.path.join(self.inputdir, "scalar.tag")
        results = Output(filename=path)
        value = results.get_eigenvalues()
        self.assertTrue(True if value is None else False)


    def common(self, reference, filename, rtol=1e-05, atol=1e-08):
        """common part of tests"""
        path = os.path.join(self.inputdir, filename)
        results = Output(filename=path)
        value = results.get_cell_volume()
        return common.type_diff(reference, value, rtol=rtol, atol=atol)


if __name__ == '__main__':
    unittest.main()

#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Tests for eigenvecout.'''


import unittest
import os.path
import sys
from dftbplus_ptools.eigenvecout import Eigenvecout
import common
from eigenvecout.eigenvec_ref import eigenvec1_ref
from eigenvecout.eigenvec_ref import eigenvec2_ref
from eigenvecout.eigenvec_ref import eigenvec3_ref


SCRIPTDIR = os.path.dirname(sys.argv[0])


class EigenvecoutTest(common.TestWithWorkDir):
    '''General tests for eigenvecout'''


    def setUp(self):
        self.inputdir = os.path.join(SCRIPTDIR, "eigenvecout")
        self.workroot = './'
        common.TestWithWorkDir.setUp(self)


    def test_eigenvec_1(self):
        """test for pattern 2,4,6"""
        path = os.path.join(self.inputdir, "eigenvec1.out")
        eigenvecclass = Eigenvecout(filename=path)
        self.assertTrue(common.type_diff(eigenvecclass.read_eigenvec(),
                                         eigenvec1_ref))


    def test_eigenvec_2(self):
        """test for pattern 1,3,5"""
        path = os.path.join(self.inputdir, "eigenvec2.out")
        eigenvecclass = Eigenvecout(filename=path)
        self.assertTrue(common.type_diff(eigenvecclass.read_eigenvec(),
                                         eigenvec2_ref))


    def test_eigenvec_3(self):
        """test for pattern 7,8,9"""
        path = os.path.join(self.inputdir, "eigenvec3.out")
        eigenvecclass = Eigenvecout(filename=path)
        self.assertTrue(common.type_diff(eigenvecclass.read_eigenvec(),
                                         eigenvec3_ref))


if __name__ == '__main__':
    unittest.main()

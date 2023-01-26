#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Tests for detailedout.'''

import unittest
import os.path
import sys
from dftbplus_ptools.detailedout import Detailedout
import common


SCRIPTDIR = os.path.dirname(sys.argv[0])


class DetailoutTest(common.TestWithWorkDir):
    '''General tests for detailedout'''

    def setUp(self):
        self.inputdir = os.path.join(SCRIPTDIR, "detailedout")
        self.workroot = './'
        common.TestWithWorkDir.setUp(self)

    def test_scc_convergence_true(self):
        """test for scc convergence (converged)"""
        path = os.path.join(self.inputdir, "scc.out")
        sccclass = Detailedout(filename=path)
        self.assertTrue(sccclass.check_for_scc_convergence())


    def test_scc_convergence_false(self):
        """test for scc convergence (not converged)"""
        path = os.path.join(self.inputdir, "geometry.out")
        sccclass = Detailedout(filename=path)
        self.assertFalse(sccclass.check_for_scc_convergence())


    def test_scc_convergence_fail(self):
        """test for scc convergence (file not existing)"""
        sccclass = Detailedout(filename=self.inputdir)
        self.assertTrue(True if sccclass.check_for_scc_convergence() is None
                        else False)


    def test_geometry_convergence_true(self):
        """test for geometry convergence (converged)"""
        path = os.path.join(self.inputdir, "geometry.out")
        geometryclass = Detailedout(filename=path)
        self.assertTrue(geometryclass.check_for_geometry_convergence())


    def test_geometry_convergence_false(self):
        """test for geometry convergence (not converged)"""
        path = os.path.join(self.inputdir, "scc.out")
        geometryclass = Detailedout(filename=path)
        self.assertFalse(geometryclass.check_for_geometry_convergence())


    def test_geometry_convergence_fail(self):
        """test for geometry convergence (file not existing)"""
        geometryclass = Detailedout(filename=self.inputdir)
        self.assertTrue(True if geometryclass.check_for_geometry_convergence()
                        is None else False)


if __name__ == '__main__':
    unittest.main()

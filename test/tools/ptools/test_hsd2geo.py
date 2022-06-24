#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Tests for hsd to geometry conversion.'''

import unittest
import os.path
import sys
from dftbplus_ptools.gen import Gen
from dftbplus_ptools.xyz import Xyz
import common


SCRIPTDIR = os.path.dirname(sys.argv[0])


class Hsd2geoTest(common.TestWithWorkDir):
    '''General tests for conversion from hsd'''


    def setUp(self):
        self.inputdir = os.path.join(SCRIPTDIR, "hsd2geo")
        self.workroot = './'
        common.TestWithWorkDir.setUp(self)


    def test_xyz(self):
        """test for conversion from hsd (xyz format)"""
        self.assertTrue(self.common(Xyz, "xyz.ref", "xyz.hsd"))

    def test_gen(self):
        """test for conversion from hsd (gen format)"""
        self.assertTrue(self.common(Gen, "gen.ref", "gen.hsd"))


    def test_explicit(self):
        """test for conversion from hsd (explicit format)"""
        self.assertTrue(self.common(Gen, "explicit.ref",
                                    "explicit.hsd", tolerance=2.6E-10))

    def test_explicit_rel(self):
        """test for conversion from hsd (explicit format, relative coords)"""
        self.assertTrue(self.common(Gen, "explicit_rel.ref",
                                    "explicit_rel.hsd"))

    def test_explicit_angstrom(self):
        """test for conversion from hsd (explicit format, unit = angstrom)"""
        self.assertTrue(self.common(Gen, "explicit.ref",
                                    "explicit_angstrom.hsd"))

    def test_explicit_pm(self):
        """test for conversion from hsd (explicit format, unit = pm)"""
        self.assertTrue(self.common(Gen, "explicit.ref",
                                    "explicit_pm.hsd"))


    def common(self, class_obj, reference, filename, tolerance=1E-10):
        """common part of tests"""
        path = os.path.join(self.inputdir, filename)
        value = class_obj.fromhsd(filename=path)
        reference_value = class_obj.fromfile(os.path.join(self.inputdir,
                                                          reference))
        return reference_value.equals(value, tolerance=tolerance)


if __name__ == '__main__':
    unittest.main()

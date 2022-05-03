#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Tests for hsdinput.'''

import unittest
import os.path
import sys
import hsd
from dftbplus_ptools.hsdinput import Changehsd
from dftbplus_ptools.gen import Gen
import common


SCRIPTDIR = os.path.dirname(sys.argv[0])


class HsdinputTest(common.TestWithWorkDir):
    '''General tests for hsd input manipulation'''


    def setUp(self):
        self.inputdir = os.path.join(SCRIPTDIR, "hsdinput")
        self.workroot = './'
        common.TestWithWorkDir.setUp(self)


    def test_skdir_str(self):
        """test for setting skdir using a string"""
        hsd_dict = Changehsd()
        hsd_dict.set_skdir("hsdinput")
        reference = Changehsd(directory=self.inputdir,
                              filename="skdir_str_ref.hsd")
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_skdir_dict(self):
        """test for setting skdir using a dictionary"""
        hsd_dict = Changehsd()
        skdir_dict = {'O-O': 'hsdinput/O-O.skf', "O-H": "hsdinput/O-H.skf",
                      'H-O': 'hsdinput/H-O.skf', "H-H": "hsdinput/H-H.skf"}
        hsd_dict.set_skdir(skdir_dict)
        reference = Changehsd(directory=self.inputdir,
                              filename="skdir_dict_ref.hsd")
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_kpts_tuple(self):
        """tests kpts given as tuple"""
        hsd_dict = Changehsd()
        hsd_dict.set_kpts((2, 2, 3))
        # str is used to not needed tempfiles and to get consistent data types
        hsd_str = hsd.dump_string(hsd_dict.get_hsd())
        hsd_dict = hsd.load_string(hsd_str)
        reference = Changehsd(directory=self.inputdir,
                              filename="kpts_tuple_ref.hsd")
        self.assertTrue(common.type_diff(reference.get_hsd(), hsd_dict))


    def test_kpts_list(self):
        """tests kpts given as list"""
        hsd_dict = Changehsd()
        hsd_dict.set_kpts([[0, 0, 0, 0.25], [0.25, 0.25, 0.25, 0.75]])
        # str is used to not needed tempfiles and to get consistent data types
        hsd_str = hsd.dump_string(hsd_dict.get_hsd())
        hsd_dict = hsd.load_string(hsd_str)
        reference = Changehsd(directory=self.inputdir,
                              filename="kpts_list_ref.hsd")
        self.assertTrue(common.type_diff(reference.get_hsd(), hsd_dict))


    def test_filling_fermi(self):
        """test for setting Fermi-filling"""
        hsd_dict = Changehsd()
        hsd_dict.set_filling("Fermi", order=5, temperature=10)
        reference = Changehsd(directory=self.inputdir,
                              filename="filling_fermi_ref.hsd")
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_filling_methfesselpaxton(self):
        """test for setting MethfesselPaxton-filling"""
        hsd_dict = Changehsd()
        hsd_dict.set_filling("MethfesselPaxton", order=5, temperature=10)
        reference = Changehsd(directory=self.inputdir,
                              filename="filling_methfesselpaxton_ref.hsd")
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_maxang_dict(self):
        """test for setting max angular momenta via dictinary"""
        hsd_dict = Changehsd()
        maxangs = {"O": "p", 'H': 's'}
        hsd_dict.set_maxang(maxangs=maxangs)
        reference = Changehsd(directory=self.inputdir,
                              filename="maxang_ref.hsd")
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_maxang_try_reading_dict(self):
        """test for reading max angular momenta with skdir as dictinary"""
        hsd_dict = Changehsd()
        skdir_dict = {'O-O': f'{self.inputdir}/O-O.skf',
                      "O-H": f"{self.inputdir}/O-H.skf",
                      'H-O': f'{self.inputdir}/H-O.skf',
                      "H-H": f"{self.inputdir}/H-H.skf"}
        hsd_dict.set_skdir(skdir_dict)
        hsd_dict.set_maxang(try_reading=["O", "H"])
        reference = Changehsd(directory=self.inputdir,
                              filename="maxang_try_reading_ref.hsd")
        reference.set_skdir(skdir_dict)
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_maxang_try_reading_str(self):
        """test for reading max angular momenta with skdir as string"""
        hsd_dict = Changehsd()
        hsd_dict.set_skdir(self.inputdir)
        hsd_dict.set_maxang(try_reading=["O", "H"])
        reference = Changehsd(directory=self.inputdir,
                              filename="maxang_try_reading_ref.hsd")
        reference.set_skdir(self.inputdir)
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_maxang_error(self):
        """tests raised error if species is in maxangs and try_reading"""
        hsd_dict = Changehsd()
        with self.assertRaises(ValueError):
            hsd_dict.set_maxang(maxangs={"O": "p"}, try_reading=["O"])


    def test_geometry(self):
        """test for setting geometry"""
        gen = Gen.fromhsd("geometry_ref.hsd", directory=self.inputdir)
        geo = gen.geometry
        hsd_dict = Changehsd()
        hsd_dict.set_geometry(geo)
        # str is used to not needed tempfiles and to get consistent data types
        hsd_str = hsd.dump_string(hsd_dict.get_hsd())
        hsd_dict = hsd.load_string(hsd_str)
        reference = hsd.load(os.path.join(self.inputdir, "geometry_ref.hsd"))
        self.assertTrue(common.type_diff(reference, hsd_dict))


    def test_get_basic_input(self):
        """test for ’get_basic_input’ to test small functions"""
        hsd_dict = Changehsd()
        hsd_dict.get_basic_input(driver="DFTB", drivermaxforce=1,
                                 drivermaxsteps=5, scc=True, scctol=1E-005,
                                 resultstag=True, unprocessed=True,
                                 forces=True)
        reference = Changehsd(directory=self.inputdir,
                              filename="basic_input_ref.hsd")
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


if __name__ == '__main__':
    unittest.main()

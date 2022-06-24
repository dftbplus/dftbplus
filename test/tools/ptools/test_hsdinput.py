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
from dftbplus_ptools.hsdinput import Hsdinput
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
        hsd_dict = Hsdinput(hamiltonian="DFTB")
        hsd_dict.set_skdir("hsdinput")
        path = os.path.join(self.inputdir, "skdir_str_ref.hsd")
        reference = Hsdinput(filename=path)
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_skdir_dict(self):
        """test for setting skdir using a dictionary"""
        hsd_dict = Hsdinput(hamiltonian="DFTB")
        skdir_dict = {'O-O': 'hsdinput/O-O.skf', "O-H": "hsdinput/O-H.skf",
                      'H-O': 'hsdinput/H-O.skf', "H-H": "hsdinput/H-H.skf"}
        hsd_dict.set_skdir(skdir_dict)
        path = os.path.join(self.inputdir, "skdir_dict_ref.hsd")
        reference = Hsdinput(filename=path)
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_kpts_tuple(self):
        """tests kpts given as tuple"""
        hsd_dict = Hsdinput(hamiltonian="DFTB")
        hsd_dict.set_kpts((2, 2, 3))
        path = os.path.join(self.inputdir, "kpts_tuple_ref.hsd")
        reference = Hsdinput(filename=path)
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_kpts_list(self):
        """tests kpts given as list"""
        hsd_dict = Hsdinput(hamiltonian="DFTB")
        hsd_dict.set_kpts([[0, 0, 0, 0.25], [0.25, 0.25, 0.25, 0.75]])
        path = os.path.join(self.inputdir, "kpts_list_ref.hsd")
        reference = Hsdinput(filename=path)
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_filling_fermi(self):
        """test for setting Fermi-filling"""
        hsd_dict = Hsdinput(hamiltonian="DFTB")
        hsd_dict.set_filling("Fermi", order=5, temperature=10)
        path = os.path.join(self.inputdir, "filling_fermi_ref.hsd")
        reference = Hsdinput(filename=path)
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_filling_methfesselpaxton(self):
        """test for setting MethfesselPaxton-filling"""
        hsd_dict = Hsdinput(hamiltonian="DFTB")
        hsd_dict.set_filling("MethfesselPaxton", order=5, temperature=10)
        path = os.path.join(self.inputdir, "filling_methfesselpaxton_ref.hsd")
        reference = Hsdinput(filename=path)
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_maxang_dict(self):
        """test for setting max angular momenta via dictinary"""
        hsd_dict = Hsdinput(hamiltonian="DFTB")
        maxangs = {"O": "p", 'H': 's'}
        hsd_dict.set_maxang(maxangs=maxangs)
        path = os.path.join(self.inputdir, "maxang_ref.hsd")
        reference = Hsdinput(filename=path)
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_maxang_try_reading_dict(self):
        """test for reading max angular momenta with skdir as dictinary"""
        hsd_dict = Hsdinput(hamiltonian="DFTB")
        skdir_dict = {'O-O': f'{self.inputdir}/O-O.skf',
                      "O-H": f"{self.inputdir}/O-H.skf",
                      'H-O': f'{self.inputdir}/H-O.skf',
                      "H-H": f"{self.inputdir}/H-H.skf",
                      "Br-Br": f"{self.inputdir}/Br-Br.skf"}
        hsd_dict.set_skdir(skdir_dict)
        hsd_dict.set_maxang(try_reading=["O", "H", "Br"])
        path = os.path.join(self.inputdir, "maxang_try_reading_ref.hsd")
        reference = Hsdinput(filename=path, hamiltonian="DFTB")
        reference.set_skdir(skdir_dict, )
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_maxang_try_reading_str(self):
        """test for reading max angular momenta with skdir as string"""
        hsd_dict = Hsdinput(hamiltonian="DFTB")
        hsd_dict.set_skdir(self.inputdir)
        hsd_dict.set_maxang(try_reading=["O", "H", "Br"])
        path = os.path.join(self.inputdir, "maxang_try_reading_ref.hsd")
        reference = Hsdinput(filename=path, hamiltonian="DFTB")
        reference.set_skdir(self.inputdir)
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_maxang_error(self):
        """tests raised error if species is in maxangs and try_reading"""
        hsd_dict = Hsdinput()
        with self.assertRaises(ValueError):
            hsd_dict.set_maxang(maxangs={"O": "p"}, try_reading=["O"])


    def test_geometry(self):
        """test for setting geometry"""
        path = os.path.join(self.inputdir, "geometry_ref.hsd" )
        gen = Gen.fromhsd(filename=path)
        geo = gen.geometry
        hsd_dict = Hsdinput()
        hsd_dict.set_geometry(geo)
        reference = hsd.load(os.path.join(self.inputdir, "geometry_ref.hsd"),
                             lower_tag_names=True)
        self.assertTrue(common.type_diff(reference, hsd_dict.get_hsd()))


    def test_get_basic_input(self):
        """test for ’get_basic_input’ to test small functions"""
        hsd_dict = Hsdinput(hamiltonian="DFTB")
        hsd_dict.get_basic_input(driver="DFTB", drivermaxforce=1,
                                 drivermaxsteps=5, scc=True, scctol=1E-005,
                                 resultstag=True, unprocessed=True,
                                 forces=True)
        path = os.path.join(self.inputdir, "basic_input_ref.hsd")
        reference = Hsdinput(filename=path)
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


    def test_geometry_inclusion(self):
        """test for including geometry files"""
        hsd_dict = Hsdinput(dictionary={})
        hsd_dict.set_geometry("geo.gen")
        path = os.path.join(self.inputdir, "geometry_inclusion.ref")
        reference = Hsdinput(filename=path)
        self.assertTrue(common.type_diff(reference.get_hsd(),
                                         hsd_dict.get_hsd()))


if __name__ == '__main__':
    unittest.main()

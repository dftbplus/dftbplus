#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Tests for dp_dos.'''

import sys
import os.path
import unittest
import common
from dptools.scripts.common import ScriptError
import dptools.scripts.dp_dos as dp_dos


SCRIPTDIR = os.path.dirname(sys.argv[0])


class DpdosTest(common.TestWithWorkDir):
    '''General tests for dp_dos'''

    def setUp(self):
        self.inputdir = os.path.join(SCRIPTDIR, 'dp_dos')
        self.workroot = './'
        common.TestWithWorkDir.setUp(self)

    def test_dos(self):
        '''DOS with broadening-function gauss'''
        infile = self.get_input('TiO2_band.out')
        reffile = self.get_input('TiO2.dat')
        outfile = self.get_output('TiO2.dat')
        cmdargs = [infile, outfile]
        dp_dos.main(cmdargs)
        self.assertTrue(common.nxy_file_equals(outfile, reffile))

    def test_dos_grid(self):
        '''DOS with specified grid separation'''
        infile = self.get_input('TiO2_band.out')
        reffile = self.get_input('TiO2_grid-separation.dat')
        outfile = self.get_output('TiO2_grid-separation.dat')
        cmdargs = ['-g', '0.03', infile, outfile]
        dp_dos.main(cmdargs)
        self.assertTrue(common.nxy_file_equals(outfile, reffile))

    def test_dos_fermi(self):
        '''DOS with broadening-function fermi and custom broadening-width'''
        infile = self.get_input('TiO2_band.out')
        reffile = self.get_input('TiO2_fermi-broadening-width0.15.dat')
        outfile = self.get_output('TiO2_fermi-broadening-width0.15.dat')
        cmdargs = ['-f', 'fermi', '-b', '0.15', infile, outfile]
        dp_dos.main(cmdargs)
        self.assertTrue(common.nxy_file_equals(outfile, reffile))

    def test_dos_mporder(self):
        '''DOS with custom order of mp function'''
        infile = self.get_input('TiO2_band.out')
        reffile = self.get_input('TiO2_mporder1.dat')
        outfile = self.get_output('TiO2_mporder1.dat')
        cmdargs = ['-f', 'mp', '-o', '1', infile, outfile]
        dp_dos.main(cmdargs)
        self.assertTrue(common.nxy_file_equals(outfile, reffile))

    def test_dos_sigmarange(self):
        '''DOS with custom broadening width sigma'''
        infile = self.get_input('TiO2_band.out')
        reffile = self.get_input('TiO2_sigma-range7.dat')
        outfile = self.get_output('TiO2_sigma-range7.dat')
        cmdargs = ['-s', '7', infile, outfile]
        dp_dos.main(cmdargs)
        self.assertTrue(common.nxy_file_equals(outfile, reffile))

    def test_dos_spinpolarized(self):
        '''DOS with spin polarization and weighted occupation'''
        infile = self.get_input('band_spin.out')
        reffile = self.get_input('band_spin-polarized_weight-occ.dat')
        outfile = self.get_output('band_spin-polarized_weight-occ.dat')
        cmdargs = ['-w', infile, outfile]
        dp_dos.main(cmdargs)
        self.assertTrue(common.nxy_file_equals(outfile, reffile))

    def test_pdos(self):
        '''PDOS with broadening-function gauss'''
        infile = self.get_input('dos_ti.1.out')
        reffile = self.get_input('band_ti1pdos.dat')
        outfile = self.get_output('band_ti1pdos.dat')
        cmdargs = ['-w', infile, outfile]
        dp_dos.main(cmdargs)
        self.assertTrue(common.nxy_file_equals(outfile, reffile))

    def test_fail_mporder_withoutmp(self):
        '''Failing due to a mporder without specifying broadening type mp.'''
        infile = self.get_input('TiO2_band.out')
        outfile = self.get_output('TiO2.dat')
        cmdargs = ['-o', '1', infile, outfile]
        with self.assertRaises(ScriptError):
            dp_dos.main(cmdargs)

    def test_fail_neg_mporder(self):
        '''Failing due to a negative order of mp function.'''
        infile = self.get_input('TiO2_band.out')
        outfile = self.get_output('TiO2.dat')
        cmdargs = ['-f', 'mp', '-o', '-1', infile, outfile]
        with self.assertRaises(ScriptError):
            dp_dos.main(cmdargs)

    def test_fail_invalid_infile(self):
        '''Failing due to invalid input file.'''
        tempname = common.get_temporary_filename(self.workroot)
        nonexisting_infile = os.path.join(self.workdir, tempname)
        outfile = self.get_output('TiO2.dat')
        cmdargs = [nonexisting_infile, outfile]
        with self.assertRaises(ScriptError):
            dp_dos.main(cmdargs)


if __name__ == '__main__':
    unittest.main()

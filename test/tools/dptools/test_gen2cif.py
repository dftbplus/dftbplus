#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2018  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Tests for gen2cif.'''

import sys
import os.path
import unittest
import tempfile
import common
from dptools.scripts.common import ScriptError
import dptools.scripts.gen2cif as gen2cif


SCRIPTDIR = os.path.dirname(sys.argv[0])


class Gen2cifTest(common.TestWithWorkDir):
    '''General tests for gen2cif'''

    def setUp(self):
        self.inputdir = os.path.join(SCRIPTDIR, 'gen2cif')
        self.workroot = './'
        common.TestWithWorkDir.setUp(self)

    def test_cluster(self):
        '''Absolute coordinates with default cubic cellsize'''
        infile = self.get_input('h2o.234.gen')
        reffile = self.get_input('h2o.234.cif')
        outfile = self.get_output('h2o.234.cif')
        cmdargs = ['-o', outfile, infile]
        gen2cif.main(cmdargs)
        self.assertTrue(common.cif_file_equals(outfile, reffile))

    def test_customcell(self):
        '''Absolute coordinates with custom cubic cellsize'''
        infile = self.get_input('h2o.234.gen')
        reffile = self.get_input('h2o.234-customcell.cif')
        outfile = self.get_output('h2o.234-customcell.cif')
        cmdargs = ['-o', outfile, '-c', '10', infile]
        gen2cif.main(cmdargs)
        self.assertTrue(common.cif_file_equals(outfile, reffile))

    def test_periodic(self):
        '''Absolute coordinates with periodic geometry'''
        infile = self.get_input('h2o.gen')
        reffile = self.get_input('h2o.cif')
        outfile = self.get_output('h2o.cif')
        cmdargs = ['-o', outfile, infile]
        gen2cif.main(cmdargs)
        self.assertTrue(common.cif_file_equals(outfile, reffile))

    def test_fraccoords(self):
        '''Fractional coordinates'''
        infile = self.get_input('gaas-frac.gen')
        reffile = self.get_input('gaas-frac.cif')
        outfile = self.get_output('gaas-frac.cif')
        cmdargs = ['-o', outfile, infile]
        gen2cif.main(cmdargs)
        self.assertTrue(common.cif_file_equals(outfile, reffile))

    def test_fail_superfluous_arguments(self):
        '''Failing due to superfluous arguments.'''
        infile = self.get_input('h2o.234.gen')
        outfile = self.get_output('h2o.234.cif')
        cmdargs = ['-o', outfile, infile, 'something']
        with self.assertRaises(ScriptError):
            gen2cif.main(cmdargs)

    def test_fail_missing_arguments(self):
        '''Failing due to missing arguments.'''
        infile = self.get_input('h2o.234.gen')
        cmdargs = ['-o', infile]
        with self.assertRaises(ScriptError):
            gen2cif.main(cmdargs)

    def test_fail_invalid_infile(self):
        '''Failing due to invalid input file.'''
        temp_file = tempfile.NamedTemporaryFile(dir=self.workroot)
        tempname = temp_file.name
        temp_file.close()
        nonexisting_infile = os.path.join(self.workdir, tempname)
        outfile = self.get_output('h2o.234.cif')
        cmdargs = ['-o', outfile, nonexisting_infile]
        with self.assertRaises(ScriptError):
            gen2cif.main(cmdargs)


if __name__ == '__main__':
    unittest.main()

#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Tests for xyz2gen.'''

import sys
import os.path
import unittest
import common
from dptools.scripts.common import ScriptError
import dptools.scripts.xyz2gen as xyz2gen


SCRIPTDIR = os.path.dirname(sys.argv[0])


class Xyz2genTest(common.TestWithWorkDir):
    '''General tests for xyz2gen'''

    def setUp(self):
        self.inputdir = os.path.join(SCRIPTDIR, 'xyz2gen')
        self.workroot = './'
        common.TestWithWorkDir.setUp(self)

    def test_cluster(self):
        '''Absolute coordinates without external lattice vectors'''
        infile = self.get_input('h2o.234.xyz')
        reffile = self.get_input('h2o.234.gen')
        outfile = self.get_output('h2o.234.gen')
        cmdargs = ['-o', outfile, infile]
        xyz2gen.main(cmdargs)
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_periodicabs_lattice(self):
        '''Absolute coordinates with external lattice vectors'''
        infile = self.get_input('h2o.xyz')
        latticefile = self.get_input('h2o.latvecs')
        reffile = self.get_input('h2o.gen')
        outfile = self.get_output('h2o.gen')
        cmdargs = ['-l', latticefile, '-o', outfile, infile]
        xyz2gen.main(cmdargs)
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_fraccoords_lattice(self):
        '''Fractional coordinates with external lattice vectors'''
        infile = self.get_input('gaas-frac.xyz')
        latticefile = self.get_input('gaas-frac.latvecs')
        reffile = self.get_input('gaas-frac.gen')
        outfile = self.get_output('gaas-frac.gen')
        cmdargs = ['-o', outfile, '-l', latticefile, '-f', infile]
        xyz2gen.main(cmdargs)
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_fail_invalid_infile(self):
        '''Failing due to invalid input file.'''
        tempname = common.get_temporary_filename(self.workroot)
        nonexisting_infile = os.path.join(self.workdir, tempname)
        outfile = self.get_output('h2o.gen')
        cmdargs = ['-o', outfile, nonexisting_infile]
        with self.assertRaises(ScriptError):
            xyz2gen.main(cmdargs)


if __name__ == '__main__':
    unittest.main()

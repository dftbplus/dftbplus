#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Tests for gen2xyz.'''

import sys
import os.path
import unittest
import common
from dptools.scripts.common import ScriptError
import dptools.scripts.gen2xyz as gen2xyz


SCRIPTDIR = os.path.dirname(sys.argv[0])


class Gen2xyzTest(common.TestWithWorkDir):
    '''General tests for gen2xyz'''

    def setUp(self):
        self.inputdir = os.path.join(SCRIPTDIR, 'gen2xyz')
        self.workroot = './'
        common.TestWithWorkDir.setUp(self)

    def test_cluster(self):
        '''Absolute coordinates without external lattice vectors'''
        infile = self.get_input('h2o.234.gen')
        reffile = self.get_input('h2o.234.xyz')
        outfile = self.get_output('h2o.234.xyz')
        cmdargs = ['-o', outfile, infile]
        gen2xyz.main(cmdargs)
        self.assertTrue(common.xyz_file_equals(outfile, reffile))

    def test_periodicabs_lattice(self):
        '''Absolute coordinates with external lattice vectors'''
        infile = self.get_input('h2o.gen')
        latticefile = self.get_input('h2o.latvecs')
        reffile = self.get_input('h2o.xyz')
        outfile = self.get_output('h2o.xyz')
        cmdargs = ['-l', latticefile, '-o', outfile, infile]
        gen2xyz.main(cmdargs)
        self.assertTrue(common.xyz_file_equals(outfile, reffile))

    def test_commandstr(self):
        '''Absolute coordinates with a comment string different to reffile'''
        infile = self.get_input('h2o.gen')
        reffile = self.get_input('h2o_comment.xyz')
        outfile = self.get_output('h2o_comment.xyz')
        cmdargs = ['-o', outfile, '-c', "something", infile]
        gen2xyz.main(cmdargs)
        self.assertTrue(common.xyz_file_equals(outfile, reffile))
        self.assertFalse(common.xyz_file_equals(outfile, reffile,\
                         check_comment=True))

    def test_fail_invalid_infile(self):
        '''Failing due to invalid input file.'''
        tempname = common.get_temporary_filename(self.workroot)
        nonexisting_infile = os.path.join(self.workdir, tempname)
        outfile = self.get_output('h2o.234.xyz')
        cmdargs = ['-o', outfile, nonexisting_infile]
        with self.assertRaises(ScriptError):
            gen2xyz.main(cmdargs)


if __name__ == '__main__':
    unittest.main()

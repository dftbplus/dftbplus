#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2018  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Tests for repeatgen.'''

import sys
import os.path
import unittest
import common
from dptools.scripts.common import ScriptError
import dptools.scripts.repeatgen as repeatgen


SCRIPTDIR = os.path.dirname(sys.argv[0])


class RepeatgenTest(common.TestWithWorkDir):
    '''General tests for repeatgen'''

    def setUp(self):
        self.inputdir = os.path.join(SCRIPTDIR, 'repeatgen')
        self.workroot = './'
        common.TestWithWorkDir.setUp(self)

    def test_cluster(self):
        '''Cluster with external lattice vector'''
        infile = self.get_input('h2o.gen')
        latticefile = self.get_input('h2o.latvecs')
        reffile = self.get_input('h2o.234.gen')
        outfile = self.get_output('h2o.234.gen')
        cmdargs = ['-o', outfile, '-l', latticefile, infile, '2', '3', '4']
        repeatgen.main(cmdargs)
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_periodic_abscoords(self):
        '''Periodic structure with absolute coordinates'''
        infile = self.get_input('gaas.gen')
        reffile = self.get_input('gaas.234.gen')
        outfile = self.get_output('gaas.234.gen')
        cmdargs = ['-o', outfile, infile, '2', '3', '4']
        repeatgen.main(cmdargs)
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_periodic_relcoords(self):
        '''Periodic structure with relative coordinates'''
        infile = self.get_input('gaas-frac.gen')
        reffile = self.get_input('gaas-frac.234.gen')
        outfile = self.get_output('gaas-frac.234.gen')
        cmdargs = ['-o', outfile, infile, '2', '3', '4']
        repeatgen.main(cmdargs)
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_stdout(self):
        '''When output goes to stdout instead of file.'''
        infile = self.get_input('gaas.gen')
        reffile = self.get_input('gaas.234.gen')
        outfile = self.get_output('gaas.234.gen')
        cmdargs = [infile, '2', '3', '4']
        with common.OutputCatcher() as output:
            repeatgen.main(cmdargs)
        outfile = output.get_as_stringio()
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_fail_cluster_without_latvecs(self):
        '''Failing due to missing lattice vectors'''
        infile = self.get_input('h2o.gen')
        outfile = self.get_output('h2o.234.gen')
        cmdargs = ['-o', outfile, infile, '2', '3', '4']
        with self.assertRaises(ScriptError):
            repeatgen.main(cmdargs)

    def test_fail_wrong_repetition(self):
        '''Failing due to wrong repetition numbers'''
        infile = self.get_input('h2o.gen')
        latticefile = self.get_input('h2o.latvecs')
        outfile = self.get_output('h2o.234.gen')
        cmdargs = ['-o', outfile, '-l', latticefile, infile, '0', '3', '4']
        with self.assertRaises(ScriptError):
            repeatgen.main(cmdargs)
            
    def test_fail_invalid_repetition(self):
        '''Failing due to invalid type of repetition numbers'''
        infile = self.get_input('h2o.gen')
        latticefile = self.get_input('h2o.latvecs')
        outfile = self.get_output('h2o.234.gen')
        cmdargs = ['-o', outfile, '-l', latticefile, infile, 'a', '3', '4']
        with self.assertRaises(ScriptError):
            repeatgen.main(cmdargs)
            
    def test_fail_missing_arguments(self):
        '''Failing due to missing arguments'''
        infile = self.get_input('h2o.gen')
        latticefile = self.get_input('h2o.latvecs')
        outfile = self.get_output('h2o.234.gen')
        cmdargs = ['-o', outfile, '-l', latticefile, infile, '3', '4']
        with self.assertRaises(ScriptError):
            repeatgen.main(cmdargs)


    def test_fail_superfluous_arguments(self):
        '''Failing due to superfluous arguments.'''
        infile = self.get_input('h2o.gen')
        latticefile = self.get_input('h2o.latvecs')
        outfile = self.get_output('h2o.234.gen')
        cmdargs = ['-o', outfile, '-l', latticefile, infile, '3', '4', '5',
                   'something']
        with self.assertRaises(ScriptError):
            repeatgen.main(cmdargs)


if __name__ == '__main__':
    unittest.main()

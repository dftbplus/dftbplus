#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Tests for straingen.'''

import sys
import os.path
import unittest
import common
from dptools.scripts.common import ScriptError
import dptools.scripts.straingen as straingen


SCRIPTDIR = os.path.dirname(sys.argv[0])


class StraingenTest(common.TestWithWorkDir):
    '''General tests for straingen'''

    def setUp(self):
        self.inputdir = os.path.join(SCRIPTDIR, 'straingen')
        self.workroot = './'
        common.TestWithWorkDir.setUp(self)

    def test_clusterStrain(self):
        '''Cluster with positive isotropic strain (using default)'''
        infile = self.get_input('h2o.gen')
        reffile = self.get_input('h2o.isostrain.gen')
        outfile = self.get_output('h2o.isostrain.gen')
        cmdargs = ['-o', outfile, infile, '10']
        straingen.main(cmdargs)
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_clusterStrainNeg(self):
        '''Cluster with negative isotropic strain (using default)'''
        infile = self.get_input('h2o.gen')
        reffile = self.get_input('h2o.negisostrain.gen')
        outfile = self.get_output('h2o.negisostrain.gen')
        cmdargs = ['-o', outfile, infile, '-10']
        straingen.main(cmdargs)
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_clusterStrainOpt(self):
        '''Cluster with positive isotropic strain
           (using default explicit component)'''
        infile = self.get_input('h2o.gen')
        reffile = self.get_input('h2o.isostrain.gen')
        outfile = self.get_output('h2o.isostrain.gen')
        cmdargs = ['-o', outfile, '-c', 'I', infile, '10']
        straingen.main(cmdargs)
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_clusterStrainNegOpt(self):
        '''Cluster with negative isotropic strain (using explicit component)'''
        infile = self.get_input('h2o.gen')
        reffile = self.get_input('h2o.negisostrain.gen')
        outfile = self.get_output('h2o.negisostrain.gen')
        cmdargs = ['-o', outfile, '-c', 'I', infile, '-10']
        straingen.main(cmdargs)
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_clusterStrainYY(self):
        '''Cluster with yy strain'''
        infile = self.get_input('h2o.gen')
        reffile = self.get_input('h2o.yy.gen')
        outfile = self.get_output('h2o.yy.gen')
        cmdargs = ['-o', outfile, '-c', 'yy', infile, '10']
        straingen.main(cmdargs)
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_clusterStrainXZ(self):
        '''Cluster with xz strain'''
        infile = self.get_input('h2o.gen')
        reffile = self.get_input('h2o.xz.gen')
        outfile = self.get_output('h2o.xz.gen')
        cmdargs = ['-o', outfile, '-c', 'xz', infile, '10']
        straingen.main(cmdargs)
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_SuperStrain(self):
        '''Supercell with positive isotropic strain (using default)'''
        infile = self.get_input('gaas.gen')
        reffile = self.get_input('gaas.isostrain.gen')
        outfile = self.get_output('gaas.isostrain.gen')
        cmdargs = ['-o', outfile, infile, '10']
        straingen.main(cmdargs)
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_FracSuperStrain(self):
        '''Supercell (fractional) with positive isotropic strain (using
        default)'''
        infile = self.get_input('gaas-frac.gen')
        reffile = self.get_input('gaas-frac.isostrain.gen')
        outfile = self.get_output('gaas-frac.isostrain.gen')
        cmdargs = ['-o', outfile, infile, '10']
        straingen.main(cmdargs)
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_stdout(self):
        '''When output goes to stdout instead of file.'''
        infile = self.get_input('h2o.gen')
        reffile = self.get_input('h2o.isostrain.gen')
        outfile = self.get_output('h2o.isostrain.gen')
        cmdargs = [infile, '10']
        with common.OutputCatcher() as output:
            straingen.main(cmdargs)
        outfile = output.get_as_stringio()
        self.assertTrue(common.gen_file_equals(outfile, reffile))

    def test_fail_invalid_option(self):
        '''Failing due to invalid option argument'''
        infile = self.get_input('h2o.gen')
        outfile = self.get_output('h2o.isostrain.gen')
        cmdargs = ['-o', outfile, '-c', 'zx', infile, '10']
        with self.assertRaises(ScriptError):
            straingen.main(cmdargs)

    def test_fail_invalid_infile(self):
        '''Failing due to invalid input file.'''
        tempname = common.get_temporary_filename(self.workroot)
        nonexisting_infile = os.path.join(self.workdir, tempname)
        outfile = self.get_output('h2o.noinfile.gen')
        cmdargs = ['-o', outfile, nonexisting_infile, '10']
        with self.assertRaises(ScriptError):
            straingen.main(cmdargs)

if __name__ == '__main__':
    unittest.main()

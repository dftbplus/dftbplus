#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2023  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Common things needed by command line scripts.'''

import numpy as np


class ScriptError(Exception):
    '''Exception thrown by command line scripts.'''


def find_auto_alignment(bandout):
    '''
    Integer occupations: finds the energy shift required to set the VBM to zero
    Fractional occupations: finds Fermi-type level

    Args:
        bandout (BandOut): representation of a band.out like file

    '''

    eigvals = bandout.eigvalarray[:, :, :, 0]
    occ = bandout.eigvalarray[:, :, :, 1]

    return np.max(eigvals[np.where(occ >= np.max(occ) / 2.0)])

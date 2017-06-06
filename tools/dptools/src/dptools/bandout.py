#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2017  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

############################################################################
# Information from band.out like files
############################################################################
import re
import numpy as np
from dptools.common import openfile

__all__ = [ "BandOut", ]

PAT_BLOCK = re.compile(r"""KPT\s*(?P<ikpt>\d+)\s+
                           (?:SPIN\s*(?P<ispin>\d+)\s+)?
                           (?:KWEIGHT\s*
                               (?P<kweight>\d+\.\d+(?:[eE][+-]\d+)?)\s+)?
                           (?P<vals>
                               (?:(?:[+-]?\d+\.\d+(?:[eE][+-]\d+)?\s+){2})+)
                           """, re.VERBOSE | re.MULTILINE)


class BandOut:
    """Representation of a band.out like file.

    Attributes:
        nspin: Number of spin channels.
        nkpt: Number of k-points.
        nstate: Number of states per k-point and spin.
        kweights: Weight of the k-points.
        eigvalarray: (nspin, nkpt, nstate, 2) shaped array containing
            the energy of the states and their occupation.
    """

    def __init__(self, kweights, eigvalarray):
        """Initializes band.out file instance.

        Args:
            kweights: Weights of the kpoints as ( nkpt, ) shaped array.
            eigvalarray: Eigenvalues and additional information (e.g occupation)
                for each spin and kpt as an (npsin, nkpt, nlevel, 2) shaped
                array.
        """
        self.kweights = kweights
        self.eigvalarray = eigvalarray
        self.nspin, self.nkpt, self.nstate = self.eigvalarray.shape[:3]

    @classmethod
    def fromfile(cls, fobj):
        """Returns a band.out representation created from a file object.

        Args:
            fobj: File like object or string with file name.
        """
        fp = openfile(fobj)
        txt = fp.read()
        fp.close()
        ispins = set()
        kweights = []
        eigvalarrays = []
        match = PAT_BLOCK.search(txt)
        while match:
            ispin = match.group("ispin")
            if ispin:
                ispins.add(int(ispin))
            else:
                ispins.add(1)
            kweight = match.group("kweight")
            if kweight:
                kweights.append(float(kweight))
            else:
                kweights.append(1.0)
            tmp = np.array(match.group("vals").split(), dtype=float)
            eigvalarrays.append(tmp.reshape((-1, 2)))
            match = PAT_BLOCK.search(txt, match.end())
        nspin = len(ispins)
        nkpt = len(eigvalarrays) / nspin
        eigvalspin = np.array(eigvalarrays).reshape((nspin, nkpt, -1, 2))
        kweights = np.array(kweights).reshape((nspin, nkpt))
        return cls(kweights, eigvalspin)

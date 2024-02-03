#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2024  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#
'''Information from band.out like files'''

from __future__ import division
import re
import numpy as np
from dptools.common import openfile

__all__ = ["BandOut"]

_PAT_BLOCK = re.compile(r"""KPT\s*(?P<ikpt>\d+)\s+
                            (?:SPIN\s*(?P<ispin>\d+)\s+)?
                            (?:KWEIGHT\s*
                                (?P<kweight>\d+\.\d+(?:[eE][+-]\d+)?)\s+)?
                            (?P<vals>
                                (?:(\d+\s+)?
                                (?:[+-]?\d+\.\d+(?:[eE][+-]\d+)?\s+){2})+)
                            """, re.VERBOSE | re.MULTILINE)


class BandOut:
    """Representation of a band.out like file.

    Attributes:
        nspin: Number of spin channels.
        nkpt: Number of k-points.
        nstate: Number of states per k-point and spin.
        kweights: Weight of the k-points.
        eigvalarray: (nspin, nkpt, nstate, 2 or 5) shaped array containing
            the energy of the states and their occupation (a single number
            or as occ_q, occ_mx, occ_my, occ_my for non-collinear spin)
    """

    def __init__(self, kweights, eigvalarray):
        """Initializes band.out file instance.

        Args:
            kweights: Weights of the kpoints as ( nkpt, ) shaped array.
            eigvalarray: Eigenvalues and additional information (e.g occupation)
                for each spin and kpt as an (npsin, nkpt, nlevel, 2) shaped
                array. Alternatively, (for non-collinear spin) 4 values for the
                additional information (occ_q, occ_mx, occ_my, occ_mz) can be
                used by passing an array of the shape (nspin, nkpt, nlevel, 5).
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
        ndatacols = 0
        for match in  _PAT_BLOCK.finditer(txt):
            ispin = match.group("ispin")
            ispins.add(int(ispin) if ispin is not None else 1)
            kweight = match.group("kweight")
            kweights.append(float(kweight) if kweight is not None else 1.0)

            vals = match.group("vals")
            tmp = np.array(vals.split(), dtype=float)
            nrows = vals.strip().count('\n') + 1
            ncols = len(tmp) // nrows

            # Optional first column might contain sequential numbering.
            # Then either 2 columns (energy, occ) follow for spin-unpol and
            # collinear cases or 5 columns (energy, occ_q, occ_mx, occ_my,
            # occ_mz) for the non-collinear case.
            if ncols not in (2, 3, 5, 6):
                raise ValueError("Bandout file contains invalid nr. of columns")
            startcol = 0 if ncols == 2 or ncols == 5 else 1
            ndatacols = 2 if ncols == 2 or ncols == 3 else 5

            tmp.shape = (nrows, ncols)
            tmp = tmp[:, startcol : startcol + ndatacols]
            eigvalarrays.append(tmp)

        nspin = len(ispins)
        nkpt = len(eigvalarrays) // nspin
        eigvalspin = np.array(eigvalarrays)
        eigvalspin.shape = (nspin, nkpt, -1, ndatacols)
        kweights = np.array(kweights)
        kweights.shape = (nspin, nkpt)
        return cls(kweights, eigvalspin)

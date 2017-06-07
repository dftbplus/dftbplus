#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2017  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

############################################################################
# Various commonly used routines.
############################################################################
import gzip

def openfile(fobj, mode=None):
    """Opens a file or a file like object.

    Args:
        fobj: File name or file like object.
        mode: File access mode (default: 'r')

    Returns:
        If a string (file name) was provided for fobj, a pointer to the
        opened file. Otherwise the original file like object.
    """
    if mode is None:
        mode = "r"
    if isinstance(fobj, str):
        if fobj.endswith(".gz"):
            fp = gzip.open(fobj, mode)
        else:
            fp = open(fobj, mode)
    else:
        fp = fobj
    return fp

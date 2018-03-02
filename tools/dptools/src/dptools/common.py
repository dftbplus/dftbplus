#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2018  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Various commonly used items of the dptools package.'''

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



class OpenFile:
    '''Represents an open file.

    It can either transparently pass through an already opened file descriptor
    (file like object) or open a file itself and close it at exit.
    '''

    def __init__(self, fobj, mode=None):
        '''Initialises an open file.

        Args:
            fobj (str or file object): File to open. If it is a string, it
                represents the name of the file, otherwise it should be a file
                like object. If it is a file name and ends with '.gz' the file
                is opened as a gzipped file.
            mode (str): Opening mode for the file (provided fobj represents a
                file name).
        '''
        self._fobj = fobj
        self._mode = mode if mode is not None else 'r'
        self._open = None
        if isinstance(fobj, str):
            if fobj.endswith(".gz"):
                self._open = gzip.open
            else:
                self._open = open
        self._fp = None


    def __enter__(self):
        if self._open is not None:
            self._fp = self._open(self._fobj, self._mode)
        else:
            self._fp = self._fobj
        return self._fp


    def __exit__(self, *args):
        if self._open is not None:
            self._fp.close()

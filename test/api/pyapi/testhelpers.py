#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

'''
An (expandable) collection of convenient methods to
run regression tests of the ctypes Python interface.
'''


import numpy as np


# relevant tags for testing the Python interface
TAGLABELS = {'freeEgy': 'mermin_energy',
             'forceTot': 'forces',
             'qOutAtGross': 'gross_atomic_charges'}


def write_autotest_tag(filename, **kwargs):
    '''Writes an autotest.tag file, containing submitted results. The
       created file can then be compared with a validated reference file.

    Args:
        filename  (str): path + name of autotest.tag file
        **kwargs (dict): keyworded dictionary to pass var-
                         iable number of keyword arguments

    '''

    # get number of chars of longest tag (used to format output)
    maxtaglen = max([len(TAGLABELS[key]) for key in kwargs])

    tagfile = open(filename, 'w')

    for key, value in kwargs.items():

        content = []
        result = np.asarray(value)
        flatres = result.flatten()
        label = TAGLABELS[key]

        # build up string to write
        content.append(label)
        content.append(' ' * (maxtaglen - len(label) + 1))

        # check if array is real
        # (because type 'real' is hardcoded at this time)
        # (the interface only outputs real quantities)
        if np.isrealobj(result):
            content.append(':real:')
        else:
            msg = 'The type "{}" of '.format(result.dtype) + \
                  'kwarg "{}" is invalid.'.format(key)
            raise TypeError(msg)

        content.append(str(result.ndim) + ':')
        content.append(','.join(
            [str(entry) for entry in reversed(result.shape)]))
        content.append('\n')

        nentries = len(flatres)
        formatstr = ' '.join([' {:24.15E}', ] * 3 + ['\n'])
        for index in range(2, nentries, 3):
            values = formatstr.format(flatres[index - 2],
                                      flatres[index - 1],
                                      flatres[index])
            content.append(values)

        # handle possible remaining entries
        nremain = nentries % 3
        formatstr = ' '.join([' {:24.15E}', ] * nremain)
        if nremain == 2:
            values = formatstr.format(flatres[-2], flatres[-1])
            content.append(values)
            content.append('\n')
        elif nremain == 1:
            values = formatstr.format(flatres[-1])
            content.append(values)
            content.append('\n')
        else:
            values = ''

        tagfile.write(''.join(content))

    tagfile.close()

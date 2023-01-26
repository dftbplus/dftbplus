#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Common items for the ptools tests framework.'''

import numpy as np
from common_setup import TestWithWorkDir


def type_diff(reference, compare, rtol=1e-05, atol=1e-08):
    """common function for comparing dictionaries

    Args:
        reference (): reference
        compare (): to compare with reference
        rtol (float): relative tolerance
        atol (float): absolute tolerance

    Returns:
        (bool): True if equal, False if not equal
    """
    if not isinstance(compare, type(reference)):
        return False

    if isinstance(reference, dict):
        if set(reference.keys())-set(compare.keys()):
            return False
        for key in reference.keys():
            if not type_diff(reference[key], compare[key], rtol=rtol,
                             atol=atol):
                return False

    elif isinstance(reference, list):
        if not len(reference) == len(compare):
            return False
        for item1, item2 in zip(reference, compare):
            if not type_diff(item1, item2, rtol=rtol, atol=atol):
                return False

    elif isinstance(reference, float):
        if not np.isclose(compare, reference, rtol=rtol, atol=atol):
            return False

    elif isinstance(reference, complex):
        if not np.isclose(compare, reference, rtol=rtol, atol=atol):
            return False

    elif isinstance(reference, np.ndarray):
        if not np.shape(reference) == np.shape(compare):
            return False
        if not reference.dtype == compare.dtype:
            return False
        if not np.allclose(compare, reference, rtol=rtol, atol=atol):
            return False

    elif isinstance(reference, str):
        if not reference.lower() == compare.lower():
            return False

    else:
        if not reference == compare:
            return False

    return True

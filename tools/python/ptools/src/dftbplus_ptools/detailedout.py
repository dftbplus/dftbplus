#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

"""
Module for interaction with the detail.out file.
"""

import os


class Detailedout:
    """class for reading detailed.out"""

    def __init__(self, directory=".", filenname="detailed.out"):
        """Initialises the Detailedout class

        Args:
            directory (str): directory of file
            filename (str): name of file
        """
        self._directory = directory
        self._filename = filenname


    def check_for_scc_convergence(self):
        """function for checking scc convergence

        Returns:
            (bool/None): 'True/False' if converged or not, 'None' if file not
                found/readable
        """
        path = os.path.join(self._directory, self._filename)
        try:
            with open(path, "r") as file:
                lines = file.readlines()
        except OSError:
            return None

        for line in lines:
            if line.startswith("SCC converged"):
                return True

        return False


    def check_for_geometry_convergence(self):
        """function for checking geometry convergence

        Returns:
            (bool/None): 'True/False' if converged or not, 'None' if file not
                found/readable
        """
        path = os.path.join(self._directory, self._filename)
        try:
            with open(path, "r") as file:
                lines = file.readlines()
        except OSError:
            return None

        for line in lines:
            if line.startswith("Geometry converged"):
                return True

        return False

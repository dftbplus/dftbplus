#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

"""
Module for interaction with the detail.out file.
"""

from __future__ import annotations


class Detailedout:
    """class for reading detailed.out"""

    def __init__(self, filename: str = "detailed.out") -> None:
        """Initialises the Detailedout class

        Args:
            filename (str): name of file
        """
        self._filename = filename


    def check_for_scc_convergence(self) -> bool | None:
        """function for checking scc convergence

        Returns:
            (bool/None): 'True/False' if converged or not, 'None' if file not
                found/readable
        """
        try:
            with open(self._filename, "r") as file:
                lines = file.readlines()
        except OSError:
            return None

        for line in lines:
            if line.startswith("SCC converged"):
                return True

        return False


    def check_for_geometry_convergence(self) -> bool | None:
        """function for checking geometry convergence

        Returns:
            (bool/None): 'True/False' if converged or not, 'None' if file not
                found/readable
        """
        try:
            with open(self._filename, "r") as file:
                lines = file.readlines()
        except OSError:
            return None

        for line in lines:
            if line.startswith("Geometry converged"):
                return True

        return False

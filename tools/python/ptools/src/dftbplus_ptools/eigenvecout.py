#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

"""
Module for interaction with the eigenvec.out file.
"""


import re


class Eigenvecout:
    """class for reading eigenvectors"""

    def __init__(self, filename="eigenvec.out"):
        """Initialises the Eigenvecout class

        Args:
            filename (str): name of file
        """
        self._filename = filename


    def read_eigenvec(self):
        """function for reading eigenvector

        Returns:
            dictionary (dict/None): dictionary containing the eigenvector or
                'None' if file not readable/found
        """

        try:
            with open(self._filename, "r") as file:
                lines = file.readlines()
        except OSError:
            return None

        dictionary = {}
        kpoint = None
        eigennum = None
        spin = None
        atom = None
        orbit = None

        pattern1 = re.compile(r"""^K-point:\s+
                              (?P<kpoint>[0-9]+)\s+
                              Eigenvector:\s+
                              (?P<eigennum>[0-9]+)\s+
                              \((?P<spin>[a-z]+)\)
                              """, re.VERBOSE)
        pattern2 = re.compile(r"""^Eigenvector:\s+
                              (?P<eigennum>[0-9]+)\s+
                              \((?P<spin>[a-z]+)\)
                              """, re.VERBOSE)
        pattern3 = re.compile(r"""^\s+(?P<atom>[0-9 a-zA-Z]+)\s+
                              (?P<orbit>s)\s+\(
                              (?P<eigenvec1>\s*[-0-9.]+),
                              (?P<eigenvec2>\s*[-0-9.]+)
                              \s*\)
                              \s+(?P<mulliken>[-0-9.]+)
                              """, re.VERBOSE)
        pattern4 = re.compile(r"""^\s+
                              (?P<atom>[0-9 a-zA-Z]+)\s+
                              (?P<orbit>s)\s+
                              (?P<eigenvec>[-0-9.]+)\s+
                              (?P<mulliken>[-0-9.]+)
                              """, re.VERBOSE)
        pattern5 = re.compile(r"""^\s+
                              (?P<orbit>[-a-zA-Z_0-9]+)\s+\(
                              (?P<eigenvec1>\s*[-0-9.]+),
                              (?P<eigenvec2>\s*[-0-9.]+)
                              \s*\)
                              \s+(?P<mulliken>[-0-9.]+)
                              """, re.VERBOSE)
        pattern6 = re.compile(r"""^\s+(?P<orbit>[-a-zA-Z_0-9]+)
                              \s+(?P<eigenvec>[-0-9.]+)
                              \s+(?P<mulliken>[-0-9.]+)
                              """, re.VERBOSE)
        pattern7 = re.compile(r"""^K-point:\s+
                              (?P<kpoint>[0-9]+)\s+
                              Eigenvector:\s+
                              (?P<eigennum>[0-9]+)
                               """, re.VERBOSE)
        pattern8 = re.compile(r"""^\s+(?P<atom>[0-9 a-zA-Z]+)\s+
                              (?P<orbit>s)\s+\(
                              (?P<eigenvec1>\s*[-0-9.]+),
                              (?P<eigenvec2>\s*[-0-9.]+)
                              \s*\)\(
                              (?P<eigenvec3>\s*[-0-9.]+),
                              (?P<eigenvec4>\s*[-0-9.]+)
                              \s*\)
                              \s+(?P<mulliken>[-0-9.]+)
                              \s+(?P<x>[-0-9.]+)
                              \s+(?P<y>[-0-9.]+)
                              \s+(?P<z>[-0-9.]+)
                              """, re.VERBOSE)
        pattern9 = re.compile(r"""^\s+
                              (?P<orbit>[-a-zA-Z_0-9]+)\s+\(
                              (?P<eigenvec1>\s*[-0-9.]+),
                              (?P<eigenvec2>\s*[-0-9.]+)
                              (?P<eigenvec3>\s*[-0-9.]+),
                              (?P<eigenvec4>\s*[-0-9.]+)
                              \s*\)
                              \s+(?P<mulliken>[-0-9.]+)
                              \s+(?P<x>[-0-9.]+)
                              \s+(?P<y>[-0-9.]+)
                              \s+(?P<z>[-0-9.]+)
                              """, re.VERBOSE)

        for line in lines:

            if pattern5:
                match = pattern5.match(line)
            else:
                match = None
            if match is not None:
                orbit = match.group("orbit")
                eigenvec1 = float(match.group("eigenvec1"))
                eigenvec2 = float(match.group("eigenvec2"))
                eigenvec = complex(eigenvec1, eigenvec2)
                mulliken = float(match.group("mulliken"))

                dictionary[kpoint][eigennum][spin][atom][orbit] = {}
                dictionary[kpoint][eigennum][spin][atom][orbit]["eigenvec"] = \
                    eigenvec
                dictionary[kpoint][eigennum][spin][atom][orbit]["mulliken"] = \
                    mulliken
                pattern6 = False
                pattern9 = False
                continue

            if pattern6:
                match = pattern6.match(line)
            else:
                match = None
            if match is not None:
                orbit = match.group("orbit")
                eigenvec = float(match.group("eigenvec"))
                mulliken = float(match.group("mulliken"))

                dictionary[kpoint][eigennum][spin][atom][orbit] = {}
                dictionary[kpoint][eigennum][spin][atom][orbit]["eigenvec"] = \
                    eigenvec
                dictionary[kpoint][eigennum][spin][atom][orbit]["mulliken"] = \
                    mulliken
                pattern5 = False
                pattern9 = False
                continue

            if pattern9:
                match = pattern9.match(line)
            else:
                match = None
            if match is not None:
                orbit = match.group("orbit")
                eigenvec1 = float(match.group("eigenvec1"))
                eigenvec2 = float(match.group("eigenvec2"))
                eigenvec_up = complex(eigenvec1, eigenvec2)
                eigenvec3 = float(match.group("eigenvec3"))
                eigenvec4 = float(match.group("eigenvec4"))
                eigenvec_down = complex(eigenvec3, eigenvec4)
                mulliken = float(match.group("mulliken"))
                x = float(match.group("x"))
                y = float(match.group("y"))
                z = float(match.group("z"))

                dictionary[kpoint][eigennum][atom][orbit] = {}
                dictionary[kpoint][eigennum][atom][orbit]["eigenvec_up"] = \
                    eigenvec_up
                dictionary[kpoint][eigennum][atom][orbit]["eigenvec_down"] = \
                    eigenvec_down
                dictionary[kpoint][eigennum][atom][orbit]["mulliken"] = \
                    mulliken
                dictionary[kpoint][eigennum][atom][orbit]["x"] = x
                dictionary[kpoint][eigennum][atom][orbit]["y"] = y
                dictionary[kpoint][eigennum][atom][orbit]["z"] = y
                pattern5 = False
                pattern6 = False
                continue

            if pattern3:
                match = pattern3.match(line)
            else:
                match = None
            if match is not None:
                atom = match.group("atom").replace(" ", "")
                orbit = match.group("orbit")
                eigenvec1 = float(match.group("eigenvec1"))
                eigenvec2 = float(match.group("eigenvec2"))
                eigenvec = complex(eigenvec1, eigenvec2)
                mulliken = float(match.group("mulliken"))

                dictionary[kpoint][eigennum][spin][atom] = {}
                dictionary[kpoint][eigennum][spin][atom][orbit] = {}
                dictionary[kpoint][eigennum][spin][atom][orbit]["eigenvec"] = \
                    eigenvec
                dictionary[kpoint][eigennum][spin][atom][orbit]["mulliken"] = \
                    mulliken
                pattern4 = False
                pattern8 = False
                continue

            if pattern4:
                match = pattern4.match(line)
            else:
                match = None
            if match is not None:
                atom = match.group("atom").replace(" ", "")
                orbit = match.group("orbit")
                eigenvec = float(match.group("eigenvec"))
                mulliken = float(match.group("mulliken"))

                dictionary[kpoint][eigennum][spin][atom] = {}
                dictionary[kpoint][eigennum][spin][atom][orbit] = {}
                dictionary[kpoint][eigennum][spin][atom][orbit]["eigenvec"] = \
                    eigenvec
                dictionary[kpoint][eigennum][spin][atom][orbit]["mulliken"] = \
                    mulliken
                pattern3 = False
                pattern8 = False
                continue

            if pattern8:
                match = pattern8.match(line)
            else:
                match = None
            if match is not None:
                atom = match.group("atom").replace(" ", "")
                orbit = match.group("orbit")
                eigenvec1 = float(match.group("eigenvec1"))
                eigenvec2 = float(match.group("eigenvec2"))
                eigenvec_up = complex(eigenvec1, eigenvec2)
                eigenvec3 = float(match.group("eigenvec3"))
                eigenvec4 = float(match.group("eigenvec4"))
                eigenvec_down = complex(eigenvec3, eigenvec4)
                mulliken = float(match.group("mulliken"))
                x = float(match.group("x"))
                y = float(match.group("y"))
                z = float(match.group("z"))

                dictionary[kpoint][eigennum][atom] = {}
                dictionary[kpoint][eigennum][atom][orbit] = {}
                dictionary[kpoint][eigennum][atom][orbit]["eigenvec_up"] = \
                    eigenvec_up
                dictionary[kpoint][eigennum][atom][orbit]["eigenvec_down"] = \
                    eigenvec_down
                dictionary[kpoint][eigennum][atom][orbit]["mulliken"] = \
                    mulliken
                dictionary[kpoint][eigennum][atom][orbit]["x"] = x
                dictionary[kpoint][eigennum][atom][orbit]["y"] = y
                dictionary[kpoint][eigennum][atom][orbit]["z"] = z
                pattern3 = False
                pattern4 = False
                continue

            if pattern1:
                match = pattern1.match(line)
            else:
                match = None
            if match is not None:
                kpoint = int(match.group("kpoint"))
                eigennum = int(match.group("eigennum"))
                spin = match.group("spin")
                if kpoint not in dictionary:
                    dictionary[kpoint] = {}
                if eigennum not in dictionary[kpoint]:
                    dictionary[kpoint][eigennum] = {}
                if spin not in dictionary[kpoint][eigennum]:
                    dictionary[kpoint][eigennum][spin] = {}
                pattern2 = False
                pattern7 = False
                continue

            if pattern2:
                match = pattern2.match(line)
            else:
                match = None
            if match is not None:
                kpoint = 1
                eigennum = int(match.group("eigennum"))
                spin = match.group("spin")
                if kpoint not in dictionary:
                    dictionary[kpoint] = {}
                if eigennum not in dictionary[kpoint]:
                    dictionary[kpoint][eigennum] = {}
                if spin not in dictionary[kpoint][eigennum]:
                    dictionary[kpoint][eigennum][spin] = {}
                pattern1 = False
                pattern7 = False
                continue

            if pattern7:
                match = pattern7.match(line)
            else:
                match = None
            if match is not None:
                kpoint = int(match.group("kpoint"))
                eigennum = int(match.group("eigennum"))
                if kpoint not in dictionary:
                    dictionary[kpoint] = {}
                if eigennum not in dictionary[kpoint]:
                    dictionary[kpoint][eigennum] = {}
                pattern1 = False
                pattern2 = False
                continue

        return dictionary

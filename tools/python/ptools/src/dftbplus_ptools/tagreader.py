#!/usr/bin/env python3
# -*- Mode: Python -*-
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#
############################################################################
#
#  tagreader -- contains a parser for the tagged output of DFTB
#
############################################################################
#
#  (The tagged data is assumed to be represented in the form as provided by
#  the taggedout module of the DFTB project. See the appropriate source code
#  for the details.)
#
############################################################################

"""contains classes for parsing tagged output of DFTB"""

import re
import numpy as np


############################################################################
# Exceptions
############################################################################

class InvalidEntry(Exception):
    """Raised if entry can't be initialized or wrong type of entry is passed"""


    def __init__(self, start=0, end=0, msg=""):
        """initializes InvalidEntry class

        Args:
            start (int): starting line of the block where the error occured
            end (int): first line after ending line of the block where the
                error occured
            msg (str): reason of the error
        """
        super().__init__()
        self.start = start
        self.end = end
        self.msg = msg


class ConversionError(Exception):
    """Raised if error occurs during conversion from string"""


############################################################################
# Conversion functors
############################################################################

class Converter():
    """Base class for string to value converters"""

    def __init__(self, nolist=False):
        """initializes Converter class

        Args:
            nolist (bool): if True, not a list, but only a single value is
                returned after conversion.
        """
        self.nolist = nolist


    def __call__(self, str_value):
        """function call operator for Converter class

        Args:
            str_value (str): string represenation of the values to convert

        Returns:
            result (list): Plain 1D list of converted values, unless nolist
                has been set to True
        """
        values = str_value.split()
        if self.nolist and len(values) > 1:
            raise ConversionError("Too many values")
        result = self.convert(values)
        if self.nolist:
            return result[0]
        return result

    @staticmethod
    def convert(values):
        """Conversion function

        Args:
            values (list): list of strings representation of values to convert

        Returns:
            values (list): list of strings representation of values to convert
        """
        return values


class FloatConverter(Converter):
    """Converts string to float"""
    @staticmethod
    def convert(values):
        """Conversion function

        Args:
            values (list): list of strings representation of values to convert

        Returns:
            result (list): list of converted values

        Raises:
            ConversionError: if unable to convert to float
        """
        result = []
        for val in values:
            try:
                result.append(float(val))
            except Exception as error:
                raise ConversionError(f"Unable to convert float '{val}'") \
                    from error
        return result


class IntConverter(Converter):
    """Converts string to integer"""
    @staticmethod
    def convert(values):
        """Conversion function

        Args:
            values (list): list of strings representation of values to convert

        Returns:
            result (list): list of converted values

        Raises:
            ConversionError: if unable to convert to int
        """
        result = []
        for val in values:
            try:
                result.append(int(val))
            except Exception as error:
                raise ConversionError(f"Unable to convert integer '{val}'") \
                    from error
        return result


class ComplexConverter(Converter):
    """Converts string to complex"""

    def __call__(self, str_value):
        """function call operator for ComplexConverter class

        Args:
            str_value (str): string represenation of the values to convert

        Returns:
            result (list): Plain 1D list of converted values, unless nolist
                has been set to True

        Raises:
            ConversionError: if unable to convert to complex because of to
                many values or odd number of values
        """
        values = str_value.split()
        if len(values) % 2:
            raise ConversionError("Odd number of values")
        if self.nolist and len(values) != 2:
            raise ConversionError("Too many values")
        result = self.convert(values)
        if self.nolist:
            return result[0]
        return result

    @staticmethod
    def convert(values):
        """Conversion function

        Args:
            values (list): list of strings representation of values to convert

        Returns:
            result (list): list of converted values

        Raises:
            ConversionError: if unable to convert to complex
        """
        result = []
        for num in range(0, len(values), 2):
            try:
                result.append(complex(float(values[num]),
                                      float(values[num+1])))
            except Exception as error:
                raise ConversionError("Unable to convert complex "
                                      f"'({values[num]},{values[num+1]})'") \
                                          from error

        return result


class LogicalConverter(Converter):
    """Converts string to logical"""

    def convert(self, values):
        """Conversion function

        Args:
            values (list): list of strings representation of values to convert

        Returns:
            result (list): list of converted values

        Raises:
            ConversionError: if unable to convert to boolean (0 or 1)
        """
        result = []
        for val in values:
            if val in ('T', 't'):
                result.append(1)
            elif val in ('F', 'f'):
                result.append(0)
            else:
                raise ConversionError(f"Unable to convert logical '{val}'")
        return result


############################################################################
# Tagged data related objects
############################################################################

class TaggedEntry():
    """Represents a tagged entry with data."""

    # Converter from string for different types
    _strToValue = {"integer": IntConverter(), "real": FloatConverter(),
                   "complex": ComplexConverter(), "logical": LogicalConverter()
                   }

    # Valid types
    _validTypes = list(_strToValue.keys())


    def __init__(self, name, dtype, rank, shape, str_value):
        """Instantiates an TaggedEntry object.

        Args:
            name (str): name of the tagged entry
            dtype (str): type of the data
            rank (int): rank of the data
            shape (tuple): shape of the data (as tuple)
            str_value (str): value as str

        Raises:
            InvalidEntry: if conversion is not possible or rank and shape
                are incompatible or number of values is invalid to shape
        """

        if dtype not in self._validTypes:
            raise InvalidEntry(msg=f"Invalid data type '{dtype}'")
        self._name = name
        self._dtype = dtype
        self._rank = rank
        self._shape = shape
        try:
            self._value = self._strToValue[dtype](str_value)
        except ConversionError as msg:
            raise InvalidEntry(msg=msg) from msg
        if (shape and len(shape) != rank) or (not shape and rank != 0):
            raise InvalidEntry(msg="Incompatible rank and shape")
        if shape and (len(self._value) != product(shape)):
            raise InvalidEntry(msg="Invalid nr. of values")

    @property
    def name(self):
        """Returns name of the entry

        Returns:
            self._name (str): name of the entry
        """
        return self._name

    @property
    def dtype(self):
        """Returns type of the data in the entry

        Returns:
            self._dtype (str): type of the data in the entry
        """
        return self._dtype

    @property
    def rank(self):
        """Returns rank of the data in the entry

        Returns:
            self._rank (int): rank of the data in the entry
        """
        return self._rank

    @property
    def shape(self):
        """Returns shape of the data in the entry

        Returns:
            self._shape (tuple): shape of the data in the entry
        """
        return self._shape

    @property
    def value(self):
        """Returns value of the data in the entry

        Returns:
            self._value (list): value of the data in the entry
        """
        return self._value

    @property
    def fvalue(self):
        """Returns value of the data in the entry as formated array

        Returns:
            result (array): value of the data in the entry
        """
        if self.dtype == 'real':
            dtype = np.float64
        elif self.dtype == 'integer':
            dtype = np.int64
        elif self.dtype == 'complex':
            dtype = np.clongdouble
        elif self.dtype == 'logical':
            dtype = bool

        # converting from column- to row-major
        if self.rank == 0:
            shape_py = (1,)
        else:
            shape_fortran = self.shape
            shape_py = shape_fortran[::-1]

        values_array = np.asarray(self.value, dtype)
        result = np.reshape(values_array, shape_py)

        return result


    def is_comparable(self, other):
        """Check if two entries are comparable

        Args:
            other (tagreader.TaggedEntry): object to compare

        Returns:
            (bool): True if comparable
        """

        return (other.name == self.name and other.dtype == self.dtype
                and other.rank == self.rank and other.shape == self.shape)


class TaggedCollection():
    """Contains a collection of tagged entries"""

    def __init__(self, entries):
        """initializes TaggedCollection class

        Args:
            entries (list): entries to add
        """

        self._entry_names = []
        self._entry_lines = []
        self._entries = []
        self.add_entries(entries)


    def add_entries(self, entries):
        """adding entries to lists in self

        Args:
            entries (list): entries to add
        """

        for entry in entries:
            tagged_line = ":".join((entry.name, entry.dtype, str(entry.rank),
                                   ",".join(map(str, entry.shape))))
            self._entry_names.append(entry.name)
            self._entry_lines.append(tagged_line)
            self._entries.append(entry)


    def get_matching_entries(self, pattern):
        """Returns entries from the collection matching a given pattern

        Args:
            pattern (re.Pattern): compiled regular expression

        Returns:
            result (list): entries matching pattern
        """

        result = []
        for ientry, entry in enumerate(self._entries):
            if pattern.match(self._entry_lines[ientry]):
                result.append(entry)

        return result


    def get_entry(self, name):
        """Returns an entry with a given name from the collection

        Args:
            name (str): name of the entry

        Returns:
            result (TaggedEntry): entry with the name "name"
        """

        try:
            ientry = self._entry_names.index(name)
        except ValueError:
            result = None
        else:
            result = self._entries[ientry]

        return result


    def del_entry(self, name):
        """Deletes the specified entry from the collection

        Args:
            name (str): name of the entry
        """

        try:
            ientry = self._entry_names.index(name)
        except ValueError:
            pass
        else:
            del self._entries[ientry]
            del self._entry_names[ientry]
            del self._entry_lines[ientry]


############################################################################
# Parses the file containing the data and returns TaggedEntries
############################################################################

class ResultParser():
    """Parser the result files containing tagged data"""

    # Pattern for lines containing the describing tag for following data
    pat_tag_line = re.compile(r"""(?P<name>[^: ]+)\s*:
                                (?P<dtype>[^:]+):
                                (?P<rank>\d):
                                (?P<shape>(?:\d+(?:,\d+)*)*)
                            """, re.VERBOSE)


    def __init__(self, file):
        """initializes ResultParser class

        Args:
            file (_io.TextIOWrapper): file like object containing tagged data
        """
        self._file = file


    def iterate_entries(self):
        """Generator for iterating over the entries of the data file.

        Returns:
            result (list): entry of read file

        Raises:
            InvalidEntry: if TaggedEntry is not creatable
        """

        name = None
        dtype = None
        rank = None
        shape = None
        value = []
        iline = 0
        itagged_line = 0
        result = []
        for line in self._file.readlines():
            iline = iline + 1
            match = self.pat_tag_line.match(line)
            if match:

                # Append data from previous tag if present
                if name:
                    try:
                        result.append(TaggedEntry(name, dtype, rank, shape,
                                          " ".join(value)))
                    except InvalidEntry as error:
                        raise InvalidEntry(itagged_line + 1, iline,
                                           msg=error.msg) from error

                name = match.group("name")
                dtype = match.group("dtype")
                rank = int(match.group("rank"))
                if rank > 0:
                    shape = tuple(map(int, match.group("shape").split(",")))
                else:
                    shape = ()
                value = []
                itagged_line = iline
            else:
                value.append(line)

        # process last entry
        if name:
            try:
                result.append(TaggedEntry(name, dtype, rank, shape,
                                          " ".join(value)))
            except InvalidEntry as error:
                raise InvalidEntry(itagged_line + 1, iline, msg=error.msg) \
                    from error

        return result

    entries = property(iterate_entries, None, None,
                       "Iterator over parsed entries")


def product(elements):
    """calculates form size of array the number of elements in array

    Args:
       elements (tuple): tuple with size information of array

    Returns:
        res (int): number of elements in array
    """
    res = 1
    for elem in elements:
        res *= elem

    return res


def results_access(filename='results.tag'):
    """returns the content of a results file as a dictionary

    Args:
        filename (str): name of file

    Retruns:
        dictionary (dict): dictionary of opend results file
    """
    dictionary = {}
    with open(filename, "r") as file:
        parser = ResultParser(file).entries
        for entry in parser:
            values = entry.fvalue
            dictionary[entry.name] = values

    return dictionary


def property_by_keyword(keyword, filename='results.tag'):
    """returns the content of a results file as a dictionary

    Args:
        keyword (str): name of the entry
        filename (str): name of file

    Retruns:
        value (array): array containing values of entry
    """
    with open(filename, "r") as file:
        tagcollection = TaggedCollection(ResultParser(file).entries)

    value = tagcollection.get_entry(keyword).fvalue

    return value

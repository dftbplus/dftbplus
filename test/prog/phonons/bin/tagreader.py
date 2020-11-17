#!/usr/bin/env python3
# -*- Mode: Python -*-
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
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

import re


############################################################################
# Exceptions
############################################################################

class InvalidEntry(Exception):
    """Raised if entry can't be initialized or wrong type of entry is passed"""
    def __init__(self, start=0, end=0, msg=""):
        """start -- starting line of the block where the error occured
        end -- first line after ending line of the block where the error occured
        msg -- reason of the error"""
        self.start = start
        self.end = end
        self.msg = msg

class ConversionError(Exception):
    """Raised if error occurs during conversion from string"""
    pass


############################################################################
# Conversion functors
############################################################################

class Converter(object):
    """Base class for string to value converters"""

    def __init__(self, nolist=False):
        """nolist -- if True, not a list, but only a single value is returned
                     after conversion.
        """
        self.nolist = nolist


    def __call__(self, strValue):
        """strValue -- string represenation of the values to convert
        return -- Plain 1D list of converted values, unless noList has been set
            to True
        """
        values = strValue.split()
        if self.nolist and len(values) > 1:
            raise ConversionError("Too many values")
        result = self.convert(values)
        if self.nolist:
            return result[0]
        else:
            return result


    def convert(self, values):
        """Conversion function
        values -- list of strings representation of values to convert
        """
        return values



class FloatConverter(Converter):
    """Converts string to float"""

    def convert(self, values):
        ll = []
        for val in values:
            try:
                ll.append(float(val))
            except Exception:
                raise ConversionError("Unable to convert float '%s'" % val)
        return ll



class IntConverter(Converter):
    """Converts string to integer"""

    def convert(self, values):
        ll = []
        for val in values:
            try:
                ll.append(int(val))
            except Exception:
                raise ConversionError("Unable to convert integer '%s'" % val)
        return ll



class ComplexConverter(Converter):
    """Converts string to complex"""

    def __call__(self, strValue):
        values = strValue.split()
        if len(values) % 2:
            raise ConversionError("Odd number of values")
        if self.nolist and len(values) != 2:
            raise ConversionError("Too many values")
        result = self.convert(values)
        if self.nolist:
            return result[0]
        else:
            return result


    def convert(self, values):
        ll = []
        for ii in range(0, len(values), 2):
            try:
                ll.append(complex(float(values[ii]), float(values[ii+1])))
            except Exception:
                raise ConversionError("Unable to convert complex '(%s,%s)'"
                                      % (values[ii], values[ii+1]))
        return ll



class LogicalConverter(Converter):
    """Converts string to logical"""

    def convert(self, values):
        ll = []
        for val in values:
            if val == 'T' or val == 't':
                ll.append(1)
            elif val == 'F' or val == 'f':
                ll.append(0)
            else:
                raise ConversionError("Unable to convert logical '%s'" % val)
        return ll



############################################################################
# Tagged data related objects
############################################################################

class TaggedEntry(object):
    """Represents a tagged entry with data."""

    # Converter from string for different types
    __strToValue = { "integer" : IntConverter(),
                                     "real"        : FloatConverter(),
                                     "complex" : ComplexConverter(),
                                     "logical" : LogicalConverter()
                                     }

    # Valid types
    __validTypes = list(__strToValue.keys())


    def __init__(self, name, type, rank, shape, strValue):
        """Instantiates an TaggedEntry object.
        name         -- name of the tagged entry
        type         -- type of the data
        rank         -- rank of the data
        shape        -- shape of the data (as tuple)
        strValue --
        """

        if not type in self.__validTypes:
            raise InvalidEntry(msg="Invalid data type '%s'" % type)
        self.__name = name
        self.__type = type
        self.__rank = rank
        self.__shape = shape
        try:
            self.__value = self.__strToValue[type](strValue)
        except ConversionError as msg:
            raise InvalidEntry(msg=msg)
        if (shape and len(shape) != rank) or (not shape and rank != 0):
            raise InvalidEntry(msg="Incompatible rank and shape")
        if shape and (len(self.__value) != product(shape)):
            raise InvalidEntry(msg="Invalid nr. of values")



    def getName(self):
        return self.__name
    name = property(getName, None, None, "name of the entry")


    def getType(self):
        return self.__type
    type = property(getType, None, None, "type of the data in the entry")


    def getRank(self):
        return self.__rank
    rank = property(getRank, None, None, "rank of the data in the entry")


    def getShape(self):
        return self.__shape
    shape = property(getShape, None, None, "shape of the data in the entry")


    def getValue(self):
        return self.__value
    value = property(getValue, None, None, "value of the data in the entry")


    def isComparable(self, other):
        """Check if two entries are comparable"""
        return (other.name == self.name and other.type == self.type
                and other.rank == self.rank and other.shape == self.shape)




class TaggedCollection(object):
    """Contains a collection of tagged entries"""

    def __init__(self, entries):
        """file -- open file like object containing collection of tagged data"""
        self.__entryNames = []
        self.__entryLines = []
        self.__entries = []
        self.addEntries(entries)


    def addEntries(self, entries):

        for entry in entries:
            taggedLine = ":".join((entry.name, entry.type, str(entry.rank),
                                   ",".join(map(str, entry.shape))))
            self.__entryNames.append(entry.name)
            self.__entryLines.append(taggedLine)
            self.__entries.append(entry)


    def getMatchingEntries(self, pattern):
        """Returns entries from the collection matching a given pattern
        pattern -- compiled regular expression
        """
        result = []
        for iEntry in range(len(self.__entries)):
            if pattern.match(self.__entryLines[iEntry]):
                result.append(self.__entries[iEntry])

        return result


    def getEntry(self, name):
        """Returns an entry with a given name from the collection
        name -- name of the entry
        """
        try:
            iEntry = self.__entryNames.index(name)
        except ValueError:
            result = None
        else:
            result = self.__entries[iEntry]

        return result


    def delEntry(self, name):
        """Deletes the specified entry from the collection
        name -- name of the entry
        """
        try:
            iEntry = self.__entryNames.index(name)
        except ValueError:
            pass
        else:
            del self.__entries[iEntry]
            del self.__entryNames[iEntry]
            del self.__entryLines[iEntry]



############################################################################
# Parses the file containing the data and returns TaggedEntries
############################################################################

class ResultParser(object):
    """Parser the result files containing tagged data"""

    # Pattern for lines containing the describing tag for following data
    patTagLine = re.compile(r"""(?P<name>[^: ]+)\s*:
                                (?P<type>[^:]+):
                                (?P<rank>\d):
                                (?P<shape>(?:\d+(?:,\d+)*)*)
                            """, re.VERBOSE)


    def __init__(self, file):
        """file -- file like object containing tagged data"""
        self.__file = file


    def iterateEntries(self):
        """Generator for iterating over the entries of the data file."""

        name = None
        type = None
        rank = None
        shape = None
        value = []
        iLine = 0
        for line in self.__file.readlines():
            iLine = iLine + 1
            failed = False
            match = self.patTagLine.match(line)
            if match:

                # Append data from previous tag if present
                if name:
                    try:
                        yield TaggedEntry(name, type, rank, shape,
                                          " ".join(value))
                    except InvalidEntry as ee:
                        raise InvalidEntry(iTaggedLine + 1, iLine, msg=ee.msg)

                name = match.group("name")
                type = match.group("type")
                rank = int(match.group("rank"))
                if rank > 0:
                    shape = tuple([ int(s)
                                    for s in match.group("shape").split(",") ])
                else:
                    shape = ()
                value = []
                iTaggedLine = iLine
            else:
                value.append(line)

        # process last entry
        if name:
            try:
                yield TaggedEntry(name, type, rank, shape, " ".join(value))
            except InvalidEntry as ee:
                raise InvalidEntry(iTaggedLine + 1, iLine, msg=ee.msg)

    entries = property(iterateEntries, None, None,
                       "Iterator over parsed entries")


def product(elements):
    res = 1
    for elem in elements:
        res *= elem
    return res

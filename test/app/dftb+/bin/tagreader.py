#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2026  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

"""Parser for the tagged output of DFTB.

The tagged data is assumed to be represented in the form as provided by the taggedout module of the
DFTB project. See the appropriate source code for the details.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from math import prod
import re
from typing import Any, IO, Iterable, Self
import numpy as np
from numpy.typing import NDArray


#
# Exceptions
#


@dataclass
class InvalidEntry(ValueError):
    """Raised if entry can't be initialized or wrong type of entry is passed

    Attributes:
        start: Starting line of the block where the error occurred.
        end: First line after ending line of the block where the error occurred.
        msg: Reason for the error.
    """

    start: int = 0
    end: int = 0
    msg: str = ""


class ConversionError(ValueError):
    """Raised if error occurs during conversion from string"""


#
# Data converters (from string to Python object)
#


class Converter(ABC):
    """Base class for string to value converters"""

    def __init__(self, nolist: bool = False):
        """Initializes the converter.

        Args:
            nolist: If True, not a list, but only a single value is returned
                after conversion.
        """
        self.nolist = nolist

    def __call__(self, strvalue: str) -> list[Any] | Any:
        """Converts a string representation into Python objects.

        Args:
            strvalue: String representation of the values to convert.

        Returns:
            List of converted values, or a single value if nolist is True.

        Raises:
            ConversionError: If nolist is True but multiple values are provided.
        """
        values = strvalue.split()
        if self.nolist and len(values) > 1:
            raise ConversionError("Too many values")
        result = self.convert(values)
        return result[0] if self.nolist else result

    @abstractmethod
    def convert(self, values: list[str]) -> list[Any]:
        """Conversion function to be implemented by subclasses.

        Args:
            values: List of strings to convert.

        Returns:
            List of converted values.
        """


class FloatConverter(Converter):
    """Converts string to float"""

    def convert(self, values: list[str]) -> list[float]:
        try:
            return [float(val) for val in values]
        except ValueError as err:
            raise ConversionError(f"Float conversion error: {err}") from err


class IntConverter(Converter):
    """Converts string to integer"""

    def convert(self, values: list[str]) -> list[int]:
        try:
            return [int(val) for val in values]
        except ValueError as err:
            raise ConversionError(f"Integer conversion error: {err}") from err


class ComplexConverter(Converter):
    """Converts string to complex."""

    def __call__(self, strvalue: str) -> list[complex] | complex:
        values = strvalue.split()
        if len(values) % 2:
            raise ConversionError("Odd number of values for complex type")
        if self.nolist and len(values) != 2:
            raise ConversionError("Too many values for single complex")
        result = self.convert(values)
        return result[0] if self.nolist else result

    def convert(self, values: list[str]) -> list[complex]:
        ll = []
        for i in range(0, len(values), 2):
            try:
                ll.append(complex(float(values[i]), float(values[i + 1])))
            except (ValueError, IndexError) as err:
                raise ConversionError(f"Unable to convert complex pair: {values[i:i+2]}") from err
        return ll


class LogicalConverter(Converter):
    """Converts string to logical"""

    _mapping = {"t": 1, "f": 0}

    def convert(self, values):
        try:
            return [self._mapping[val.lower()] for val in values]
        except KeyError as err:
            raise ConversionError(f"Logical conversion error {err}") from err


#
# Tagged data related objects
#


class TaggedEntry:
    """Represents a tagged entry with data."""

    _converters = {
        "integer": IntConverter(),
        "real": FloatConverter(),
        "complex": ComplexConverter(),
        "logical": LogicalConverter(),
    }

    _NUMPY_DTYPE_MAP = {
        "integer": np.int64,
        "real": np.float64,
        "complex": np.complex128,
        "logical": np.bool_,
    }


    def __init__(self, name: str, dtype: str, rank: int, shape: tuple[int, ...], strvalue: str):
        """Instantiates a TaggedEntry object.

        Args:
            name: Name of the tagged entry.
            dtype: Type of the data (e.g., 'real', 'integer').
            rank: Rank of the data.
            shape: Shape of the data.
            strvalue: Raw string representation.

        Raises:
            InvalidEntry: If data or metadata is inconsistent.
        """
        if dtype not in self._converters:
            raise InvalidEntry(msg=f"Invalid data type '{dtype}'")

        self._name = name
        self._dtype = dtype
        self._rank = rank
        self._shape = shape

        try:
            self._value = self._converters[dtype](strvalue)
        except ConversionError as err:
            raise InvalidEntry(msg=str(err)) from err

        if (shape and len(shape) != rank) or (not shape and rank != 0):
            raise InvalidEntry(msg="Incompatible rank and shape")

        if shape and (len(self._value) != prod(shape)):
            raise InvalidEntry(msg="Invalid number of values for specified shape")

    @property
    def name(self) -> str:
        """Name of the tagged entry"""
        return self._name

    @property
    def dtype(self) -> str:
        """Data type of the tagged entry"""
        return self._dtype

    @property
    def rank(self) -> int:
        """Rank of the data in the tagged entry"""
        return self._rank

    @property
    def shape(self) -> tuple[int, ...]:
        """Shape of the data in the tagged entry

        The shape is the one used in Fortran (using the column major convention).
        """
        return self._shape

    @property
    def value(self) -> list[Any]:
        """Flat list containing the values of the data in the tagged entry"""
        return self._value

    @property
    def value_array(self) -> NDArray[Any]:
        """Value as a numpy array.

        The Fortran shape is reversed in the returned array, but the data is not transposed.  That
        means, the data is aligned in memory exaclty in the same order as in Fortran, but you
        traverse the consecutive elements in memory with the Python (row-major) indexing convention.

        """
        arrayshape = self._shape[::-1] if self._shape else None
        return np.array(self._value, dtype=self._NUMPY_DTYPE_MAP[self._dtype]).reshape(arrayshape)

    def is_comparable(self, other: Self) -> bool:
        """Check if metadata matches."""
        return (
            isinstance(other, TaggedEntry)
            and self.name == other.name
            and self.dtype == other.dtype
            and self.rank == other.rank
            and self.shape == other.shape
        )


class TaggedCollection:
    """Contains a collection of tagged entries"""

    def __init__(self, entries: Iterable[TaggedEntry]):
        """Initializes a collection of tagged items.

        Args:
            entries: List of entries to store.
        """
        self._entries: dict[str, TaggedEntry] = {}
        self._header_cache: dict[str, str] = {}
        self.add_entries(entries)

    def add_entries(self, entries: Iterable[TaggedEntry]) -> None:
        """Adds new entries to the collection of tagged items.

        Args:
            entries: List of entries to add.
        """
        for entry in entries:
            self._entries[entry.name] = entry
            header = ":".join(
                (entry.name, entry.dtype, str(entry.rank), ",".join(map(str, entry.shape)))
            )
            self._header_cache[entry.name] = header

    def get_matching_entries(self, pattern: re.Pattern) -> list[TaggedEntry]:
        """Returns entries where the serialized tag matches the pattern.

        Args:
            pattern: Compiled regular expression pattern.

        Returns:
            List of the matching entries.
        """
        results = []
        for name, header in self._header_cache.items():
            if pattern.match(header):
                results.append(self._entries[name])
        return results

    def get_entry(self, name: str) -> TaggedEntry | None:
        """Returns an entry with a given name from the collection.

        Args:
            name: Name of the entry to look for.

        Returns
            The entry with the given name or None if not found.
        """
        return self._entries.get(name)

    def del_entry(self, name: str) -> None:
        """Removes an entriy with the given name.

        Args:
            name: Name of the entry to remove
        """
        self._entries.pop(name, None)
        self._header_cache.pop(name, None)


class TaggedDataParser:
    """Parser for tagged data"""

    # Pattern for tagged headers
    TAG_PATTERN = re.compile(
        r"^(?P<name>[^: ]+)\s*:(?P<type>[^:]+):(?P<rank>\d):(?P<shape>(?:\d+(?:,\d+)*)*)$"
    )

    def __init__(self, fobj: IO[str]):
        """Initializes with an open text file."""
        self._fobj = fobj

    @property
    def entries(self):
        """Generator for iterating over the entries of the data file."""

        name, dtype, rank, shape = None, None, None, None
        value_buffer = []
        tag_line_num = 0
        line_num = 0
        for line_num, line in enumerate(self._fobj, start=1):
            match = self.TAG_PATTERN.match(line)
            if match:
                # New tag found -> reading of prev. entry finished: create entry with cached data
                if name:
                    yield self._new_entry(
                        name, dtype, rank, shape, value_buffer, tag_line_num, line_num
                    )
                    name, dtype, rank, shape = None, None, None, None
                    value_buffer = []
                name = match.group("name")
                dtype = match.group("type")
                rank = int(match.group("rank"))
                shape = tuple(int(s) for s in match.group("shape").split(",")) if rank > 0 else ()
                tag_line_num = line_num
            else:
                value_buffer.append(line)

        # Iteration finished: create last entry with cached data
        if name:
            yield self._new_entry(
                name, dtype, rank, shape, value_buffer, tag_line_num, line_num + 1
            )

    def _new_entry(self, name, dtype, rank, shape, buffer, start, end) -> TaggedEntry:
        """Generates a new TaggedEntry or raises InvalidEntry exception if not possible."""
        try:
            return TaggedEntry(name, dtype, rank, shape, " ".join(buffer))
        except InvalidEntry as err:
            raise InvalidEntry(start, end, msg=err.msg) from err

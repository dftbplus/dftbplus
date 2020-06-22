#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

############################################################################
#
# Simple iterator for reading non-empty and non-comment lines from
# a file like object
#
############################################################################

class FreeSepString(str):
    """Free Separated String
    New string class: extends the built-in's strip() for arbitary
    separators not only whitespace chars."""


    def strip(self, strSeparator = None):
        """Return itself without the leading or trailing separators.
        
        Arguments:
            strSeparator -- separator string (default: None -- whitespace
                                            characters)

        Exceptions:
            TypeError -- if strSeparator is not string instance.
        """
        if strSeparator:
            if not (isinstance(strSeparator, str)):
                raise TypeError("expected a character buffer object")
            splitted = self.split(strSeparator)
            empties = [ len(x) == 0 for x in splitted ]
            try:
                i1 = empties.index(0)
            except ValueError:
                return ""
            rev =empties[:]
            rev.reverse()
            i2 = len(empties) - rev.index(0)
            return strSeparator.join(splitted[i1:i2])
        else:
            return str.strip(self)

                

class UncommLines(object):
    """Simple iterator which returnes the nonempty lines in a file. Comment
    strings and text behind them is ignored. Whitespace chars (\t \n etc.)
    before/after the first/last non-whitespace char in the lines are stripped.
    """
        
    def __init__(self, file, iStartLine = 0, iEndLine = 0, strComment = "#",
                             strSeparator = None, returnLineNr = None):
        """Constructor for UncommLines

        Arguments:
                file         -- a file like object with a readline() method.
                iStartLine   -- first nonempty line to return (default: 0)
                iEndLine     -- first nonempty line not to return
                    (default: 0 -- all lines are returned)
                strComment   -- comment string
                strSeparator -- separator string (default: None -- whitespace
                        chars like in the split() method of strings)
                returnLineNr  -- if line nr. of current line should be returned

        Exceptions:
                AttributeError -- if file has no readline() method.
        """

        if hasattr(file, 'readline'):
            self.__file = file
        else:
            raise AttributeError("Specified file has no readline() method")
        try:
            self.__iStartLine = iStartLine
            self.__iEndLine = iEndLine
        except (ValueError, TypeError):
            raise TypeError("Range indexes has to be of int or long type")
        self.__strComment = strComment
        self.__strSeparator = strSeparator
        self.__iLine = 0
        self.__iReadLine = 0
        self.__returnLineNr = returnLineNr


    def __iter__(self):
        return self


    def next(self):
        """Returns next item.

        Exceptions:
                StopIteration    -- if there is no next item
        """

        while (not self.__iEndLine) or (self.__iLine < self.__iEndLine):
            try:
                line = self.__file.readline()
                self.__iReadLine += 1
            except IOError:
                pass
            else:
                if line:
                    tmp = (line.split(self.__strComment, 1)[0]).split("\n")[0]
                    tmp2 = FreeSepString(tmp).strip(self.__strSeparator)
                    if len(tmp2):
                        self.__iLine += 1
                        if self.__iLine > self.__iStartLine:
                            if self.__returnLineNr:
                                return (tmp2, self.__iReadLine - 1)
                            else:
                                return tmp2
                else:
                    break
        raise StopIteration
                        

    def __next__(self):
        return self.next()

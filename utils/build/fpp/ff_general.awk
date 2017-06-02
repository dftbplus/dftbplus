#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2017  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

#
# Folds long lines caused by cpp macro substitution in F90 source code.
#
# Puts #line directives before each newly created continuation line if
# its first char is not part of a quotation.
#
#
BEGIN {
  maxLen = 130;              # max. len of lines (without cont. chars)
  FS = "!";                  # comment char
  maxContLines = 39;         # max. nr. of continuation lines
  contChar ="&";             # continuation character
  patStmtSep = ";[ \t]*$";   # pattern for last nonblank char in the line
                             # being a statement separator
}


# Determines the quoted sections in a given string.
# @param line   IN  Line to investigate.
# @param qBegin OUT List of the beginnnig positions (first pos. after the
#                   opening quote).
# @param qEnd   OUT List of the ending positions (first pos. after the closing
#                   quote).
# @return           Nr. of quoted sections found.
# @note The subroutine is not aware of quote escaping, like \" or \'.
#
function getQuotedSections(line, qBegin, qEnd,

			   ll1, iSect, i1, i2, iStart)
{
  ll1 = length(line) + 1;
  iSect = 0;
  i1 = index(line, "\"");
  i2 = index(line, "'");
  iStart = 1;
  while (i1 != 0 || i2 != 0) {
    #  Opening quote
    if ((i2 == 0) || (i1 != 0 && i1 < i2)) {
      i2 = 0;
      iSect++;
      iStart += i1;                                 # iStart += (i1 - 1) + 1
      qBegin[iSect] = iStart;
      i1 = index(substr(line, iStart), "\"");
      if (i1 != 0) {
	iStart += i1;
	qEnd[iSect] = iStart;
	i1 = index(substr(line, iStart), "\"");
	i2 = index(substr(line, iStart), "'");
      }
      else {
	qEnd[iSect] = ll1;
      }
    }
    # Opening apostrophe
    else {
      i1 = 0;
      iSect++;
      iStart += i2;
      qBegin[iSect] = iStart;
      i2 = index(substr(line, iStart), "'");
      if (i2 != 0) {
	iStart += i2;
	qEnd[iSect] = iStart;
	i1 = index(substr(line, iStart), "\"");
	i2 = index(substr(line, iStart), "'");
      }
      else {
	qEnd[iSect] = ll1;
      }
    }
  }
  return iSect;
}



# Returns the end position of a quoted section containing a given position.
# @param iPos     IN Position to investigate.
# @param nrQSects IN Nr. of quited sections.
# @param qBegin   IN List of beginning positions.
# @param qEnd     IN List of ending positions.
# @return            Position of the next closing quote if iPos falls in a
#                    quoted section; othewise 0.
#
function quoteEnd(iPos, nrQSects, qBegin, qEnd,

		  ii)
{
  if (!nrQSects)
    return 0;
  for (ii = 0; ii < nrQSects && iPos >= qBegin[ii+1]; ii++) ;
  if (ii && iPos < qEnd[ii]) {
    return qEnd[ii];
  }
  else {
    return 0;
  }
}



# Splits a string in pieces with a given maximal length.
# @param text     IN  String to process.
# @param maxLen   IN  Maximal length of one piece.
# @param nrQSects IN  Nr. of quoted sections in the string.
# @param qBegin   IN  List of beginning positions of quotedd sections.
# @param qEnd     IN  List of ending positions of quoted sections.
# @param splitted OUT List of the splitted pieces.
# @param quited   OUT List of boolens (one for each splitted piece). True,
#                     it the first char. of the piece falls in a quited section,
#                     false otherwise.
# @return             Nr. of pieces.
# @note The function tries to minimize the nr. of pieces with quoted first char.
#
function splitString(text, maxLen, nrQSects, qBegin, qEnd, splitted, quoted,

		     ii, nrCuts, splitPos, quoteEndPos, lenText) 
{
  lenText = length(text);
  nrCuts = 1;
  splitPos[nrCuts] = 1;
  quoted[nrCuts] = 0;
  quoteEndPos = 0;
  ii = 1 + maxLen;
  while (ii <= lenText || (quoteEndPos > 0 && quoteEndPos <= lenText)) {
    if (!quoteEndPos) {
      quoteEndPos = quoteEnd(ii, nrQSects, qBegin, qEnd);
    }
    if (quoteEndPos && quoteEndPos <= ii) {
      nrCuts++;
      splitPos[nrCuts] = quoteEndPos;
      quoted[nrCuts] = 0;
      ii = quoteEndPos + maxLen;
      quoteEndPos = 0;
    }
    else {
      nrCuts++;
      splitPos[nrCuts] = ii;
      quoted[nrCuts] = (quoteEndPos != 0);
      ii += maxLen;
    }
  }
  splitPos[nrCuts + 1] = lenText + 1;
  for (ii = 1; ii <= nrCuts; ii++) {
    splitted[ii] = substr(text, splitPos[ii], splitPos[ii+1] - splitPos[ii]);
  }
  return nrCuts;
}



# Returns the continuation character for a line
# @param  splitted  Array containing the lines.
# @param  quoted    Array containing flag for quoted beginning for each line
# @param  iLine     Nr. of line to analyze.
# @return      Global continuation character unless the line ends with an
#              unquoted statement separator.
function getContChar(splitted, quoted, iLine)
{
  if (!quoted[iLine + 1] && match(splitted[iLine], patStmtSep)) {
    return "";
  }
  else {
    return contChar;
  }
}



# Line directives -> update shift between virtual and real linenrs.
#
/\# *[0-9]*/ {
  nrFields = split($0, A, " +");
  lineNr = A[2];
  fileName = A[3];
  #invocNr = A[4];
  realLineNr = NR + 1;
  print;
  next;
}



# General line
#
{ 
  # remove trailing space
  longline = $0;
  sub("[ \t]+$", "", longline);

  # if line is too long or contains "& &" strings -> process it
  if (length(longline) > maxLen || match(longline, "&[ \t]*&")) {

    # Remove leading space
    sub("^[ \t]+", "", longline);

    # Split line in statement and comment part (take care of quoted "!")
    nrQSects = getQuotedSections(longline, qBegin, qEnd);
    iComment = 0;
    iStart = 0;
    iExc = index(longline, "!");
    while (iExc) {
      iExc += iStart
      if (!quoteEnd(iExc, nrQSects, qBegin, qEnd)) {
	iComment = iExc;
	break;
      }
      iStart = iExc;
      iExc = index(substr(longline, iStart + 1), "!");
    }
    if (iComment) {
      line = substr(longline, 1, iComment - 1);
      #comment = substr(longline, iComment);
    }
    else {
      line = longline;
      #comment = "";
    }
    sub("[ \t]+$", "", line);

    # Remove "&  &" strings from non-comment part of the line
    iCont = match(line, "&[ \t]*&");
    unchanged = 1;
    if (iCont) {
      iContReal = 1;
      newline = "";
      while (iCont) {
	if (quoteEnd(iCont, nrQSects, qBegin, qEnd)) {
	  newline = sprintf("%s%s", newline, 
			    substr(line, iContReal, iCont - 1 + RLENGTH));
	}
	else {
	  newline = sprintf("%s%s", newline, 
			    substr(line, iContReal, iCont - 1));
	  unchanged = 0;
	}
	iContReal += RLENGTH + iCont - 1;
	iCont = match(substr(line, iContReal), "&[ \t]*&");
      }
      line = sprintf("%s%s", newline, substr(line, iContReal));
    }

    # If "& &" was a comment and line does not need folding -> print line
    if (unchanged && length(longline) <= maxLen) {
      print longline
      next;
    }

    # Statement part longer than maxlen -> fold it
    if (length(line) > maxLen) {
      lineNr = NR - realLineNr + lineNr;
      realLineNr = NR;
      nrLines = splitString(line, maxLen, nrQSects, qBegin, qEnd,
			    splitted, quoted);
      if (nrLines > maxContLines) {
	printf("!!! FOLDING ERROR:\nNr. of folded lines > 39!\n");
	exit 1;
      }
      oldContChar = getContChar(splitted, quoted, 1);
      printf("%s%s\n", splitted[1], oldContChar);
      for (i = 2; i < nrLines; i++) {
	if (!quoted[i]) {
	  printf("# %d %s\n", lineNr, fileName);
	}
	curContChar = getContChar(splitted, quoted, i);
	printf("%s%s%s\n", oldContChar, splitted[i], curContChar);
	oldContChar = curContChar;
      }
      if (!quoted[nrLines]) {
	printf("# %d %s\n", lineNr, fileName);
      }
      printf("%s%s\n", oldContChar, splitted[nrLines]);
      printf("# %d %s\n", lineNr + 1, fileName);
    }
    # Statement part shorter than maxlen -> print it
    else {
      print line;
    }
  }

  # Line is not too long -> print
  else {
    print longline;
  }
}

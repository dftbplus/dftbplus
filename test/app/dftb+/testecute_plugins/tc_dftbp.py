#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
"""Implements testecute actions for dftbplus."""

import glob
import re
from testecute import registry, StepResult, StepStatus


@registry.register_action("dftbp.assert_absent")
def assert_absent(
    patterns: list[str],
    globs: list[str],
    require_glob_match: bool = True,
    flags: list[str] | None = None,
):
    """Checks files for the occurance of substrings.

    Args:
        patterns: The patterns (regular expressions) to search for.
        globs: List of file globs to search in.
        require_glob_match: If True, at least one file must match each glob pattern.
        flags: List of regex flags to use (e.g., 'IGNORECASE', 'MULTILINE').
    """
    flags = flags or []
    re_flags = 0
    for flag in flags:
        try:
            re_flags |= re.RegexFlag[flag.upper()]
        except KeyError:
            return StepResult(StepStatus.ERROR, message=f"Invalid regex flag: '{flag}'")
    regexes = [re.compile(pattern, re_flags) for pattern in patterns]

    files = []
    for pattern in globs:
        matched_files = glob.glob(pattern)
        if require_glob_match and not matched_files:
            return StepResult(
                StepStatus.ERROR, message=f"No files matched the glob pattern '{pattern}'."
            )
        files.extend(matched_files)

    clean_files: list[str] = []
    affected_files: list[tuple[str, list[str]]] = []
    for file in files:
        matching_patterns = []
        with open(file, "r", encoding="utf-8") as f:
            content = f.read()
        for regex in regexes:
            if regex.search(content):
                matching_patterns.append(regex.pattern)
        if matching_patterns:
            affected_files.append((file, matching_patterns))
        else:
            clean_files.append(file)

    details = {"affected_files": affected_files, "clean_files": clean_files}
    if affected_files:
        return StepResult(
            StepStatus.FAILED, "Some of the files contains matching content", details=details
        )
    return StepResult(StepStatus.OK, details=details)

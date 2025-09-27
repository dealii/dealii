#!/usr/bin/env python3

## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------
"""Search the ./include directory and print out a summary of each deprecated
object, class, function, etc. Compatile with emacs' compilation-mode.

You may need to run

    git fetch --tags upstream

to make the tag (i.e., which release corresponds to which commit) printing work.
"""

from __future__ import print_function

import datetime
import errno
import os
import re
import subprocess
import sys


def decode_and_split(string):
    """Decode a string as unicode and split its lines. Works with python 2 and
    python 3.
    """
    if sys.version_info.major > 2:
        return str(string, "utf8").splitlines()
    else:
        return unicode(string, errors="replace").splitlines()


class DeprecatedDeclaration(object):
    """A record of a deprecated declaration. The input line containing the
    deprecation usually comes from some grep-like program and is assumed to have
    a file name, line number, and contain DEAL_II_DEPRECATED, e.g.,

    ./include/deal.II/base/utilities.h:457: bool job_supports_mpi() DEAL_II_DEPRECATED;

    is one such valid line. The file name and line number are used to query git
    for information regarding when the deprecation notice appeared. The line
    itself is printed out later.
    """

    def __init__(self, line):
        assert "DEAL_II_DEPRECATED" in line
        # in case the line contains extra ':'s, split input based on where grep
        # puts the line number:
        line_number_re = re.search(":[0-9]+:", line)
        if line_number_re:
            line_number_range = (
                line_number_re.regs[0][0] + 1,
                line_number_re.regs[0][1] - 1,
            )
            self.file_name = line[0 : line_number_range[0] - 1]
            self.line_n = int(line[slice(*line_number_range)])
            self.deprecation_line = line[line_number_range[1] + 1 :].strip()
        else:
            raise ValueError(
                "The given line does not contain a line number of "
                "the form, e.g., :42:."
            )

        git_log_output = subprocess.check_output(
            ["git", "blame", "-p", "-L", "{0},{0}".format(self.line_n), self.file_name]
        )

        self.git_log_output = decode_and_split(git_log_output)
        self.commit_hash = self.git_log_output[0].split(" ")[0]
        self.commit_summary = self.git_log_output[9][len("summary ") :]
        self.epoch_time = int(self.git_log_output[7][len("commiter-time ") :])
        self.output_time = datetime.datetime.fromtimestamp(
            self.epoch_time, datetime.UTC
        )

        git_tag_output = subprocess.check_output(
            ["git", "tag", "--contains", self.commit_hash, "-l", "v[0-9].[0-9].[0-9]"]
        )
        git_tag_output = decode_and_split(git_tag_output)
        if git_tag_output:
            # matched tags must start with 'v[0-9]': skip the v
            self.release = git_tag_output[0][1:]
        else:
            self.release = ""

    def print_record(self):
        """Print a description of the record to stdout. The 'Release' field
        corresponds to the release number in which the declaration was first
        marked as deprecated.
        """
        print("File: " + self.file_name + ":" + str(self.line_n) + ": ")
        print("Line: " + self.deprecation_line)
        print("Date: " + str(self.output_time))
        print("Message: " + self.commit_summary)
        print("Hash: " + self.commit_hash)
        print("Release: " + self.release)


def main():
    # check that we are in the top directory
    if not all(
        (
            os.path.isdir(directory)
            for directory in ["./include", "./tests", "./examples"]
        )
    ):
        raise Exception(
            "This script must be called from the top-level directory of deal.II."
        )
    result = subprocess.check_output(
        ["grep", "-R", "--line-number", "DEAL_II_DEPRECATED", "./include/"]
    )
    result = decode_and_split(result)

    deprecated_declarations = list()
    for line in result:
        try:
            deprecated_declarations.append(DeprecatedDeclaration(line))
        except subprocess.CalledProcessError:  # ignore errors coming from git
            pass

    deprecated_declarations.sort(key=lambda u: u.epoch_time)

    try:
        for record in deprecated_declarations:
            record.print_record()
            print("")

        sys.stderr.flush()
        sys.stdout.flush()
    except IOError as pipe_error:
        if pipe_error.errno == errno.EPIPE:
            # this is OK, it usually just means the user piped this script to
            # head or some other screen capturing program that ignores some of
            # its input
            pass


if __name__ == "__main__":
    main()

#!/usr/bin/env python
## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2019 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

# This script checks if the files in the deal.II repository are encoded with
# valid UTF8.
#
# This script should be invoked from the root folder of the deal.II
# repository:
# contrib/utilities/check_encoding.py
from __future__ import print_function

import itertools
import io
import os
import sys


def filename_generator(suffix):
    for root, _, file_names in os.walk("./"):
        for file_name in file_names:
            if file_name.endswith(suffix):
                if root == "./":
                    yield root + file_name
                else:
                    yield root + "/" + file_name


filenames = itertools.chain(
    filename_generator(".h"), filename_generator(".cc"), filename_generator(".html")
)

return_code = 0
for filename in filenames:
    file_handle = io.open(filename, encoding="utf-8")
    try:
        file_handle.read()
    except UnicodeDecodeError:
        print(filename + " is not encoded with UTF-8")
        return_code = 1
    finally:
        file_handle.close()

sys.exit(return_code)

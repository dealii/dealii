#!/usr/bin/env python
## ---------------------------------------------------------------------
##
## Copyright (C) 2019 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

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


filenames = itertools.chain(filename_generator(".h"),
                            filename_generator(".cc"),
                            filename_generator(".html"))

return_code = 0
for filename in filenames:
    file_handle = io.open(filename, encoding='utf-8')
    try:
        file_handle.read()
    except UnicodeDecodeError:
        print(filename + ' is not encoded with UTF-8')
        return_code = 1
    finally:
        file_handle.close()

sys.exit(return_code)

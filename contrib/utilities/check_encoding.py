#!/usr/bin/env python3
## ---------------------------------------------------------------------
##
## Copyright (C) 2019 by the deal.II authors
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
# This script should be invoked from the root folder of the deal.ii
# repository:
# contrib/utilities/check_encoding.py

import glob
import itertools

filenames = itertools.chain(glob.iglob('**/*.h', recursive=True),
                            glob.iglob('**/*.cc', recursive=True),
                            glob.iglob('**/*.html', recursive=True))

for filename in filenames:
    try:
        with open(filename, encoding='utf-8') as file:
            file.read()
    except:
        raise Exception(filename + ' is not encoded is not encoded with UTF-8')

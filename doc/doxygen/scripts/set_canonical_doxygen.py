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

# This script sets the canonical webpage for all the webpages of the
# doxygen documentation. For more information see
# https://en.wikipedia.org/wiki/Canonical_link_element
#
# This script is invoked by CMake. To add the canonical link to the webpages of
# the deal.ii repository you can use the script
# contrib/utilities/set_canonical_webpages.py

import glob

filenames = glob.iglob('**/*html', recursive=True)

for filename in filenames:
    file = open(filename, 'r')
    new_text_data = str()
    if '<link rel="canonical"' not in file.read():
        file.seek(0)
        for line in file.readlines():
            # Do not add the canonical link twice
            canonical_link_added = False
            new_text_data += line
            if (not canonical_link_added) and ('<head>' in line):
                new_text_data += ('<link rel="canonical" href="https://www.dealii.org/current/doxygen/deal.II/'
                                  + filename + '" />' + '\n')
                canonical_link_added = True
        file.close()

        # Truncate the file and write the new text
        file = open(filename, 'w+')
        file.write(new_text_data)
        file.close()

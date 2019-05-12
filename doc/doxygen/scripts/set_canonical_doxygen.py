#!/usr/bin/env python
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

import os

LINK_START = ('<link rel="canonical" href="https://www.dealii.org/current'
              + '/doxygen/')

def filename_generator():
    for root, _, file_names in os.walk("./"):
        for file_name in file_names:
            if file_name.endswith(".html"):
                if root == "./":
                    yield root + file_name
                else:
                    yield root + "/" + file_name

for filename in filename_generator():
    # there is no relevant content in the header
    if filename == "./header.html":
        pass
    new_text_data = str()
    with open(filename, 'r') as file_handle:
        if '<link rel="canonical"' not in file_handle.read():
            file_handle.seek(0)
            for line in file_handle.readlines():
                # Do not add the canonical link twice
                canonical_link_added = False
                new_text_data += line
                if not canonical_link_added and '<head>' in line:
                    assert filename[:2] == "./"
                    new_text_data += LINK_START + filename[2:] + '" />\n'
                    canonical_link_added = True
    if new_text_data:
        # Truncate the file and write the new text
        with open(filename, 'w+') as file_handle:
            file_handle.write(new_text_data)

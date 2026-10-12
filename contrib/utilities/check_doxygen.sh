#!/bin/bash

## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2020 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

if [ ! -f contrib/utilities/checkdoxygen.py ]; then
  echo "*** This script must be run from the top-level directory of deal.II."
  exit 1
fi


find doc examples include \( -name "*.h" -o -name "*.dox" \) -print0 | xargs -0 -n 1 contrib/utilities/checkdoxygen.py

#!/bin/sh
## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

# 
# This is a little script that checks if the git history of the 
# feature branch is linear, 
# i.e. no merges of branch 'master' into that branch are present.
# 
# NOTE: This script will do nothing in branch 'master' since the log will always be empty
#       

if [[ $(git log master.. --merges) ]]; then
  echo "There are merge commits in this feature branch!"
  exit 1
fi

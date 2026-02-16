## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# Configuration for the ARPACK library:
#

set(FEATURE_ARPACK_DEPENDS LAPACK)

configure_feature(ARPACK)
set(DEAL_II_ARPACK_WITH_PARPACK ${ARPACK_WITH_PARPACK})

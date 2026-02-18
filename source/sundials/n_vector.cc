// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



#include <deal.II/base/config.h>

#include <deal.II/sundials/n_vector.h>
#include <deal.II/sundials/n_vector.templates.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <sundials/sundials_nvector.h>
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
#    include <sundials/sundials_context.h>
#  endif

DEAL_II_NAMESPACE_OPEN

// We don't build the .inst file if deal.II isn't configured
// with SUNDIALS, but doxygen doesn't know that and tries to find that
// file anyway for parsing -- which then of course it fails on. So
// exclude the following from doxygen consideration.
#  ifndef DOXYGEN
#    include "sundials/n_vector.inst"
#  endif

DEAL_II_NAMESPACE_CLOSE

#endif

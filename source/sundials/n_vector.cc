// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/base/config.h>

#include <deal.II/sundials/n_vector.h>
#include <deal.II/sundials/n_vector.templates.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <sundials/sundials_nvector.h>
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
#    include <sundials/sundials_context.h>
#  endif

DEAL_II_NAMESPACE_OPEN

#  include "sundials/n_vector.inst"

DEAL_II_NAMESPACE_CLOSE

#endif

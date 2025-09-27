// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
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

#ifdef DEAL_II_TRILINOS_WITH_SACADO

#  include <deal.II/differentiation/ad/sacado_number_types.h>

DEAL_II_NAMESPACE_OPEN

/*---------------------- Explicit Instantiations ----------------------*/

#  include "differentiation/ad/sacado_number_types.inst1"
#  ifdef DEAL_II_TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD
#    include "differentiation/ad/sacado_number_types.inst2"
#  endif

DEAL_II_NAMESPACE_CLOSE

#endif

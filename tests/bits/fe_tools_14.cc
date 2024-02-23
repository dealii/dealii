// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include "../tests.h"

#include "fe_tools_common.h"

// check
//   FE::hp_constraints_are_implemented ()
// a bit like hp_constraints_are_implemented, but with a different
// set of elements



template <int dim>
void
check_this(const FiniteElement<dim> &fe1, const FiniteElement<dim> &fe2)
{
  deallog << (fe1.hp_constraints_are_implemented() ? "true" : "false")
          << std::endl;
  deallog << (fe2.hp_constraints_are_implemented() ? "true" : "false")
          << std::endl;
}

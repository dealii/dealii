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


#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "dof_tools_common.h"

// check
//   FE::hp_constraints_are_implemented
// a bit like fe_tools_14, but works on a different set of elements



template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler)
{
  deallog << dof_handler.get_fe().get_name() << ": "
          << (dof_handler.get_fe().hp_constraints_are_implemented() ? "true" :
                                                                      "false")
          << std::endl;
}

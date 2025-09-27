// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2022 by the deal.II authors
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
//   DoFTools::count_boundary_dofs



template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler)
{
  // no other args
  deallog << dof_handler.n_boundary_dofs() << std::endl;

  // with std::map
  std::map<types::boundary_id, const Function<dim> *> fm;
  fm[0] = nullptr;
  deallog << dof_handler.n_boundary_dofs(fm) << std::endl;

  // with std::set
  const std::set<types::boundary_id> s = {0};
  deallog << dof_handler.n_boundary_dofs(s) << std::endl;
}

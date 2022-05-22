// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


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

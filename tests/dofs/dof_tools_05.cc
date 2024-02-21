// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2023 by the deal.II authors
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

#include "dof_tools_common.h"

// check
//   DoFTools::extract_boundary_dofs



template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler)
{
  ComponentMask component_select(dof_handler.get_fe().n_components(), true);
  IndexSet      boundary_dofs(dof_handler.n_dofs());

  // first with all components
  {
    boundary_dofs =
      DoFTools::extract_boundary_dofs(dof_handler, component_select);
    output_bool_vector(boundary_dofs);
  }

  // next with only every second
  // component
  for (unsigned int i = 1; i < component_select.size(); i += 2)
    component_select.set(i, false);
  {
    boundary_dofs =
      DoFTools::extract_boundary_dofs(dof_handler, component_select);
    output_bool_vector(boundary_dofs);
  }

  // third further restrict to
  // boundary indicator 0
  {
    const std::set<types::boundary_id> boundary_ids = {0};
    boundary_dofs = DoFTools::extract_boundary_dofs(dof_handler,
                                                    component_select,
                                                    boundary_ids);
    output_bool_vector(boundary_dofs);
  }
}

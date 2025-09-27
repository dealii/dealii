// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2020 by the deal.II authors
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
//   DoFTools::extract_hanging_node_constraints



template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler)
{
  const types::global_dof_index n_dofs = dof_handler.n_dofs();

  const IndexSet hanging_node_dofs =
    DoFTools::extract_hanging_node_dofs(dof_handler);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  for (types::global_dof_index dof = 0; dof < n_dofs; ++dof)
    if (hanging_node_dofs.is_element(dof))
      AssertThrow(constraints.is_constrained(dof), ExcInternalError());

  AssertThrow(hanging_node_dofs.n_elements() == constraints.n_constraints(),
              ExcInternalError());
  output_bool_vector(hanging_node_dofs);
}

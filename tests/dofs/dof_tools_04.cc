// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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


#include "../tests.h"

#include "dof_tools_common.h"
#include "dof_tools_common_fake_hp.h"

// check
//   DoFTools::extract_hanging_node_constraints



template <typename DoFHandlerType>
void
check_this(const DoFHandlerType &dof_handler)
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

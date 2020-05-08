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
#include "dof_tools_common_parallel.h"

// check
//  DoFTools::extract_dofs as in dof_tools_12 in parallel for fewer elements
// The output of this file is for every process the (boolean) vector of
// extracted locally owned DoFs (corresponding to local DoF indices), not the
// extracted vector of the global DoF, which has to be mapped with the locally
// owned IndexSet first.

template <typename DoFHandlerType>
void
check_this(const DoFHandlerType &dof_handler)
{
  std::vector<bool> mask(dof_handler.get_fe_collection().n_components(), false);

  // only select first component
  mask[0] = true;
  output_bool_vector(DoFTools::extract_dofs(dof_handler, ComponentMask(mask)));

  // also select last component
  mask.back() = true;
  output_bool_vector(DoFTools::extract_dofs(dof_handler, ComponentMask(mask)));
}

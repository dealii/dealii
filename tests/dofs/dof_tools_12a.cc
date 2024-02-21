// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
#include "dof_tools_common_parallel.h"

// check
//  DoFTools::extract_dofs as in dof_tools_12 in parallel for fewer elements
// The output of this file is for every process the (boolean) vector of
// extracted locally owned DoFs (corresponding to local DoF indices), not the
// extracted vector of the global DoF, which has to be mapped with the locally
// owned IndexSet first.

template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler)
{
  std::vector<bool> mask(dof_handler.get_fe_collection().n_components(), false);

  // only select first component
  mask[0] = true;
  output_bool_vector(DoFTools::extract_dofs(dof_handler, ComponentMask(mask)));

  // also select last component
  mask.back() = true;
  output_bool_vector(DoFTools::extract_dofs(dof_handler, ComponentMask(mask)));
}

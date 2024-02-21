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
//   DoFTools::extract_dofs



template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler)
{
  std::vector<bool> mask(dof_handler.get_fe(0).n_components(), false);

  // only select first component
  mask[0] = true;
  output_bool_vector(DoFTools::extract_dofs(dof_handler, ComponentMask(mask)));

  // also select last component
  mask.back() = true;
  output_bool_vector(DoFTools::extract_dofs(dof_handler, ComponentMask(mask)));
}

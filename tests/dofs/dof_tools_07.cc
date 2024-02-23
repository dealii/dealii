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


#include "../tests.h"

#include "dof_tools_common.h"

// check
//   DoFTools::count_dofs_per_fe_component



template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler)
{
  const std::vector<types::global_dof_index> n_dofs =
    DoFTools::count_dofs_per_fe_component(dof_handler);
  for (unsigned int i = 0; i < n_dofs.size(); ++i)
    deallog << n_dofs[i] << ' ';
  deallog << std::endl;
}

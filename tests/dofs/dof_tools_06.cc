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
//   DoFTools::extract_subdomain_dofs



template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler)
{
  std::vector<bool> dofs(dof_handler.n_dofs());

  for (unsigned int level = 0;
       level < dof_handler.get_triangulation().n_levels();
       ++level)
    {
      DoFTools::extract_subdomain_dofs(dof_handler, level, dofs);
      output_bool_vector(dofs);
    }
}

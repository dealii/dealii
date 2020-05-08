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
//   DoFTools::count_dofs_per_fe_component



template <typename DoFHandlerType>
void
check_this(const DoFHandlerType &dof_handler)
{
  const std::vector<types::global_dof_index> n_dofs =
    DoFTools::count_dofs_per_fe_component(dof_handler);
  for (unsigned int i = 0; i < n_dofs.size(); ++i)
    deallog << n_dofs[i] << " ";
  deallog << std::endl;
}

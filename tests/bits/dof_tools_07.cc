// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include "../tests.h"
#include "dof_tools_common.h"

// check
//   DoFTools::count_dofs_per_component


std::string output_file_name = "output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  std::vector<types::global_dof_index> n_dofs(dof_handler.get_fe().n_components());
  DoFTools::count_dofs_per_component (dof_handler,
                                      n_dofs);
  for (unsigned int i=0; i<n_dofs.size(); ++i)
    deallog << n_dofs[i] << " ";
  deallog << std::endl;
}

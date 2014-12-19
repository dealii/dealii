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
//   DoFTools::extract_boundary_dofs


std::string output_file_name = "output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  std::vector<bool> component_select (dof_handler.get_fe().n_components(),
                                      true);
  std::vector<bool> boundary_dofs (dof_handler.n_dofs());

  // first with all components
  {
    DoFTools::extract_boundary_dofs (dof_handler,
                                     component_select,
                                     boundary_dofs);
    output_bool_vector (boundary_dofs);
  }

  // next with only every second
  // component
  for (unsigned int i=1; i<component_select.size(); i+=2)
    component_select[i] = false;
  {
    DoFTools::extract_boundary_dofs (dof_handler,
                                     component_select,
                                     boundary_dofs);
    output_bool_vector (boundary_dofs);
  }

  // third further restrict to
  // boundary indicator 0
  {
    std::set<types::boundary_id> boundary_ids;
    boundary_ids.insert (0);
    DoFTools::extract_boundary_dofs (dof_handler,
                                     component_select,
                                     boundary_dofs,
                                     boundary_ids);
    output_bool_vector (boundary_dofs);
  }
}

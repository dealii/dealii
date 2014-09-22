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
//   DoFTools::extract_dofs


std::string output_file_name = "output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  std::vector<bool> selected_dofs (dof_handler.n_dofs());
  std::vector<bool> mask (dof_handler.get_fe().n_components(), false);

  // only select first component
  mask[0] = true;
  DoFTools::extract_dofs (dof_handler, ComponentMask(mask), selected_dofs);
  output_bool_vector (selected_dofs);

  // also select last component
  mask.back() = true;
  DoFTools::extract_dofs (dof_handler, ComponentMask(mask), selected_dofs);
  output_bool_vector (selected_dofs);
}

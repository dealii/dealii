//----------------------------  dof_tools_12.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_tools_12.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.cc"

// check
//   DoFTools::extract_dofs


std::string output_file_name = "dof_tools_12.output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  std::vector<bool> selected_dofs (dof_handler.n_dofs());
  std::vector<bool> mask (dof_handler.get_fe().n_components(), false);

                                   // only select first component
  mask[0] = true;
  DoFTools::extract_dofs (dof_handler, mask, selected_dofs);
  output_bool_vector (selected_dofs);

                                   // also select last component
  mask.back() = true;
  DoFTools::extract_dofs (dof_handler, mask, selected_dofs);
  output_bool_vector (selected_dofs);
}

//----------------------------  dof_tools_4.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_tools_4.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.cc"

// check
//   DoFTools::extract_hanging_node_constraints


std::string output_file_name = "dof_tools_04.output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  std::vector<bool> hanging_node_dofs (dof_handler.n_dofs());
  DoFTools::extract_hanging_node_dofs (dof_handler,
                                       hanging_node_dofs);
  output_bool_vector (hanging_node_dofs);
}

//----------------------------  dof_tools_7.cc  ---------------------------
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
//----------------------------  dof_tools_7.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.cc"

// check
//   DoFTools::count_dofs_per_component


std::string output_file_name = "dof_tools_07.output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  std::vector<unsigned int> n_dofs(dof_handler.get_fe().n_components());
  DoFTools::count_dofs_per_component (dof_handler,
                                      n_dofs);
  for (unsigned int i=0; i<n_dofs.size(); ++i)
    deallog << n_dofs[i] << " ";
  deallog << std::endl;
}

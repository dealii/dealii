//----------------------------  dof_tools_6.cc  ---------------------------
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
//----------------------------  dof_tools_6.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.cc"

// check
//   DoFTools::extract_subdomain_dofs


std::string output_file_name = "dof_tools_06.output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  std::vector<bool> dofs(dof_handler.n_dofs());
  
  for (unsigned int level=0; level<dof_handler.get_tria().n_levels(); ++level)
    {
      DoFTools::extract_subdomain_dofs (dof_handler,
                                        level,
                                        dofs);
      output_bool_vector (dofs);
    }
}

//----------------------------  dof_tools_14.cc  ---------------------------
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
//----------------------------  dof_tools_14.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.cc"
#include <lac/vector.h>

// check
//   DoFTools::count_boundary_dofs



std::string output_file_name = "dof_tools_14.output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
                                   // no other args
  deallog << dof_handler.n_boundary_dofs() << std::endl;

                                   // with FunctionMap
  typename FunctionMap<dim>::type fm;
  fm[0] = 0;
  deallog << dof_handler.n_boundary_dofs(fm) << std::endl;
  
                                   // with std::set
  std::set<unsigned char> s;
  s.insert (0);
  deallog << dof_handler.n_boundary_dofs(s) << std::endl;
  
}

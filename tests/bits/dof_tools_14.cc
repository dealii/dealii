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
#include <deal.II/lac/vector.h>

// check
//   DoFTools::count_boundary_dofs



std::string output_file_name = "output";


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
  std::set<types::boundary_id> s;
  s.insert (0);
  deallog << dof_handler.n_boundary_dofs(s) << std::endl;

}

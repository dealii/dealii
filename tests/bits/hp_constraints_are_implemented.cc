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
//   FE::hp_constraints_are_implemented
// a bit like fe_tools_14, but works on a different set of elements



std::string output_file_name = "output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  deallog << dof_handler.get_fe().get_name()
          << ": "
          << (dof_handler.get_fe().hp_constraints_are_implemented() ? "true" : "false")
          << std::endl;
}

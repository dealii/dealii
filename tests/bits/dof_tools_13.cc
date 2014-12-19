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
//   DoFTools::distribute_cell_to_dof_vector



std::string output_file_name = "output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  // this doesn't make much sense if
  // the element is not primitive
  if (dof_handler.get_fe().is_primitive() == false)
    return;

  Vector<double> cell_data (dof_handler.get_tria().n_active_cells());
  for (unsigned int i=0; i<cell_data.size(); ++i)
    cell_data(i) = i;

  // distribute to first component
  Vector<double> dof_data (dof_handler.n_dofs());
  DoFTools::distribute_cell_to_dof_vector (dof_handler,
                                           cell_data,
                                           dof_data);

  // output every third element
  for (unsigned int i=0; i<dof_data.size(); i+=3)
    deallog << dof_data(i) << " ";
  deallog << std::endl;


  // distribute to last component. note that
  // there will still be data left from the
  // first component.
  DoFTools::distribute_cell_to_dof_vector (dof_handler,
                                           cell_data,
                                           dof_data,
                                           dof_handler.get_fe().n_components()-1);
  for (unsigned int i=0; i<dof_data.size(); i+=3)
    deallog << dof_data(i) << " ";
  deallog << std::endl;
}

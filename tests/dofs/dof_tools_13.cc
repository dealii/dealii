// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "dof_tools_common.h"
#include "dof_tools_common_fake_hp.h"

// check
//   DoFTools::distribute_cell_to_dof_vector



template <typename DoFHandlerType>
void
check_this(const DoFHandlerType &dof_handler)
{
  // this doesn't make much sense if
  // the element is not primitive
  if (dof_handler.get_fe().is_primitive() == false)
    return;

  Vector<double> cell_data(dof_handler.get_triangulation().n_active_cells());
  for (unsigned int i = 0; i < cell_data.size(); ++i)
    cell_data(i) = i;

  // distribute to first component
  Vector<double> dof_data(dof_handler.n_dofs());
  DoFTools::distribute_cell_to_dof_vector(dof_handler, cell_data, dof_data);

  // output every third element
  for (unsigned int i = 0; i < dof_data.size(); i += 3)
    deallog << dof_data(i) << " ";
  deallog << std::endl;


  // distribute to last component. note that
  // there will still be data left from the
  // first component.
  DoFTools::distribute_cell_to_dof_vector(
    dof_handler, cell_data, dof_data, dof_handler.get_fe().n_components() - 1);
  for (unsigned int i = 0; i < dof_data.size(); i += 3)
    deallog << dof_data(i) << " ";
  deallog << std::endl;
}

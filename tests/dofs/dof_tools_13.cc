// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "dof_tools_common.h"

// check
//   DoFTools::distribute_cell_to_dof_vector



template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler)
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
    deallog << dof_data(i) << ' ';
  deallog << std::endl;


  // distribute to last component. note that
  // there will still be data left from the
  // first component.
  DoFTools::distribute_cell_to_dof_vector(
    dof_handler, cell_data, dof_data, dof_handler.get_fe().n_components() - 1);
  for (unsigned int i = 0; i < dof_data.size(); i += 3)
    deallog << dof_data(i) << ' ';
  deallog << std::endl;
}

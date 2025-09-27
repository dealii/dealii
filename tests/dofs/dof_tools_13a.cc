// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2022 by the deal.II authors
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

// this is a variant of dof_tools_13. it turned out that the function
// tested here didn't zero out its argument up front, instead adding
// to the previous content. this is not what we intended, so test with
// a preset vector and make sure that we get the same result as in
// dof_tools_13.



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

  IndexSet component_dofs;
  {
    std::vector<bool> component_mask(dof_handler.get_fe().n_components(),
                                     false);
    component_mask[0] = true;
    component_dofs =
      DoFTools::extract_dofs(dof_handler, ComponentMask(component_mask));

    for (unsigned int i = 0; i < dof_data.size(); ++i)
      if (component_dofs.is_element(i) == true)
        dof_data(i) = i + 1;
      else
        dof_data(i) = 0;
  }

  DoFTools::distribute_cell_to_dof_vector(dof_handler, cell_data, dof_data);
  // output every third element
  for (unsigned int i = 0; i < dof_data.size(); i += 3)
    deallog << dof_data(i) << ' ';
  deallog << std::endl;

  // check that no other values were
  // set
  for (unsigned int i = 0; i < dof_data.size(); ++i)
    if (component_dofs.is_element(i) == false)
      AssertThrow(dof_data(i) == 0, ExcInternalError());


  // distribute to last component. by
  // default we distribute to
  // component zero

  // preset the vector again to make
  // sure that the function zeroes out
  // previous content.
  {
    std::vector<bool> component_mask(dof_handler.get_fe().n_components(),
                                     false);
    component_mask.back() = true;
    component_dofs =
      DoFTools::extract_dofs(dof_handler, ComponentMask(component_mask));
    for (unsigned int i = 0; i < dof_data.size(); ++i)
      if (component_dofs.is_element(i) == true)
        dof_data(i) = i + 1;
      else
        dof_data(i) = 0;
  }
  DoFTools::distribute_cell_to_dof_vector(
    dof_handler, cell_data, dof_data, dof_handler.get_fe().n_components() - 1);
  for (unsigned int i = 0; i < dof_data.size(); i += 3)
    deallog << dof_data(i) << ' ';
  deallog << std::endl;

  // check that no other values were
  // set
  for (unsigned int i = 0; i < dof_data.size(); ++i)
    if (component_dofs.is_element(i) == false)
      AssertThrow(dof_data(i) == 0, ExcInternalError());
}

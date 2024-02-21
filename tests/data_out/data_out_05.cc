// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// like the _01 test, except that we declare output component fields
// as vectors if possible. since VTK is the only format that currently
// supports this, use it as the only output format. (deal.II
// intermediate also supports vector-valued components, but we already
// test for this in the _04 tests)
//
// attach the same data twice to see if we can produce output with
// more than one vector field

#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"

#include "data_out_common.h"



template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler,
           const Vector<double>  &v_node,
           const Vector<double>  &v_cell)
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      dof_handler.get_fe().n_components(),
      DataComponentInterpretation::component_is_scalar);
  // if possible, declare the last
  // dim components as vectors
  if (dof_handler.get_fe().n_components() >= dim)
    for (unsigned int i = dof_handler.get_fe().n_components() - dim;
         i < dof_handler.get_fe().n_components();
         ++i)
      data_component_interpretation[i] =
        DataComponentInterpretation::component_is_part_of_vector;

  deallog << "Checking " << dof_handler.get_fe().get_name() << std::endl
          << "Component mask: ";
  for (unsigned int i = 0; i < dof_handler.get_fe().n_components(); ++i)
    deallog << (data_component_interpretation[i] ==
                    DataComponentInterpretation::component_is_scalar ?
                  "scalar " :
                  "vector ");
  deallog << std::endl;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(
    v_node,
    std::vector<std::string>(dof_handler.get_fe().n_components(),
                             "node_data_1"),
    DataOut<dim>::type_dof_data,
    data_component_interpretation);
  data_out.add_data_vector(
    v_node,
    std::vector<std::string>(dof_handler.get_fe().n_components(),
                             "node_data_2"),
    DataOut<dim>::type_dof_data,
    data_component_interpretation);
  data_out.add_data_vector(v_cell, "cell_data", DataOut<dim>::type_cell_data);
  data_out.build_patches();

  data_out.write_vtk(deallog.get_file_stream());
}
